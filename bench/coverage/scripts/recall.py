# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Recall on the held-out runs: of the runs labelled with a layout (or a kit family), what the
binary said. The validation levels (tally.py) measure the other side, the precision of a HIGH call;
a run called below HIGH, ambiguous or not at all counts neither for nor against a level, so this
table shows how often that happens.

Each held-out run with a label independent of the reads (tally.independent) is grouped by that
label: one layout, or a kit family when the label names only the family. Per group:
  runs, studies     held-out runs and the studies they come from
  high_labelled     HIGH, naming the labelled layout (or a layout of the labelled family), or right
                    by adjudication (adjudication.tsv)
  high_other        HIGH, naming another kit: every one is adjudicated disputed (none wrong)
  medium_low        a match at MEDIUM or LOW, of any layout
  ambiguous         the reads could not separate two layouts
  no_match          no layout fits
  refused           the input was refused (not FASTQ, an SRA archive, malformed)
Runs whose label names layouts of no single family are left out; so are read-derived labels.

Writes results/recall.tsv and the README's recall table (between its markers); --check compares
both instead of writing.

Usage: python recall.py <calls.tsv> <adjudication.tsv> <recall.tsv> <README.md> [--check]
"""

import collections
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from invariant import invariant  # noqa: E402
import tally  # noqa: E402

COLUMNS = ["label", "read_kind", "runs", "studies", "high_labelled", "high_other", "medium_low",
           "ambiguous", "no_match", "refused"]
OUTCOMES = COLUMNS[4:]
BEGIN, END = "<!-- recall:begin -->", "<!-- recall:end -->"


def group_of(truth):
    """The label a run is grouped under: its one layout, or the family all its layouts share."""
    if len(truth) == 1:
        return truth[0]
    families = {tally.FAMILY_OF.get(t) for t in truth}
    return f"{families.pop()} (family)" if len(families) == 1 and None not in families else None


def outcome_of(call, adjudications):
    status = call["status"].split(":")[0]
    if status == "refused":
        return "refused"
    invariant(status == "ok", f"{call['run']}: status {call['status']!r} is neither ok nor refused")
    if call["verdict"] == "no_match":
        return "no_match"
    if call["verdict"] == "ambiguous":
        return "ambiguous"
    invariant(call["verdict"] == "match", f"{call['run']}: verdict {call['verdict']!r}")
    if call["confidence"] != "high":
        return "medium_low"
    result, _, _ = tally.outcome(call, adjudications)
    invariant(result != "wrong", f"{call['run']}: a wrong HIGH call (tally.py reports it)")
    return "high_labelled" if result == "right" else "high_other"


def count(calls, adjudications):
    counts = collections.defaultdict(lambda: {"studies": set(), **{o: 0 for o in OUTCOMES}})
    for call in calls:
        truth = tuple(t for t in call["truth"].split(tally.SEPARATOR) if t)
        if not truth or not tally.independent(call["referee"]):
            continue
        group = group_of(truth)
        if group is None:
            continue
        entry = counts[(group, call["read_kind"])]
        entry["studies"].add(call["study"] or f"run:{call['run']}")
        entry[outcome_of(call, adjudications)] += 1
    for (group, kind), e in counts.items():
        invariant(sum(e[o] for o in OUTCOMES) >= len(e["studies"]), f"{group} {kind}: runs < studies")
    return counts


def render_tsv(counts, source):
    out = ["# Recall on held-out runs: what the binary said, by the run's independent label "
           "(bench/coverage/scripts/recall.py).", f"# From: {source}", "\t".join(COLUMNS)]
    for (group, kind), e in sorted(counts.items()):
        runs = sum(e[o] for o in OUTCOMES)
        out.append("\t".join(map(str, [group, kind, runs, len(e["studies"]), *(e[o] for o in OUTCOMES)])))
    return "\n".join(out) + "\n"


def pct(n, d):
    return f"{100 * n / d:.0f}%" if d else "-"


def render_readme(counts):
    """The README table: groups of 20 runs or more, by runs; smaller groups are in recall.tsv."""
    lines = ["| Label | Read | Runs (studies) | HIGH, as labelled | HIGH, another kit (disputed) "
             "| MEDIUM or LOW | Ambiguous | No call |", "|---|---|---|---|---|---|---|---|"]
    shown = [(k, e, sum(e[o] for o in OUTCOMES)) for k, e in counts.items()]
    for (group, kind), e, runs in sorted(shown, key=lambda x: (-x[2], x[0])):
        if runs < 20:
            continue
        no_call = e["no_match"] + e["refused"]
        lines.append(f"| {group} | {kind} | {runs} ({len(e['studies'])}) | {pct(e['high_labelled'], runs)} "
                     f"| {pct(e['high_other'], runs)} | {pct(e['medium_low'], runs)} "
                     f"| {pct(e['ambiguous'], runs)} | {pct(no_call, runs)} |")
    return "\n".join(lines)


def with_table(readme, table):
    invariant(readme.count(BEGIN) == 1 and readme.count(END) == 1, "the README has one recall table")
    head, rest = readme.split(BEGIN)
    _, tail = rest.split(END)
    return f"{head}{BEGIN}\n{table}\n{END}{tail}"


def main(calls_path, adjudication_path, out_path, readme_path, check=False):
    calls = tally.rows(calls_path)
    adjudications = {r["run"]: r for r in tally.rows(adjudication_path)}
    counts = count(calls, adjudications)
    tsv = render_tsv(counts, tally.comment(calls_path).lstrip("# "))
    with open(readme_path) as f:
        readme = f.read()
    want = with_table(readme, render_readme(counts))
    if check:
        with open(out_path) as f:
            invariant(f.read() == tsv, f"{out_path} is stale: rerun recall.py without --check")
        invariant(readme == want, f"{readme_path}'s recall table is stale: rerun recall.py")
        print(f"{out_path} and the README's recall table match")
        return
    with open(out_path, "w") as f:
        f.write(tsv)
    with open(readme_path, "w") as f:
        f.write(want)
    print(render_readme(counts))


if __name__ == "__main__":
    if len(sys.argv) not in (5, 6) or sys.argv[5:] not in ([], ["--check"]):
        sys.exit(__doc__)
    main(*sys.argv[1:5], check=sys.argv[5:] == ["--check"])
