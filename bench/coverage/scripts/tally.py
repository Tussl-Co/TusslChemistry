# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Turns a coverage run (calls.tsv) into each layout's validation: of the held-out runs the binary
called that layout at HIGH, how many were right, at what resolution, from how many studies, how
many are disputed and how many wrong. Writes data/validation.tsv, which the binary's table must
match (a test there checks it); --check compares instead of writing.

What is counted is the precision of a HIGH call, the question a user asks of a verdict ("it says
L at HIGH: how often is that wrong?"). Only HIGH calls count; MEDIUM and LOW are not vouched for.

  exact_runs    HIGH calls of the layout on runs whose truth names exactly that layout
  exact_studies the studies those runs come from (a run without a known study is its own)
  family_runs   HIGH calls on runs whose truth admits the layout among others (kit family only)
  read_label_runs HIGH calls right only against a label derived from the reads themselves (Arc's
                barcode-list matching, linker positions): the same kind of signal tchem uses, so
                not an independent check, and not counted in exact_runs or exact_studies (A4,
                adopted 9 Oct 2026; bench/VALIDATION.md, section 4)
  disputed      calls whose label disagrees and only the reads' own evidence supports (adjudicated)
  wrong         calls whose truth does not admit the layout, unless adjudicated otherwise

Levels and bounds are the binary's alone (src/layouts/validation.rs): this script counts, and the
README's coverage table is rendered from data/validation.tsv by a test there.

Usage: python tally.py <calls.tsv> <adjudication.tsv> <validation.tsv> [--check]
"""

import collections
import csv
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

SEPARATOR = " | "
COLUMNS = ["layout", "read_kind", "exact_runs", "exact_studies", "family_runs", "read_label_runs",
           "disputed", "wrong"]
# Referees whose label comes from the reads' own evidence (sources.py, adjudication.tsv): a call
# right against one of these shows two read matchers agree, not that the lab's record agrees.
READ_DERIVED = ("Arc's", "linker positions in its reads")
# Kit families whose labels rarely name the version or bead type, though the reads separate the
# layouts. The family record asks: of the HIGH calls naming any layout of a family, how many runs
# are that family? Same names as src/layouts/validation.rs FAMILIES (tests/tidy.rs checks).
FAMILIES = {
    "BD Rhapsody": ("BD Rhapsody (original V1 beads)", "BD Rhapsody (Enhanced Beads)",
                    "BD Rhapsody (Enhanced V2/V3 beads)"),
    "Parse Evercode WT": ("Parse Evercode WT", "Parse Evercode WT v3"),
    "SeekGene SeekOne DD": ("SeekGene SeekOne DD 3'", "SeekGene SeekOne DD 5'"),
    "Singleron GEXSCOPE": ("Singleron GEXSCOPE V1", "Singleron GEXSCOPE V2",
                           "Singleron GEXSCOPE V3"),
}
FAMILY_OF = {layout: family for family, layouts in FAMILIES.items() for layout in layouts}
def rows(path):
    with open(path, newline="") as f:
        lines = [line for line in f if not line.startswith("#")]
    return list(csv.DictReader(lines, delimiter="\t"))


def comment(path):
    with open(path) as f:
        return f.readline().strip()


def independent(referee):
    """Whether a referee's label is independent of the reads (not READ_DERIVED)."""
    invariant(referee, "every label names its referee")
    return not referee.startswith(READ_DERIVED)


def outcome(call, adjudications):
    """right/disputed/wrong for one HIGH call, the resolution it was right at, and whether the
    label that decided it is independent of the reads."""
    decided = adjudications.get(call["run"])
    if decided is not None:
        truth = tuple(t for t in decided["truth"].split(SEPARATOR) if t)
        invariant(decided["outcome"] != "right" or call["layout"] in truth, f"{call['run']}: adjudicated right, but the truth given does not admit {call['layout']}")
        resolution = "exact" if len(truth) == 1 else "family"
        return decided["outcome"], resolution, independent(decided["referee"])
    truth = tuple(t for t in call["truth"].split(SEPARATOR) if t)
    result = "right" if call["layout"] in truth else "wrong"
    return result, call["resolution"], independent(call["referee"])


def tally(calls, adjudications):
    counts = collections.defaultdict(lambda: {"exact": 0, "studies": set(), "family": 0,
                                              "read_label": 0, "disputed": 0, "wrong": []})
    used = set()
    for call in calls:
        if call["status"] != "ok" or call["confidence"] != "high":
            continue
        invariant(call["verdict"] == "match" and call["layout"], f"{call['run']}: HIGH without a layout")
        result, resolution, from_record = outcome(call, adjudications)
        if call["run"] in adjudications:
            used.add(call["run"])
        entry = counts[(call["layout"], call["read_kind"])]
        if result == "right" and not from_record:
            entry["read_label"] += 1
        elif result == "right" and resolution == "exact":
            entry["exact"] += 1
            entry["studies"].add(call["study"] or f"run:{call['run']}")
        elif result == "right":
            entry["family"] += 1
        elif result == "disputed":
            entry["disputed"] += 1
        else:
            invariant(result == "wrong", f"{call['run']}: unknown outcome {result}")
            entry["wrong"].append(call)
    stale = set(adjudications) - used
    invariant(not stale, f"adjudicated runs no longer called HIGH: {sorted(stale)}")
    return counts


def family_outcome(call, adjudications):
    """right/disputed/wrong at the family level for a HIGH call naming a layout of a family, and
    whether its label is independent of the reads. Right when the label admits any layout of the
    same family: a disputed exact call (a BD label naming other beads) is right here."""
    family = FAMILY_OF[call["layout"]]
    decided = adjudications.get(call["run"])
    truth = tuple(t for t in call["truth"].split(SEPARATOR) if t)
    referee = call["referee"]
    if decided is not None and decided["outcome"] == "right":
        truth = tuple(t for t in decided["truth"].split(SEPARATOR) if t)
        referee = decided["referee"]
    in_family = any(FAMILY_OF.get(t) == family for t in truth)
    if in_family:
        return "right", independent(referee)
    if decided is not None and decided["outcome"] == "disputed":
        return "disputed", independent(decided["referee"])
    return "wrong", independent(referee)


def tally_families(calls, adjudications):
    """The family record: one row per family and read kind with any HIGH call naming it."""
    counts = collections.defaultdict(lambda: {"exact": 0, "studies": set(), "family": 0,
                                              "read_label": 0, "disputed": 0, "wrong": []})
    for call in calls:
        if call["status"] != "ok" or call["confidence"] != "high" or call["layout"] not in FAMILY_OF:
            continue
        result, from_record = family_outcome(call, adjudications)
        entry = counts[(FAMILY_OF[call["layout"]], call["read_kind"])]
        if result == "right" and not from_record:
            entry["read_label"] += 1
        elif result == "right":
            entry["exact"] += 1
            entry["studies"].add(call["study"] or f"run:{call['run']}")
        elif result == "disputed":
            entry["disputed"] += 1
        else:
            invariant(result == "wrong", f"{call['run']}: unknown family outcome {result}")
            entry["wrong"].append(call)
    invariant(all(f in FAMILIES for f, _ in counts), "every family row names a family")
    return counts


def family_path(out_path):
    """data/validation_family.tsv beside data/validation.tsv."""
    head, tail = os.path.split(out_path)
    invariant(tail.endswith(".tsv"), f"{out_path}: a .tsv")
    return os.path.join(head, tail[:-len(".tsv")] + "_family.tsv")


def render(counts, source, what="each layout"):
    out = [f"# Validation of {what}'s HIGH calls on held-out runs (bench/coverage).",
           f"# From: {source}",
           "\t".join(COLUMNS)]
    for (layout, kind), e in sorted(counts.items()):
        invariant(e["exact"] >= len(e["studies"]), f"{layout}: more studies than runs")
        out.append("\t".join(map(str, [layout, kind, e["exact"], len(e["studies"]), e["family"],
                                       e["read_label"], e["disputed"], len(e["wrong"])])))
    return "\n".join(out) + "\n"


def report(counts):
    print(f"{'layout':44} {'kind':5} {'exact':>5} {'stud':>4} {'fam':>4} {'read':>4} {'disp':>4} "
          f"{'wrong':>5}")
    for (layout, kind), e in sorted(counts.items()):
        studies, wrong = len(e["studies"]), len(e["wrong"])
        print(f"{layout:44} {kind:5} {e['exact']:5} {studies:4} {e['family']:4} {e['read_label']:4} "
              f"{e['disputed']:4} {wrong:5}")
    for (layout, _), e in sorted(counts.items()):
        for call in e["wrong"]:
            print(f"WRONG  {layout} on {call['source']} {call['run']}: truth {call['truth']!r}")


def main(calls_path, adjudication_path, out_path, check=False):
    calls = rows(calls_path)
    adjudications = {r["run"]: r for r in rows(adjudication_path)}
    invariant(all(r["outcome"] in ("right", "disputed", "wrong") for r in adjudications.values()),
              "every adjudication is right, disputed or wrong")
    counts = tally(calls, adjudications)
    report(counts)
    families = tally_families(calls, adjudications)
    print("kit families (a HIGH call naming any layout of the family; independent labels only):")
    report(families)
    source = comment(calls_path).lstrip("# ")
    outputs = [(out_path, render(counts, source)), (family_path(out_path), render(families, source, "each kit family"))]
    if check:
        for path, text in outputs:
            with open(path) as f:
                if f.read() != text:
                    sys.exit(f"{path} is stale: rerun tally.py without --check")
        print(f"{out_path} and its family table match (the README tables are checked by cargo test)")
        return
    for path, text in outputs:
        with open(path, "w") as f:
            f.write(text)
    print("now: TCHEM_BLESS=1 cargo test the_readme_states  (renders the README's table)")


if __name__ == "__main__":
    if len(sys.argv) not in (4, 5) or sys.argv[4:] not in ([], ["--check"]):
        sys.exit(__doc__)
    main(*sys.argv[1:4], check=sys.argv[4:] == ["--check"])
