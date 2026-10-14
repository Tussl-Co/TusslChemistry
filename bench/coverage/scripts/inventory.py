# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""What each bench's catalogue holds, and what the coverage run takes from it: the first step of
one validation procedure for every kit (ROADMAP.md, "Trust"). It reads committed files and git
only, never the network, and states only what a file or commit records. Where no record exists it
says so; it never infers a status from a mention.

1. Catalogues: columns, `split`, runs and studies by split, `examined.tsv`, and how many coverage
   cases sources.py draws from each.
2. Examined records, both ways: every run in an examined.tsv is named in the bench's RESEARCH.md,
   and every catalogue row marked `examined` has its study in examined.tsv. Design rows of those
   catalogues, and whether RESEARCH.md names them.
3. Statuses held in code (sources.DESIGN, small_cases, NOT_CLAIMED), each cited commit checked.
4. Rule 2: the drawn cases set aside because their study gave a design run for a layout they
   test (sources.shaped_by_design), and the counts with and without them (tally.py's own tally).
5. sources.py's cases now against the last coverage run's calls.tsv, and the vendor gold set.

`split=design` in a drawn catalogue is set by the draw rule (bench/heldout/scripts/catalogue.py)
or by moving a run after a fix. This script does not check that those reads shaped the layout.

Usage: python inventory.py
"""

import collections
import glob
import os
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sources  # noqa: E402

BENCH = sources.BENCH
REPO = os.path.dirname(BENCH)
tally = sources._load("coverage_tally", os.path.join(BENCH, "coverage", "scripts", "tally.py"))


def rows(path):
    return list(sources.rows(path))

CALLS = os.path.join(BENCH, "coverage", "results", "calls.tsv")
ADJUDICATION = os.path.join(BENCH, "coverage", "results", "adjudication.tsv")
VENDOR_RESULTS = os.path.join(BENCH, "vendor_gold", "results", "results_quantumscale.tsv")
COMMIT = re.compile(r"(?<=[(, ])([0-9a-f]{7,40})(?=[),])")  # "(81da1b8)", "(2d0140e, validate.sh)"
NO_DATA = "/nonexistent"  # sources.cases() only joins paths under it; nothing is read from it
# The sources coverage reads, by bench name: filled from sources.cases() in main, so a catalogue
# coverage starts reading is reported without editing this file.
SOURCE_OF = {}
LONG_DESIGN = os.path.join(BENCH, "long_read", "results", "catalogue.tsv")
LONG_HELDOUT = os.path.join(BENCH, "long_read", "results", "catalogue_heldout.tsv")


def columns(path):
    with open(path, newline="") as f:
        for line in f:
            if not line.startswith("#"):
                return line.rstrip("\n").split("\t")
    return []


def study_of(row):
    """The study accession, as sources.long_read_cases cuts it ("PRJNA1 / GSE2: title" -> PRJNA1)."""
    return re.split(r"[ :/]", row.get("study", ""), maxsplit=1)[0]


def named(accession, text):
    return re.search(rf"(?<![A-Za-z0-9]){re.escape(accession)}(?![0-9])", text) is not None


def commit_exists(ref):
    return subprocess.run(["git", "-C", REPO, "cat-file", "-e", f"{ref}^{{commit}}"],
                          capture_output=True).returncode == 0


def long_layout(kit):
    return tuple(layout for pattern, layout in sources.LONG_KIT if re.match(pattern, kit))[:1]


def catalogues():
    """(bench, path, rows, split of a row) for every catalogue, the long-read pair included."""
    out = []
    for path in sorted(glob.glob(f"{BENCH}/*/results/catalogue.tsv")):
        bench = path.split(os.sep)[-3]
        if bench in ("vendor_gold", "long_read"):
            continue  # datasets read through results*.tsv; the long-read pair is added below
        out.append((bench, path, rows(path), lambda r: r.get("split") or None))
    out.append(("long_read (design)", LONG_DESIGN, rows(LONG_DESIGN), lambda r: "design"))
    out.append(("long_read", LONG_HELDOUT, rows(LONG_HELDOUT), lambda r: None))
    return out


# 1 and 2 ------------------------------------------------------------------------------------------

def print_catalogues(cats, cases_by_source):
    print("## 1. Catalogues\n")
    print("| Bench | split column | Runs (studies) by split | examined.tsv runs (studies) | "
          "Coverage reads it | Coverage cases |")
    print("|---|---|---|---|---|---|")
    for bench, path, data, split in cats:
        by = collections.defaultdict(lambda: [0, set()])
        for r in data:
            by[split(r) or "(no split)"][0] += 1
            by[split(r) or "(no split)"][1].add(study_of(r))
        text = "; ".join(f"{k} {n} ({len(s)})" for k, (n, s) in sorted(by.items()))
        examined_path = os.path.join(os.path.dirname(path), "examined.tsv")
        examined = rows(examined_path) if os.path.exists(examined_path) else None
        ex = "none" if examined is None else f"{len(examined)} ({len({e['study'] for e in examined})})"
        source = SOURCE_OF.get(bench)
        cases = len(cases_by_source.get(source, [])) if source else 0
        print(f"| {bench} | {'yes' if 'split' in columns(path) else 'no'} | {text} | {ex} | "
              f"{'yes' if source else 'no'} | {cases} |")
    print()


def print_examined(cats):
    print("## 2. Examined and design records\n")
    print("| Bench | examined.tsv runs named in RESEARCH.md | `examined` rows whose study is in "
          "examined.tsv | design rows named in RESEARCH.md | examined studies held out |")
    print("|---|---|---|---|---|")
    unknown = []
    for bench, path, data, split in cats:
        examined_path = os.path.join(os.path.dirname(path), "examined.tsv")
        if not os.path.exists(examined_path):
            if any(split(r) in ("held_out", None) for r in data) and bench in SOURCE_OF:
                unknown.append(bench)
            continue
        examined = rows(examined_path)
        research_path = os.path.join(BENCH, bench, "RESEARCH.md")
        research = open(research_path).read() if os.path.exists(research_path) else ""
        in_research = sum(named(e["run"], research) for e in examined)
        studies = {e["study"] for e in examined}
        marked = [r for r in data if split(r) == "examined"]
        backed = sum(study_of(r) in studies for r in marked)
        design = [r for r in data if split(r) == "design"]
        design_named = sum(named(r["run"], research) for r in design)
        leaked = sorted(studies & {study_of(r) for r in data if split(r) == "held_out"})
        print(f"| {bench} | {in_research} / {len(examined)} | {backed} / {len(marked)} | "
              f"{design_named} / {len(design)} | {len(leaked)} |")
    print()
    print("Benches whose held-out rows feed coverage but which have no record of what was read "
          f"while their layouts were designed (no examined.tsv): {', '.join(unknown)}. "
          "Whether other runs of theirs were read is unknown.\n")


# 3 ------------------------------------------------------------------------------------------------

def print_code_statuses(cats):
    where = {}
    for bench, _, data, split in cats:
        for r in data:
            where.setdefault(r.get("run", ""), f"{bench} ({split(r) or 'no split'})")
    print("## 3. Statuses held in code, not in a catalogue (`sources.py`)\n")
    print("| Where | Run | Says | Commit cited | Commit exists | In a catalogue |")
    print("|---|---|---|---|---|---|")
    entries = [("DESIGN", run, why) for run, why in sources.DESIGN.items()]
    entries += [(f"small_cases ({c.source})", c.run, f"held out: {c.referee}")
                for c in sources.small_cases(NO_DATA)]
    for kind, run, says in entries:
        refs = COMMIT.findall(says)
        exists = ", ".join("yes" if commit_exists(c) else "NO" for c in refs) or "-"
        print(f"| {kind} | {run} | {says} | {', '.join(refs) or 'none'} | {exists} | "
              f"{where.get(run, 'no')} |")
    print(f"\nNOT_CLAIMED kits (no layout is right for them): {', '.join(sorted(sources.NOT_CLAIMED))}\n")


# 4 ------------------------------------------------------------------------------------------------

def counts(calls, drop):
    kept = [c for c in calls if c["run"] not in drop]
    adjudications = {r["run"]: r for r in rows(ADJUDICATION) if r["run"] not in drop}
    return {k: (e["exact"], len(e["studies"]))
            for k, e in tally.tally(kept, adjudications).items()}


def print_design_overlap(drawn):
    """Rule 2's effect: the drawn cases set aside because their study gave a design run for a layout
    they test (sources.shaped_by_design), and the counts with and without them."""
    calls = rows(CALLS)
    by_run = {c["run"]: c for c in calls}
    hits = sources.shaped_by_design(drawn, sources.design_studies())
    print("## 4. Rule 2: drawn cases from a study that gave a design run for the same layout\n")
    print("| Source | Runs | Study | Design run (record) | Call at the last coverage run |")
    print("|---|---|---|---|---|")
    groups = collections.defaultdict(list)
    for case, what in hits:
        groups[(case.source, case.study, what)].append(case.run)
    for (source, study, what), runs in sorted(groups.items()):
        said = collections.Counter(
            f"{by_run[r]['confidence'] or by_run[r]['verdict']} {by_run[r]['layout']}".strip()
            for r in runs if r in by_run)
        call_text = "; ".join(f"{n}x {s}" for s, n in sorted(said.items()))
        print(f"| {source} | {len(runs)} | {study} | {what} | {call_text} |")
    print()
    drop = {case.run for case, _ in hits}
    before, after = counts(calls, set()), counts(calls, drop)
    print("Exact-kit counts with and without them (tally.py on the last coverage run's calls.tsv, "
          "nothing re-run):\n")
    print("| Layout | Reads | Exact runs now | Studies now | Runs without | Studies without |")
    print("|---|---|---|---|---|---|")
    for key in sorted(before):
        a = after.get(key, (0, 0))
        if before[key] != a:
            print(f"| {key[0]} | {key[1]} | {before[key][0]} | {before[key][1]} | {a[0]} | {a[1]} |")
    print()


# 5 ------------------------------------------------------------------------------------------------

def print_cross_check(cases_by_source):
    called = collections.Counter(c["source"] for c in rows(CALLS))
    with open(CALLS) as f:
        header = f.readline().strip().lstrip("# ")
    print("## 5. Cases now against the last coverage run\n")
    print(f"Last run: {header}\n")
    print("| Source | Cases sources.py yields now | Rows in calls.tsv |")
    print("|---|---|---|")
    for source in sorted(set(cases_by_source) | set(called)):
        now, then = len(cases_by_source.get(source, [])), called.get(source, 0)
        print(f"| {source} | {now} | {then}{'' if now == then else ' (differs)'} |")
    read = rows(VENDOR_RESULTS)
    names = {r["dataset"] for r in read}
    others = sorted(p for p in glob.glob(f"{BENCH}/vendor_gold/results/results*.tsv")
                    if p != VENDOR_RESULTS)
    missing = {os.path.basename(p): len({r["dataset"] for r in rows(p)} - names) for p in others}
    why = collections.Counter()
    for r in read:
        if r["dataset"] in sources.DESIGN:
            why["a design run (sources.DESIGN)"] += 1
        elif r["status"] != "ok":
            why[f"status {r['status']}"] += 1
        elif sources.vendor_expected(r["vendor"], r["chemistry"], r["dataset"]) is None:
            why["no expected layout from the vendor's wording"] += 1
    print(f"\nVendor gold: sources.py reads {os.path.basename(VENDOR_RESULTS)} ({len(names)} "
          "datasets). Datasets in the other results files missing from it: "
          + ", ".join(f"{k} {n}" for k, n in sorted(missing.items())) + ". Left out of coverage: "
          + ", ".join(f"{n} {k}" for k, n in sorted(why.items())) + ".\n")


def main():
    drawn = sources.drawn(NO_DATA)
    cases_by_source = collections.defaultdict(list)
    for case in sources.cases(NO_DATA):
        cases_by_source[case.source].append(case)
    SOURCE_OF.update({name: name for name in cases_by_source})
    invariant(SOURCE_OF, "coverage reads some source")
    cats = catalogues()
    invariant(len(cats) >= 12, f"only {len(cats)} catalogues found")
    print_catalogues(cats, cases_by_source)
    print_examined(cats)
    print_code_statuses(cats)
    print_design_overlap(drawn)
    print_cross_check(cases_by_source)


if __name__ == "__main__":
    if len(sys.argv) != 1:
        sys.exit(__doc__)
    main()
