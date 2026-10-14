# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Recounts every number in bench/second_opinion/README.md from the result files.

Usage: python summarize.py <results dir>
"""

import collections
import csv
import sys

SANGER = {"10x_3p_v2_or_5p_v1_v2": {"10x 3' v2", "10x 5' v1/v2"}, "10x_3p_v3_family": {"10x 3' v3/v3.1"},
          "10x_3p_v4_gemx": {"10x 3' v4 (GEM-X)"}, "10x_5p_v3_gemx": {"10x 5' v3 (GEM-X)"},
          "10x_multiome_arc_gex": {"10x Multiome GEX"}}
REFEREE_SPECIES = ("Homo sapiens", "Mus musculus")


def rows(p):
    return list(csv.DictReader(open(p), delimiter="\t"))


def labels(r):
    return {x.strip() for x in r["label"].split("|")} if r["label"] else set()


def end(layout):
    return "5'" if "5'" in layout else ("3'" if ("3'" in layout or "Multiome" in layout) else "")


def cellranger_layouts(name):
    """Cell Ranger's chemistry name -> the tchem layouts it is consistent with."""
    n = name.lower()
    if "arc" in n or "multiome" in n:
        return {"10x Multiome GEX"}
    if "3'" in name or "3p" in n:
        for v, lay in (("v4", "10x 3' v4 (GEM-X)"), ("v3", "10x 3' v3/v3.1"), ("v2", "10x 3' v2")):
            if v in n:
                return {lay}
    if "5'" in name or "5p" in n:
        return {"10x 5' v3 (GEM-X)"} if "v3" in n else {"10x 5' v1/v2"}
    return set()


def side(call, label, consistent):
    if not consistent:
        return "no call"
    c, l = call in consistent, bool(label & consistent)
    return "both" if c and l else "tchem" if c else "label" if l else "neither"


def main(d):
    rs = rows(f"{d}/second_opinions.tsv")
    print(f"disputed HIGH calls: {len(rs)}; reads on disk: {sum(r['heads'] == 'True' for r in rs)}")
    tenx = [r for r in rs if r["call"].startswith("10x") and labels(r) and all(l.startswith("10x") for l in labels(r))]
    print(f"\nSanger, 10x label and 10x call ({len(tenx)}):",
          dict(collections.Counter(side(r["call"], labels(r), SANGER.get(r["sanger_id"], set())) for r in tenx)))
    nc = [r for r in tenx if not SANGER.get(r["sanger_id"])]
    why = collections.Counter("variable read lengths" if "variable read lengths" in r["sanger_error"]
                              else "Flex (no Flex list)" if r["call"] == "10x Flex v1"
                              else r["sanger_error"][:60] for r in nc)
    print("  no call, why:", dict(why))
    other = [r for r in rs if r not in tenx]
    print("\nSanger, other disputes:", dict(collections.Counter(
        ("call 10x" if r["call"].startswith("10x") else "call other vendor") + " / " +
        ("label 10x" if labels(r) and all(l.startswith("10x") for l in labels(r)) else
         "label other vendor" if labels(r) else "no label claimed") + " -> " +
        sanger_finding(r) for r in other)))
    both = [r for r in tenx if r["sanger_id"] == "10x_3p_v2_or_5p_v1_v2" and side(r["call"], labels(r), SANGER[r["sanger_id"]]) == "both"]
    print(f"\n  of the {len(both)} it cannot split (one list), strand referee:",
          dict(collections.Counter(strand_side(r) for r in both)))
    ends = [r for r in rs if len({end(l) for l in labels(r)}) == 1 and end(r["call"])
            and "" not in {end(l) for l in labels(r)} and end(r["call"]) not in {end(l) for l in labels(r)}]
    print(f"\nStrand referee, 3' vs 5' disputes ({len(ends)}):", dict(collections.Counter(strand_side(r) for r in ends)))
    print("  species:", dict(collections.Counter(r["species"] or "not stated" for r in ends)))
    cr = rows(f"{d}/cellranger.tsv")
    print(f"\nCell Ranger, sampled human 10x disputes ({len(cr)}):", dict(collections.Counter(
        side(r["call"], labels(r), cellranger_layouts(r["cellranger"])) if r["status"] == "ok" else "failed: " + r["message"][:60]
        for r in cr)))
    print("  chemistries:", dict(collections.Counter(r["cellranger"] for r in cr)))
    secs = sorted(int(r["seconds"]) for r in cr if r["status"] == "ok")
    if secs:
        print(f"  wall seconds per run: median {secs[len(secs) // 2]}, range {secs[0]}-{secs[-1]}")


def sanger_finding(r):
    """What Sanger's checker said about a run's files: a 10x list, no 10x list, or a refusal."""
    if r["sanger_id"]:
        return "finds a 10x list"
    if "No file matched" in r["sanger_error"] or "No GEX-compatible" in r["sanger_error"]:
        return "finds no 10x list"
    return "refuses the files (variable read lengths)" if "variable read lengths" in r["sanger_error"] \
        else "refuses the files (" + r["sanger_error"][:40] + ")"


def strand_side(r):
    if r["species"] not in REFEREE_SPECIES:
        return "no referee for the species"
    if not r["strand_end"]:
        return "undetermined"
    return "tchem" if r["strand_end"] == end(r["call"]) else "label"


if __name__ == "__main__":
    main(sys.argv[1])
