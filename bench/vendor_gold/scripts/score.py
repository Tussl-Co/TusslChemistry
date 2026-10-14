# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Scores the vendor gold set: each call against the chemistry the vendor states. The mapping from
the vendor's wording to the expected layout was fixed before any result was seen (f1d2ea1); its
wording parser was then corrected for three reading errors found when auditing the scored rows
(see the comments in `expected`). The truth (each vendor's own words) is unchanged. Kits the
tool does not support yet (3' v1) are scored only on not being called confidently as something
else. 10x Multiome GEX and Flex were in that group until they were added (8 Oct 2026); the same
wording now expects them. Every Flex dataset in the catalogue is Flex v1 (737K-fixed-rna-profiling,
an 8-base probe barcode); GEM-X Flex v2 has its own list.

Usage: python score.py <results.tsv>
"""

import collections
import csv
import re
import sys

UNSUPPORTED = "unsupported"


def expected(vendor, chemistry, dataset=""):
    c = chemistry.lower()
    v = vendor.lower()
    if "10x" in v:
        # Flex before Multiome: Flex guides mention "Feature Barcode", which the frozen version's
        # bare "arc" (for the ARC list) wrongly matched, as it did "GemCode barcode".
        if "flex" in c or "fixed rna" in c or "probe" in c:
            return "10x Flex v1"
        if "multiome" in c or re.search(r"\barc\b|arc-v1", c):
            return "10x Multiome GEX"
        if "gemcode" in c or "i7 gemcode barcode" in c:
            return UNSUPPORTED + " (3' v1)"
        # 10x sold 5' gene expression v1 and v1.1 as "Single Cell V(D)J Reagent Kits", and only
        # 5' kits had a v1.1. A 5 must be followed by a prime mark, not a space: guide numbers
        # such as "CG000185 RevB" end in 5.
        five = re.search(r"5['’ʹ]|5 ?prime|\b5p\b|five|v\(d\)j|v1\.1", c)
        gemx = "gem-x" in c or "gemx" in c
        if five:
            return "10x 5' v3 (GEM-X)" if gemx or re.search(r"\bv3\b", c) else "10x 5' v1/v2"
        if gemx or re.search(r"\bv4\b", c):
            return "10x 3' v4 (GEM-X)"
        if re.search(r"\bv3(\.1)?\b", c):
            return "10x 3' v3/v3.1"
        if re.search(r"\bv2\b", c):
            return "10x 3' v2"
        if re.search(r"\bv1\b", c):
            return UNSUPPORTED + " (3' v1)"
        return None
    if "bd" in v:
        # Added 8 Oct 2026 with the Enhanced V2/V3 layout: BD's newer beads, sold with WTA Next
        # and OMICS-One, and BD's own "Enhanced Beads V3" wording.
        if re.search(r"omics-one|wta next|tcr/bcr next|enhanced beads? v[23]", c):
            return "BD Rhapsody (Enhanced V2/V3 beads)"
        return "BD Rhapsody (Enhanced Beads)" if "enhanced" in c else "BD Rhapsody (original V1 beads)"
    if "parse" in v:
        return "Parse Evercode WT"
    if "scale" in v:
        # Added 8 Oct 2026 with the QuantumScale layout: ScaleBio names its QuantumScale test sets
        # in the dataset name ("quantumV1", "QS-") rather than the chemistry column.
        if re.search(r"quantum|\bqs-", c + " " + dataset.lower()):
            return "ScaleBio QuantumScale RNA"
        return "ScaleBio RNA"
    if "singleron" in v:
        for version in ("V1", "V2", "V3"):
            if re.search(rf"\b{version.lower()}\b", c):
                return f"Singleron GEXSCOPE {version}"
    return None


# The wording rules, pinned: each line is (vendor, chemistry, dataset, expected). Includes the three
# reading errors found in Phase A's audit and the mappings added since; a change to `expected`
# that breaks one stops the scorer before it scores anything.
WORDING_CASES = [
    ("10x Genomics", "Flex Gene Expression", "", "10x Flex v1"),
    ("10x Genomics", "Feature Barcode technology for Cell Surface Protein, Fixed RNA", "", "10x Flex v1"),
    ("10x Genomics", "Single Cell Multiome ATAC + Gene Expression v1", "", "10x Multiome GEX"),
    ("10x Genomics", "GemCode barcode, Single Cell 3' v1", "", UNSUPPORTED + " (3' v1)"),
    ("10x Genomics", "Chromium Single Cell V(D)J Reagent Kits", "", "10x 5' v1/v2"),
    ("10x Genomics", "Feature Barcoding (CG000185 RevB) Single Cell 3' v3", "", "10x 3' v3/v3.1"),
    ("10x Genomics", "GEM-X Single Cell 5' v3", "", "10x 5' v3 (GEM-X)"),
    ("BD Rhapsody", "BD OMICS-One WTA Next Assay", "", "BD Rhapsody (Enhanced V2/V3 beads)"),
    ("BD Rhapsody", "WTA Demo Datasets(Enhanced Beads)", "", "BD Rhapsody (Enhanced Beads)"),
    ("BD Rhapsody", "WTA Demo Datasets(V1 Beads)", "", "BD Rhapsody (original V1 beads)"),
    ("ScaleBio", "ScaleBio RNA", "ScaleBio_test_202506_quantumV1", "ScaleBio QuantumScale RNA"),
    ("ScaleBio", "ScaleBio RNA", "ScaleBio_test_202312_rnaV1.1", "ScaleBio RNA"),
    ("Singleron", "Singleron GEXSCOPE V3", "", "Singleron GEXSCOPE V3"),
]


def wording_sound():
    for vendor, chemistry, dataset, want in WORDING_CASES:
        got = expected(vendor, chemistry, dataset)
        if got != want:
            raise AssertionError(f"wording rule broken: {chemistry!r} -> {got!r}, want {want!r}")


def main(path):
    wording_sound()
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    tally = collections.Counter()
    for r in rows:
        want = expected(r["vendor"], r["chemistry"], r["dataset"])
        if r["status"] != "ok":
            outcome = "not run"
        elif want is None:
            outcome = "unscored (chemistry wording not mapped)"
        elif want.startswith(UNSUPPORTED):
            confident = r["verdict"] == "match" and r["confidence"] == "high"
            outcome = f"{want}: " + ("CONFIDENTLY CALLED " + r["layout"] if confident else "not confidently called (pass)")
        elif r["verdict"] == "match" and r["layout"] == want:
            outcome = f"right ({r['confidence']})"
        elif r["verdict"] == "match":
            outcome = f"WRONG ({r['confidence']}): called {r['layout']}"
        else:
            outcome = f"no call ({r['verdict']}: {r['candidates']})"
        key = outcome if not outcome.startswith(("unsupported", "no call")) else outcome.split(":")[0]
        tally[key] += 1
        print(f"{r['dataset'][:44]:44} {r['vendor'][:8]:8} {(want or r['chemistry'])[:26]:26} {outcome}")
    print("\n" + "\n".join(f"{n:4}  {k}" for k, n in tally.most_common()))


if __name__ == "__main__":
    main(sys.argv[1])
