# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draw the frozen Arc benchmark selection from the labelled pool, stratified by
(Arc's read-derived whitelist x Arc's metadata-derived 3'/5' label), oversampling rare kits.

Usage: python select.py <pool.tsv> <seed> <out selection.tsv>
"""

import csv
import random
import sys

# (arc_whitelist, tech_10x) -> how many to draw; every available row when fewer exist.
STRATA = [
    ("3M-february-2018", "3_prime_gex", 20),
    ("737K-august-2016", "3_prime_gex", 15),
    ("737K-august-2016", "5_prime_gex", 15),
    ("3M-5pgex-jan-2023", "5_prime_gex", 12),
    ("3M-3pgex-may-2023", "3_prime_gex", 12),
    # Metadata disagrees with the reads' whitelist family: candidates for "Arc's label is wrong".
    ("3M-february-2018", "5_prime_gex", 6),
    ("3M-5pgex-jan-2023", "3_prime_gex", 4),
    # Out of scope: score only "did not confidently claim a wrong layout".
    ("737K-arc-v1", "multiome", 5),
    ("unassigned", None, 5),
    # Real 10x look-alikes that are not gene expression: receptor and feature-barcode libraries.
    ("737K-august-2016", "vdj", 4),
    ("3M-february-2018", "feature_barcoding", 4),
]

PLAN = """# Frozen before running (7 Oct 2026). Truth axes, scored separately:
#   whitelist: Arc's STARsolo-selected list, recovered from h5ad cell barcodes (read-derived)
#   end (3'/5'): Arc's tech_10x from SRAgent (an LLM reading SRA metadata; may itself be wrong)
# Verdict mapping: see VERDICT_MAP in run_bench.py; AMBIGUOUS 3' v2/5' v2 -> (737K, ambiguous).
# Outcome categories after inspecting every disagreement:
#   agree | tool wrong | Arc wrong | data unusable (archive damage) | out of scope
# Out-of-scope rows (multiome, unassigned) pass when no 10x GEX layout is claimed at HIGH.
"""


def main(pool_path, seed, out_path):
    rows = [r for r in csv.DictReader(open(pool_path), delimiter="\t") if r["arc_whitelist"] != "error"]
    rng = random.Random(seed)
    chosen = []
    for wl, tech, n in STRATA:
        cands = [r for r in rows if r["arc_whitelist"] == wl and (tech is None or r["tech_10x"] == tech)]
        rng.shuffle(cands)
        chosen += cands[:n]
        print(f"{wl:20} {str(tech):12} available {len(cands):4}  drawn {min(n, len(cands))}")
    with open(out_path, "w") as out:
        out.write(PLAN + f"# seed {seed}, pool {pool_path.rsplit('/', 1)[-1]}\n")
        out.write("srx\ttech_10x\tarc_whitelist\tobs_count\n")
        for r in chosen:
            out.write(f"{r['srx']}\t{r['tech_10x']}\t{r['arc_whitelist']}\t{r['obs_count']}\n")
    print(len(chosen), "selected")


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
