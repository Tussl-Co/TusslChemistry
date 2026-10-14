# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draw the frozen head-to-head selection from both labelled pools, stratified by (Arc's
read-derived whitelist x Arc's metadata label), excluding every sample of the first Arc benchmark
(we tuned on those, so this set is held out).

Usage: python select.py <pool.tsv> <pool_recent.tsv> <first_selection.tsv> <seed> <out.tsv>
"""

import csv
import random
import sys

# (arc_whitelist, tech_10x, how many); every available row when fewer exist.
STRATA = [
    ("3M-february-2018", "3_prime_gex", 180),
    ("737K-august-2016", "5_prime_gex", 110),
    ("737K-august-2016", "3_prime_gex", 80),
    ("3M-3pgex-may-2023", "3_prime_gex", 30),
    ("3M-5pgex-jan-2023", "5_prime_gex", 30),
    # Out of scope by Arc's label: no tool should confidently claim a 10x GEX layout unless the
    # library really is GEX (adjudicated from ENA titles).
    ("3M-february-2018", "feature_barcoding", 20),
    ("737K-august-2016", "feature_barcoding", 10),
    ("737K-august-2016", "vdj", 20),
    ("unassigned", "multiome", 15),
    ("3M-february-2018", "multiome", 5),
    ("737K-august-2016", "other", 10),
    ("3M-february-2018", "other", 10),
]

PLAN = """# Head-to-head plan, frozen before running (7 Oct 2026).
# Tools: tussl-chemistry (commit in the selection's git history) vs Sanger cellgeni
#   infer_10x_chemistry.py @ a8611621 (2026-06-26), --dry-run --allow-missing-whitelists, same files.
#   Neither tool holds the 10x v1, Multiome GEX or ATAC lists.
# Held out: no sample from bench/arc/results/selection.tsv.
# Truth axes, scored separately:
#   whitelist: Arc's STARsolo-selected list from h5ad barcodes (read-derived, independent of both tools)
#   end (3'/5'): adjudicated from evidence independent of both tools, decided by these rules:
#     T = ENA experiment_title/library_name/sample_title says 5' (5'|5 prime|5prime|5p|five) or 3' alone;
#     N = nuclear-gene strand: R2 reads matching spliced transcripts of ~40 highly expressed nuclear
#         genes, sense vs antisense; a call needs 3x and 0.2% of reads (sense = 3', antisense = 5').
#     truth = T and N when they agree, or the one present; conflict or neither = undetermined.
#     Arc's tech_10x (an LLM reading metadata) is reported against, not used as truth.
#   library type: ENA titles naming TCR|BCR|VDJ -> vdj; ADT|HTO|CITE|hashtag|antibody|feature|CRISPR|
#     sgRNA|guide -> feature; else gex. Out-of-scope by Arc but gex by title counts as in scope.
# Metrics: whitelist accuracy (in-scope GEX); 3'/5' resolution rate and accuracy against truth;
#   no-call rate and cause; confident (HIGH) wrong calls; CPU seconds per sample, as each tool ships.
# Mapping: ours VERDICT_MAP (bench/arc/scripts/run_bench.py); Sanger 10x_3p_v2_or_5p_v1_v2 -> (737K, unresolved),
#   10x_3p_v3_family -> (3M-feb, 3'), 10x_3p_v4_gemx -> (3pgex, 3'), 10x_5p_v3_gemx -> (5pgex, 5').
"""


def main(pool, recent, first, seed, out_path):
    held = {l.split("\t", 1)[0] for l in open(first) if not l.startswith("#")}
    rows = {}
    for path in (pool, recent):
        for r in csv.DictReader(open(path), delimiter="\t"):
            if r["arc_whitelist"] != "error" and r["srx"] not in held:
                rows.setdefault(r["srx"], r)
    rows = list(rows.values())
    rng = random.Random(seed)
    chosen = []
    for wl, tech, n in STRATA:
        cands = sorted((r for r in rows if r["arc_whitelist"] == wl and r["tech_10x"] == tech), key=lambda r: r["srx"])
        rng.shuffle(cands)
        chosen += cands[:n]
        print(f"{wl:20} {tech:18} available {len(cands):5}  drawn {min(n, len(cands))}")
    with open(out_path, "w") as out:
        out.write(PLAN + f"# seed {seed}\n")
        out.write("srx\ttech_10x\tarc_whitelist\tobs_count\n")
        for r in chosen:
            out.write(f"{r['srx']}\t{r['tech_10x']}\t{r['arc_whitelist']}\t{r['obs_count']}\n")
    print(len(chosen), "selected")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), sys.argv[5])
