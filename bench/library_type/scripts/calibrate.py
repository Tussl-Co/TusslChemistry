# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Library-type signals on every head-to-head sample's cDNA read, to set thresholds from data:
  vdj:     reads with >= 2 sampled 25-mers from TCR/BCR constant-region transcripts (either strand)
  guide:   reads carrying a 12-mer of the sgRNA scaffold (either strand)
  feature: reads carrying a 10x feature-barcode capture sequence or the TotalSeq-A polyA-after-barcode shape

  hk:      reads with >= 2 sampled 25-mers from ~40 housekeeping transcripts (either strand)
  r2top:   share of reads whose first 15 bases are among the 200 commonest (tag panels are tiny)

Usage: python calibrate.py <h2h results.tsv> <adjudication.tsv> <vdj_constant.fa> <nuclear_genes.fa> <out.tsv>
"""
import csv, gzip, os, sys
from concurrent.futures import ProcessPoolExecutor

K, STEP, READS = 25, 5, 50_000
COMP = str.maketrans("ACGTN", "TGCAN")
SCAFFOLD = "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
# The optimised (F+E) scaffold used by many CROP-seq/Perturb-seq vectors.
SCAFFOLD_FE = "GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
CAPTURE = ["GCTTTAAGGCCGGTCCTAGCAA", "GCTCACCTATTAGCGGCTAAGG"]
rc = lambda s: s.translate(COMP)[::-1]

def kmers_of(fasta):
    seqs, cur = [], []
    for line in open(fasta):
        if line.startswith(">"):
            if cur: seqs.append("".join(cur)); cur = []
        else: cur.append(line.strip())
    if cur: seqs.append("".join(cur))
    s = {x[i:i + K] for x in seqs for i in range(len(x) - K + 1)}
    return s | {rc(k) for k in s}

def init(fasta, hk_fasta):
    global VDJ, SCAF, CAP, HK
    VDJ = kmers_of(fasta)
    HK = kmers_of(hk_fasta)
    SCAF = {x[i:i + 12] for x in (SCAFFOLD, SCAFFOLD_FE) for i in range(len(x) - 11)}
    SCAF |= {rc(k) for k in SCAF}
    CAP = CAPTURE + [rc(c) for c in CAPTURE]

def one(row):
    path = os.path.expandvars(row["cdna"])
    import collections
    n = vdj = guide = feat = hk = 0
    starts = collections.Counter()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 != 1: continue
            s = line.strip(); n += 1
            vdj += sum(s[j:j + K] in VDJ for j in range(0, len(s) - K + 1, STEP)) >= 2
            guide += any(s[j:j + 12] in SCAF for j in range(0, len(s) - 11, 3))
            feat += any(c in s for c in CAP) or s[15:25] == "AAAAAAAAAA"
            hk += sum(s[j:j + K] in HK for j in range(0, len(s) - K + 1, STEP)) >= 2
            starts[s[:15]] += 1
            if n >= READS: break
    n = max(n, 1)
    top = sum(c for _, c in starts.most_common(200))
    return (row["srx"], f"{100 * vdj / n:.2f}", f"{100 * guide / n:.2f}", f"{100 * feat / n:.2f}",
            f"{100 * hk / n:.2f}", f"{100 * top / n:.2f}")

def main(results, adjudication, fasta, hk_fasta, out):
    adj = {r["srx"]: r for r in csv.DictReader(open(adjudication), delimiter="\t")}
    rows = [r for r in csv.DictReader(open(results), delimiter="\t") if r["condition"] == "ok" and r.get("cdna")]
    with ProcessPoolExecutor(8, initializer=init, initargs=(fasta, hk_fasta)) as ex, open(out, "w") as o:
        o.write("srx\ttech_10x\tlib_type_title\ttitle\tvdj_pct\tguide_pct\tfeature_pct\thk_pct\tr2top_pct\n")
        for srx, v, g, f, h, t in ex.map(one, rows, chunksize=4):
            a = adj.get(srx, {})
            o.write(f"{srx}\t{a.get('tech_10x', '')}\t{a.get('lib_type', '')}\t{a.get('title', '')[:90]}\t{v}\t{g}\t{f}\t{h}\t{t}\n"); o.flush()

if __name__ == "__main__":
    main(*sys.argv[1:6])
