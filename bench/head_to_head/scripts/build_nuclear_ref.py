# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Builds the nuclear-gene strand reference used to adjudicate 3' vs 5', independently of both tools
(neither uses nuclear genes): every annotated transcript of ~40 highly and ubiquitously expressed
human genes, spliced from GENCODE exons on GRCh38 and written in the gene's sense orientation.

Usage: python build_nuclear_ref.py <gencode.gtf.gz> <GRCh38.fa (with .fai)> <out.fa> [GENE,GENE,...]
  The optional list replaces the default genes (e.g. TCR/BCR constant regions for library typing).
"""

import collections
import gzip
import re
import sys

GENES = {
    "MALAT1", "NEAT1", "ACTB", "GAPDH", "B2M", "TMSB4X", "TMSB10", "FTL", "FTH1", "EEF1A1", "EEF2",
    "TPT1", "PTMA", "MYL6", "CFL1", "PFN1", "OAZ1", "UBA52", "NACA", "HNRNPA1", "VIM", "LGALS1",
    "S100A6", "SRGN", "TXNIP", "CD74", "HLA-A", "HLA-B", "RPLP0", "RPLP1", "RPL3", "RPL10", "RPL13",
    "RPL41", "RPS2", "RPS12", "RPS18", "RPS27A", "H3-3B", "ACTG1",
}
COMP = str.maketrans("ACGTN", "TGCAN")


def fetch(fa, fai, chrom, start, end):
    """1-based inclusive [start, end] from an indexed FASTA."""
    _, length, offset, line_bases, line_width = fai[chrom]
    s = start - 1
    seek = offset + (s // line_bases) * line_width + s % line_bases
    fa.seek(seek)
    need = end - start + 1
    raw = fa.read(need + need // line_bases + 2).replace(b"\n", b"")
    return raw[:need].decode().upper()


def main(gtf, ref, out, genes=None):
    wanted = set(genes.split(",")) if genes else GENES
    fai = {l.split("\t")[0]: list(map(int, l.split("\t")[1:5])) and [l.split("\t")[0]] + list(map(int, l.split("\t")[1:5]))
           for l in open(ref + ".fai")}
    exons = collections.defaultdict(list)
    strand_of, gene_of = {}, {}
    with gzip.open(gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if c[2] != "exon":
                continue
            gene = re.search(r'gene_name "([^"]+)"', c[8]).group(1)
            if gene not in wanted:
                continue
            tx = re.search(r'transcript_id "([^"]+)"', c[8]).group(1)
            exons[tx].append((c[0], int(c[3]), int(c[4])))
            strand_of[tx], gene_of[tx] = c[6], gene
    with open(ref, "rb") as fa, open(out, "w") as o:
        for tx, ex in sorted(exons.items()):
            ex.sort(key=lambda e: e[1])
            seq = "".join(fetch(fa, fai, ch, s, e) for ch, s, e in ex)
            if strand_of[tx] == "-":
                seq = seq.translate(COMP)[::-1]
            o.write(f">{tx} {gene_of[tx]}\n{seq}\n")
    found = {gene_of[t] for t in exons}
    print(f"{len(exons)} transcripts from {len(found)} genes; missing: {sorted(wanted - found)}")


if __name__ == "__main__":
    main(*sys.argv[1:5])
