# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the 10x single-cell ATAC catalogue (Multiome ATAC and scATAC) from ENA, by metadata
alone, before any run is read.

Rules, fixed before selection (8 Oct 2026; research: ../RESEARCH.md):
- ATAC-seq runs whose protocol names a 10x word (10x Genomics, Chromium, Cell Ranger, a bare
  10x that is not a buffer) or Multiome, and no other vendor's or method's kit (bare "BD" and
  "Bio-Rad" name sorters and reagents in genuine 10x protocols, so they do not refuse);
- multiome_atac when the protocol names Multiome or Cell Ranger ARC (or its kit numbers), else
  scatac when it names single-cell/nucleus ATAC; scATAC v1, v1.1 and v2 share one list and one
  read structure, so one kit;
- the cell barcode is the i5 index read, which ENA keeps only in 3- or 4-file runs: exactly one
  file of 16, 24 or 26 bases (the barcode, Multiome's with an 8-base spacer), exactly two of at
  least 30 (the genomic reads), any other an 8- or 10-base sample index; a 3-file run with a bare
  accession file (ENA's unpaired reads) has lost its barcode read and is left out;
- studies examined during research (results/examined.tsv) are never held out.

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw: per-study cap, design split, examined studies)

EXAMINED = os.path.join(HERE, "..", "results", "examined.tsv")
WORDS = ["*chromium*", "*10x*", "*multiome*", "*scATAC*", "*snATAC*", "*cellranger*",
         "*cell*ranger*", "*single*cell*ATAC*", "*single*nucle*"]
QUERY = ('library_strategy="ATAC-seq" AND ('
         + " OR ".join(f'library_construction_protocol="{w}"' for w in WORDS) + ")")
TENX = re.compile(r"10x genomics|chromium|cell ?ranger|\b10x\b(?! ?(buffer|td|pbs|tagment))", re.I)
OTHER = re.compile(r"rhapsody|evercode|parse bio|scale ?bio|surecell|ddseq|hydrop|sci-?atac|"
                   r"dscatac|dnbelab|\bmgi\b|pipseq|snare-seq|share-seq|paired-seq|s3-atac|txci|"
                   r"spear-atac|issaac", re.I)
MULTIOME = re.compile(r"multiome|multi-ome|atac\s*\+\s*(gene expression|gex)|cell ?ranger[- ]arc|"
                      r"1000283|1000285", re.I)
SCATAC = re.compile(r"single[- ](cell|nucle(i|us)) atac|\bs[cn]-?atac|chromium.{0,60}atac|"
                    r"cell ?ranger[- ]atac", re.I)
BARCODE_LENGTHS = {16, 24, 26}
INDEX_LENGTHS = {8, 10}
GENOMIC_MIN = 30


def snippet(protocol, match):
    start = max(0, match.start() - 40)
    return " ".join(protocol[start:match.end() + 40].split())


def classify(protocol):
    if OTHER.search(protocol) or not (TENX.search(protocol) or MULTIOME.search(protocol)):
        return None
    for kit, rule in (("multiome_atac", MULTIOME), ("scatac", SCATAC)):
        if (match := rule.search(protocol)):
            return kit, snippet(protocol, match)
    return None


def shape(names, lengths):
    """One barcode file, two genomic files, any rest sample index; no bare accession file."""
    bare = any(re.fullmatch(r"[SED]RR\d+\.fastq\.gz", n.rsplit("/", 1)[-1]) for n in names)
    barcodes = sum(n in BARCODE_LENGTHS for n in lengths)
    genomic = sum(n >= GENOMIC_MIN for n in lengths)
    rest_are_index = all(n in BARCODE_LENGTHS or n >= GENOMIC_MIN or n in INDEX_LENGTHS
                         for n in lengths)
    return not bare and barcodes == 1 and genomic == 2 and rest_are_index


def main(out_path):
    examined = catalogue.examined_studies(EXAMINED)
    rows = catalogue.draw(QUERY, classify, shape, files=(3, 4), examined_studies=examined)
    catalogue.write(rows, out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
