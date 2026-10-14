# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the Bio-Rad ddSEQ / Illumina SureCell WTA 3' catalogue from ENA, by metadata alone.

Rules, fixed before selection (8 Oct 2026; research: ../RESEARCH.md):
- a study qualifies when its protocol names SureCell or ddSEQ, and not ATAC, Bio-Rad's newer
  ddSEQ Single-Cell 3' RNA-Seq kit (another structure) or another vendor's kit;
- read shape: R1 66 to 80 bases (phase 0-5, three 6-base blocks, two 15-base linkers, ACG, UMI 8,
  GAC, polyT) and R2 at least 50;
- research read the first run of every paired study on ENA, so every study is examined
  (results/examined.tsv): the catalogue holds a design run and no held-out run. ddSEQ needs
  unseen data from elsewhere before any validation is claimed.

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw: per-study cap, design split, examined studies)

EXAMINED = os.path.join(HERE, "..", "results", "examined.tsv")


def snippet(protocol, match):
    start = max(0, match.start() - 40)
    return " ".join(protocol[start:match.end() + 40].split())

QUERY = 'library_construction_protocol="*ddSEQ*" OR library_construction_protocol="*SureCell*"'
KIT = re.compile(r"SureCell|ddSEQ", re.I)
ATAC = re.compile(r"atac", re.I)
NEWER_KIT = re.compile(r"ddSEQ Single[- ]Cell 3.{0,2} RNA[- ]Seq|Omnition|deconvolution oligo", re.I)


def classify(protocol):
    if ATAC.search(protocol) or NEWER_KIT.search(protocol):
        return None
    if catalogue.names_other_vendor(protocol, "biorad"):
        return None
    match = KIT.search(protocol)
    return ("ddseq_wta3", snippet(protocol, match)) if match else None


def shape(r1_length, r2_length):
    return 66 <= r1_length <= 80 and r2_length >= 50


def main(out_path):
    examined = catalogue.examined_studies(EXAMINED)
    catalogue.write(catalogue.draw(QUERY, classify, shape, examined_studies=examined), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
