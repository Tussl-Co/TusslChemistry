# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the Takara ICELL8 3' DE catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026; research: ../RESEARCH.md):
- only ICELL8 3' DE carries its cell barcode in the reads (full-length ICELL8 and Shasta use index
  pairs, out of scope); a study qualifies when its protocol names ICELL8 and the 3' DE kit, and
  not ATAC, CUT&Tag, small RNA, TCR, 5', full-length wording or another vendor's kit;
- read shape: R1 21 to 30 bases (barcode 11 + UMI 10 or 14, then polyT) and R2 at least 50;
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


def snippet(protocol, match):
    start = max(0, match.start() - 40)
    return " ".join(protocol[start:match.end() + 40].split())

QUERY = 'library_construction_protocol="*icell8*"'
ICELL8 = re.compile(r"\bi-?cell ?8\b|wafergen", re.I)
REFUSE = re.compile(r"atac|cut ?& ?tag|small rna|pscsr|\btcr\b|5\s*['′’]", re.I)
FULL_LENGTH = re.compile(r"full[- ]length|smart-?seq (pro )?icell8|smart-?seq pro|shasta total", re.I)
THREE_PRIME = re.compile(r"3\s*['′’]?\s*-?\s*DE\b|\b30 DE\b|3['′’]\s*(end|differential)|"
                         r"differential expression|3['′’]? transcriptome", re.I)


def classify(protocol):
    if not ICELL8.search(protocol) or REFUSE.search(protocol) or FULL_LENGTH.search(protocol):
        return None
    if catalogue.names_other_vendor(protocol, "icell8"):
        return None
    match = THREE_PRIME.search(protocol)
    return ("icell8_3de", snippet(protocol, match)) if match else None


def shape(r1_length, r2_length):
    return 21 <= r1_length <= 30 and r2_length >= 50


def main(out_path):
    examined = catalogue.examined_studies(EXAMINED)
    catalogue.write(catalogue.draw(QUERY, classify, shape, examined_studies=examined), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
