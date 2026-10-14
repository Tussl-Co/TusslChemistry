# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the CEL-seq2, SORT-seq and MARS-seq catalogue from ENA, by metadata alone, before any run
is read.

Rules, fixed before selection (8 Oct 2026; research: ../RESEARCH.md, which read lengths only):
- RNA-seq runs whose protocol names SORT-seq (sortseq), CEL-seq2 with its "2" (celseq2) or
  MARS-seq (marsseq), first match in that order; a bare CEL-seq is v1 and left out, as are
  scDam&T/DamID, VASA, HIV and reporter "sort-seq", Hi-C, methylome, ATAC and protocols naming
  another vendor's kit;
- read shape, in either file order (tchem detects swapped files): the shorter read 12 to 40 bases
  and the longer at least 35. Equal long pairs (76/76) are left out: lengths cannot show the
  barcode survived. (Research proposed 14 as SORT-seq's and MARS-seq's floor and 50 as MARS-seq's
  long read; one rule for all is kept, since the draw's shape does not know the kit.)

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

QUERY = ('(library_construction_protocol="*CEL*seq*" OR library_construction_protocol="*SORT*seq*"'
         ' OR library_construction_protocol="*MARS*seq*") AND library_strategy="RNA-Seq"')
REFUSE = re.compile(r"dam ?& ?t|damid|vasa|\bhiv\b|degenerate|sort_seq|reporter|hi-?c\b|pbat|"
                    r"methylom|atac", re.I)
KITS = [
    ("sortseq", re.compile(r"\bsort[- ]?seq\b|single cell discoveries|scdiscoveries", re.I)),
    ("celseq2", re.compile(r"\bcel[- ]?seq[- ]?2\b", re.I)),
    ("marsseq", re.compile(r"\bmars[- ]?seq", re.I)),
]


def classify(protocol):
    if REFUSE.search(protocol) or any(v.search(protocol) for v in catalogue.VENDORS.values()):
        return None
    for kit, rule in KITS:
        if (match := rule.search(protocol)):
            return kit, snippet(protocol, match)
    return None


def shape(r1_length, r2_length):
    short, long = sorted((r1_length, r2_length))
    return 12 <= short <= 40 and long >= 35


def main(out_path):
    examined = catalogue.examined_studies(EXAMINED)
    catalogue.write(catalogue.draw(QUERY, classify, shape, examined_studies=examined), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
