# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the 10x Visium catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026; research: ../RESEARCH.md):
- kits follow the barcode list a run's reads can show, which follows the slide, not the assay:
  visium_v1slide (Visium v1/v2 slides, list visium-v1/v2: fresh-frozen and FFPE without CytAssist),
  visium_cytassist (CytAssist slides, lists visium-v4/v5), visium_hd (two-part HD barcodes);
- a study qualifies when its protocol names Visium or CytAssist, names no other vendor's kit and
  no 10x single-cell kit (Chromium, GEM-X, Next GEM, Multiome, Flex: such studies file
  single-cell libraries beside the spatial ones), and does not cite both the fresh-frozen
  (CG000239) and an FFPE (CG000407, CG000495) user guide;
- read shape: R1 at least 28 (barcode 16 + UMI 12; HD's R1 is 43) and R2 at least 40;
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

QUERY = 'library_construction_protocol="*visium*" OR library_construction_protocol="*cytassist*"'
VISIUM = re.compile(r"visium|cytassist", re.I)
SINGLE_CELL = re.compile(r"chromium|gem-x|single[- ]cell (3|5|gene expression|immune)|next ?gem|"
                         r"multiome|\bflex\b", re.I)
FRESH_FROZEN_GUIDE = re.compile(r"CG000239")
FFPE_GUIDE = re.compile(r"CG000407|CG000495")
HD = re.compile(r"visium\s*hd|visium high[- ]definition", re.I)
CYTASSIST = re.compile(r"cytassist", re.I)


def classify(protocol):
    named = VISIUM.search(protocol)
    if not named or catalogue.names_other_vendor(protocol, "10x") or SINGLE_CELL.search(protocol):
        return None
    if FRESH_FROZEN_GUIDE.search(protocol) and FFPE_GUIDE.search(protocol):
        return None
    for kit, rule in (("visium_hd", HD), ("visium_cytassist", CYTASSIST)):
        if (match := rule.search(protocol)):
            return kit, snippet(protocol, match)
    return "visium_v1slide", snippet(protocol, named)


def shape(r1_length, r2_length):
    return r1_length >= 28 and r2_length >= 40


def main(out_path):
    examined = catalogue.examined_studies(EXAMINED)
    catalogue.write(catalogue.draw(QUERY, classify, shape, examined_studies=examined), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
