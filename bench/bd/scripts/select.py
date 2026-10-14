# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the BD Rhapsody catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026):
- runs first public from 2024 (the Enhanced V2/V3 bead era) whose library protocol names
  Rhapsody and no other vendor's kit;
- archive metadata almost never states the bead version, so every row is one kit, bd_rhapsody,
  scored as the BD family; bead versions were settled on BD's own demo runs, so no archive run
  is a design run;
- a read shape: R1 at least 46 bases (an Enhanced bead's cell label and UMI) and R2 at least 40;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study.

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw)

QUERY = 'library_construction_protocol="*Rhapsody*" AND first_public>=2024-01-01'
RHAPSODY = re.compile(r"rhapsody", re.I)


def classify(protocol):
    if catalogue.names_other_vendor(protocol, "bd"):
        return None
    match = RHAPSODY.search(protocol)
    start = max(0, match.start() - 40)
    return "bd_rhapsody", " ".join(protocol[start:match.end() + 60].split())


def shape(r1_length, r2_length):
    return r1_length >= 46 and r2_length >= 40


def main(out_path):
    catalogue.write(catalogue.draw(QUERY, classify, shape, design=False), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
