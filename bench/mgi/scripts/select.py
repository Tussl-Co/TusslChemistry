# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the MGI DNBelab C-series catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026):
- a study qualifies when its library protocol names DNBelab and either C4 (mgi_c4) or the
  high-throughput kit (mgi_c_ht: "HT" or "high-throughput" beside DNBelab or C-series);
  ATAC protocols and protocols naming both are left out;
- a gene-expression read shape: R1 26 to 40 bases (barcode read: 10 + 10 barcode, 10 UMI) and
  R2 at least 50 (cDNA). Added after the first draw, before any detection ran: its C4 design run
  had lost its barcode read on export (two 150-base cDNA reads) and its HT design run was the
  kit's oligo library (R1 20, R2 30), so metadata alone cannot draw this catalogue;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study,
  the earliest study's first run per kit for design.

Revision, 9 Oct 2026, before any run of the added studies was read or run through the tool. An
audit of ENA's 92 DNBelab studies found that the 8 Oct rules kept 32 and refused genuine C-series
gene-expression studies on wording alone. A protocol the 8 Oct rules classify keeps that label.
Only a protocol they refuse goes on to the clauses below:
- other vendor: a protocol that names another vendor's kit too (catalogue.names_other_vendor) is
  refused. The study may hold libraries of either kit. This applies to the new clauses only, so the
  8 Oct catalogue, including its C4 design study, stays as drawn;
- H, hyphens: U+2010 to U+2015 and U+2212 read as "-" before matching. One study wrote
  "High‑throughput" with U+2011, and the HT pattern missed it;
- B, C4 as the device: a protocol names both C4 and HT. If every C4 mention is the droplet device
  ("DNBelab C4 /DNBelabTaiM4 system", "C4 platform", "C4 chip"), the HT name is the reagent kit,
  and the study is mgi_c_ht. If C4 names a kit or prep set, the study stays refused;
- G, generic C-series: "DNBelab C Series" or bare "DNBelab C", with no C4 or HT token, is
  mgi_c_series. The protocol does not say which C-series RNA kit was used;
- T, the TaiM4 instrument: "C-TaiM4" / "C-TaiM 4", with no C4 or HT token, is mgi_c_series too;
- mgi_c_series shares the one layout ("MGI DNBelab C4 / C-series HT") and its design runs. It has
  no design run of its own: its rows are all held out. The 8 Oct design runs stay as they are;
- ATAC refusal is unchanged (clause A, "ATAC only when an MGI ATAC kit is named", was not adopted).

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
import catalogue  # noqa: E402  (shared ENA draw: paired runs, per-study cap, design split)
from invariant import invariant  # noqa: E402

QUERY = 'library_construction_protocol="*DNBelab*"'
C4 = re.compile(r"DNBelab[^.;]{0,30}?\bC4\b", re.I)
HT = re.compile(r"DNBelab[^.;]{0,40}?(?:\bHT\b|high[- ]throughput)|C[- ]?series[^.;]{0,20}?"
                r"(?:\bHT\b|high[- ]throughput)", re.I)
ATAC = re.compile(r"atac", re.I)
# 9 Oct 2026 clauses.
HYPHENS = re.compile("[‐-―−]")
C4_DEVICE = re.compile(r"C4\s*(?:/\s*(?:DNBelab\s*)?TaiM\s*4\s*)?(?:system|device|platform|chip|"
                       r"channel|instrument|microfluidic|scRNA (?:chip|slide)|droplet)", re.I)
GENERIC = re.compile(r"DNBelab(?:TM|™)?\s+C(?:[- ]?Series)?\b(?![- ]?(?:4|TaiM))", re.I)
TAIM = re.compile(r"C[- ]?TaiM\s?4?", re.I)
DESIGN_KITS = ("mgi_c4", "mgi_c_ht")  # kits with their own design run; mgi_c_series has none


def snippet(protocol, match):
    start = max(0, match.start() - 40)
    return " ".join(protocol[start:match.end() + 40].split())


def classify_8oct(protocol):
    """The 8 Oct rules, unchanged."""
    if ATAC.search(protocol):
        return None
    c4, ht = C4.search(protocol), HT.search(protocol)
    if c4 and not ht:
        return "mgi_c4", snippet(protocol, c4)
    if ht and not c4:
        return "mgi_c_ht", snippet(protocol, ht)
    return None


def c4_is_device(protocol):
    """Whether every C4 match names the droplet device, not a kit."""
    return all(C4_DEVICE.match(protocol, m.end() - 2) for m in C4.finditer(protocol))


def classify(protocol):
    named = classify_8oct(protocol)
    if named is not None:
        return named
    if catalogue.names_other_vendor(protocol, "mgi"):
        return None
    p = HYPHENS.sub("-", protocol)  # H
    named = classify_8oct(p)
    if named is not None:
        return named
    if ATAC.search(p):
        return None
    c4, ht = C4.search(p), HT.search(p)
    if c4 and ht:
        return ("mgi_c_ht", snippet(p, ht)) if c4_is_device(p) else None  # B
    if c4 or ht:
        return None  # unreachable: classify_8oct took these
    generic = GENERIC.search(p) or TAIM.search(p)  # G, T
    return ("mgi_c_series", snippet(p, generic)) if generic else None


def gene_expression_shape(r1_length, r2_length):
    return 26 <= r1_length <= 40 and r2_length >= 50


def split_without_c_series_design(rows):
    """The shared draw gives every kit a design run; mgi_c_series shares the layout's design runs,
    so its rows are all held out."""
    for row in rows:
        if row["kit"] not in DESIGN_KITS and row["split"] == "design":
            row["split"] = "held_out"
    for kit in {r["kit"] for r in rows}:
        designs = sum(r["kit"] == kit and r["split"] == "design" for r in rows)
        invariant(designs == (1 if kit in DESIGN_KITS else 0), f"{kit}: {designs} design runs")
    return rows


def draw():
    return split_without_c_series_design(catalogue.draw(QUERY, classify, gene_expression_shape))


def main(out_path):
    catalogue.write(draw(), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
