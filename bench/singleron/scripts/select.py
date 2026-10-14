# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the Singleron GEXSCOPE catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (9 Oct 2026):
- a study qualifies when its library protocol names Singleron or GEXSCOPE and no other vendor's
  kit;
- a protocol naming the GEXSCOPE or Singleron MicroBead kit (not any "microbead": Miltenyi's
  cell-separation MicroBeads matched the first version of this rule, 9 Oct 2026) is
  singleron_microbead: GEXSCOPE MicroBead is structurally Drop-seq and its expected call is
  Drop-seq (src/layouts/singleron.rs);
- otherwise a version explicitly attached to a Singleron or GEXSCOPE kit ("GEXSCOPE V2",
  "GEXSCOPE Single Cell RNA Library Kit V3", "SCOPE-chip ... V1") names singleron_v1, _v2 or
  _v3. A version of anything else does not count: CeleScope or other software, an instrument
  (Singleron Matrix v1.0.1), a reagent (Separation Solution V2), another Singleron product
  (AccuraCode), a catalogue number;
- otherwise the study is the family, singleron (any of the three layouts);
- evidence is a 120-character window of the protocol around the matched words, whitespace
  collapsed after slicing;
- a read shape: R1 at least 70 bases (V1's barcode structure needs 69, then polyT starts) and
  R2 at least 50;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study. No
  run is design: the layouts were designed on runs recorded in no catalogue. Every study is drawn
  alike, including the three that today's validation labelled from their reads alone.

Drawn first in scratch on 9 Oct 2026 (the draft also wrote refused.tsv and versions.tsv); the
MicroBead rule was narrowed after that draw and before any run was read, which relabelled 5
rows from 3 studies to the family and drew nothing new.

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw: paired runs, per-study cap, splits)

QUERY = ('library_construction_protocol="*Singleron*"'
         ' OR library_construction_protocol="*GEXSCOPE*"')
WINDOW = 120
MICROBEAD = re.compile(r"(?:gexscope|singleron)[\s\u00ae\u2122]{0,3}micro ?bead", re.I)
FAMILY = re.compile(r"singleron|gexscope|scope-?chip", re.I)
# Up to a clause end; a full stop inside a version number ("scopeV3.0.1") does not end it.
GAP = r"(?:[^.;]|\.(?=\d))"


def token(group):
    return rf"\bv\s?(?P<{group}>[1-3])\b(?!\.\d)"


VERSION = re.compile(
    # the kit word first: "GEXSCOPE ... Library Kit V2", "Singleron ... (kit V2)", "SCOPE-chip V1"
    rf"(?:gexscope|singleron|scope-?chip){GAP}{{0,70}}?(?:\bkit\b|chip)[^.;()]{{0,25}}?{token('a')}"
    rf"|(?:gexscope|scope-?chip)[\s®™-]{{0,3}}{token('b')}"
    # the version first: "Library Kit Cell/Tissue for Matrix V2 protocol GEXSCOPE"
    rf"|\bkit\b[^.;()]{{0,30}}?{token('c')}{GAP}{{0,40}}?(?:gexscope|singleron)",
    re.I)
GROUPS = ("a", "b", "c")
# Singleron products other than GEXSCOPE: a version beside them is theirs.
NOT_GEXSCOPE = re.compile(r"accuracode", re.I)


def window(protocol, start, end, keep=None):
    """WINDOW characters around protocol[start:end], whitespace collapsed after slicing; when
    the span is longer than WINDOW, the window ends at `keep` (the version token)."""
    if end - start > WINDOW and keep is not None:
        start = max(0, keep - WINDOW)
    else:
        start = max(0, (start + end) // 2 - WINDOW // 2)
    start = min(start, max(0, len(protocol) - WINDOW))
    return " ".join(protocol[start:start + WINDOW].split())


def versions(protocol):
    """Counted version matches: (version, its position, match), leaving out those beside another product."""
    found = []
    for m in VERSION.finditer(protocol):
        if NOT_GEXSCOPE.search(protocol[max(0, m.start() - 40):m.end()]):
            continue
        group = next(g for g in GROUPS if m.group(g))
        found.append((m.group(group), m.start(group), m))
    return found


def classify(protocol):
    if catalogue.names_other_vendor(protocol, "singleron"):
        return None
    bead = MICROBEAD.search(protocol)
    if bead:
        return "singleron_microbead", window(protocol, bead.start(), bead.end())
    found = versions(protocol)
    named = {v for v, _, _ in found}
    if len(named) == 1:
        v, at, m = found[0]
        return f"singleron_v{v}", window(protocol, m.start(), m.end(), at + 1)
    family = FAMILY.search(protocol)
    if family is None:  # the query matched, so this is unreachable
        return "singleron", window(protocol, 0, 0)
    return "singleron", window(protocol, family.start(), family.end())


def shape(r1_length, r2_length):
    return r1_length >= 70 and r2_length >= 50


def main(out_path):
    catalogue.write(catalogue.draw(QUERY, classify, shape, design=False), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
