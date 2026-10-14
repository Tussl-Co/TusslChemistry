# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the PIPseq (Fluent, now Illumina Single Cell 3' RNA Prep) catalogue from ENA, by
metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026):
- a study qualifies when its library protocol names PIPseq or Fluent, exactly one chemistry
  version (v3, v4 or v5), and no other vendor's kit;
- a read shape: R1 at least 45 bases (four barcode parts, linkers and UMI) and R2 at least 50;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study,
  the earliest study's first run per version for design.

Revised 9 Oct 2026, before any of the affected studies was run through the tool: the version
token also reads "version" as well as "v", any decimals ("v4.05"), an optional "plus" ("V4plus",
"v4.0 PLUS"), and Fluent's kit catalogue number FBS-SCR-T<size>-<version>-V<rev> wherever it
appears (it names the kit by itself). Reason: an audit of ENA runs whose protocol mentions
PIPseq found three single-version v4 studies the old token missed: PRJNA1228638 ("Pipseq T2
V4plus"), PRJNA1278736 ("PIPseq 3' T20 version 4.0") and PRJNA1217330 ("FBS-SCR-T20-4-V4.05",
past a period that ended the window). On every study the query returns, no other
classification changes.

Revised again 9 Oct 2026 (second draw), before any of the added studies was run through the tool:
a study whose protocol states no version token, but names PIPseq V by name ("PIPseq V",
"PIP-Seq V", "PIPseq-V", "PIPseqV5", Fluent's "Pipseek V5") or Illumina's renamed kit ("Single
Cell 3' RNA Prep"), is v5; the query also asks for "Single Cell 3*RNA Prep". Reason: the PIPseq V
research (bench/pipseq/research_v5) found single-kit V studies the token cannot see, among them
PRJNA1284893 ("PIP-Seq V T20 kit"), PRJNA1477052 ("PIPseq-V T10 kit") and PRJNA1481985 ("Illumina
Single Cell 3' RNA Prep, T10 kit"). A study that states one version keeps it, so no study already
drawn changes kit: PRJNA1472464 ("V4 T20 PIP-seq kit (formerly Fluent Bio, now the Illumina Single
Cell 3'RNA Prep Kit)") stays v4. Rows from this draw are appended for new studies only
(SELECTION.md).

Usage: python select.py <out.tsv>
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw)

QUERY = ('library_construction_protocol="*PIPseq*" OR library_construction_protocol="*PIP-seq*"'
         ' OR library_construction_protocol="*Fluent*Bio*"'
         ' OR library_construction_protocol="*Single Cell 3*RNA Prep*"')
TOKEN = r"(?:v|version)\s?([3-5])(?:\.\d+)?(?:\s?plus)?\b"
VERSION = re.compile(r"(?:pip-?seq|fluent)[^.;]{0,60}?\b" + TOKEN +
                     r"|\b" + TOKEN + r"[^.;]{0,40}?pip-?seq"
                     r"|\bFBS-SCR-T\d+-([3-5])-V\d", re.I)
# PIPseq V by name, read only when no version token is stated (the second 9 Oct revision).
V_NAMED = re.compile(r"pip-?seq(?:\u2122|tm)?[ -]?V(?:\s?5(?:\.\d+)?)?(?![\w.])"
                     r"|pip-?seq\W{0,2}V5|pipseek\s?V5|single[- ]cell 3\W{0,3}\s?RNA Prep", re.I)


def classify(protocol):
    if catalogue.names_other_vendor(protocol, "pipseq"):
        return None
    matches = list(VERSION.finditer(protocol))
    found = {m.group(1) or m.group(2) or m.group(3) for m in matches}
    if not found:
        named = V_NAMED.search(protocol)
        if named is None:
            return None
        start = max(0, named.start() - 40)
        return "pipseq_v5", " ".join(protocol[start:named.end() + 40].split())
    if len(found) != 1:
        return None
    match = matches[0]
    start = max(0, match.start() - 40)
    return f"pipseq_v{found.pop()}", " ".join(protocol[start:match.end() + 40].split())


def shape(r1_length, r2_length):
    return r1_length >= 45 and r2_length >= 50


def main(out_path):
    catalogue.write(catalogue.draw(QUERY, classify, shape), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
