# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the 10x Flex (Fixed RNA Profiling) catalogue from ENA, by metadata alone, before any run
is read.

Rules, fixed before selection (9 Oct 2026; revised the same day, before any read was run,
after the first draw showed that ENA's "*Fixed*RNA*" matched almost nothing and that a bare
"Flex" admitted QuantStudio, Qubit and E8 Flex medium in 10x 3' protocols):
- the query: "Flex" anywhere, or "Fixed" with "10x" or "Chromium";
- a study qualifies when its library protocol names no other vendor's kit, names no Visium,
  CytAssist or spatial capture (Flex is a single-cell assay; Visium probe kits must not count),
  has a 10x context (10x, Chromium, Fixed RNA Profiling, Single Cell Gene Expression Flex or
  GEM-X Flex) and a Flex statement within one sentence: "Fixed RNA Profiling", "GEM-X Flex",
  "Flex" followed within 30 characters by an assay word (gene expression, assay, kit, probe,
  sample prep, single cell, chemistry, v1, v2), or "10x"/"Chromium" followed within 40
  characters by "Flex". A "Flex" within 25 characters after an instrument name (QuantStudio,
  Qubit, StepOne, LightCycler, ViiA, CFX), or followed within 15 by "medium"/"media", does not
  count;
- the kit is 10x_gemx_flex_v2 when the protocol says GEM-X Flex, Flex v2 or Flex Gene Expression
  v2, and 10x_flex_v1 otherwise;
- a read shape: both reads at least 28 bases and one at least 45. The barcode read (16-base
  barcode, 12-base UMI at its start) may be sequenced long, to 150 or more, and either file may
  carry it; the tool detects which. (Revised the same day: "R1 exactly 28" refused Flex runs
  with a long barcode read that the tool had called correctly.)
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study. No
  run is design: Flex's design data came from a 10x demo. No study was examined;
- merged: every row of the earlier hand-checked catalogue (GEO records, 28-base R1s seen then)
  whose study the draw does not hold, as held_out. The `origin` column says which (draw or
  old_catalogue). Given a catalogue already in this format, its old_catalogue rows are taken,
  so a redraw over its own output keeps them.

Usage: python select.py <out.tsv> [<old catalogue.tsv>]
"""

import csv
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw)

QUERY = ('(library_construction_protocol="*Flex*") OR (library_construction_protocol="*Fixed*"'
         ' AND (library_construction_protocol="*10x*"'
         ' OR library_construction_protocol="*Chromium*"))')
OLD = os.path.join(os.path.dirname(__file__), "..", "results", "catalogue.tsv")
OUT_FIELDS = catalogue.OUT_FIELDS + ["origin"]
SPATIAL = re.compile(r"visium|cytassist|spatial[- ]capture", re.I)
CONTEXT = re.compile(r"10x|chromium|fixed rna profiling|single cell gene expression flex"
                     r"|gem-x flex", re.I)
# A Flex statement, within one sentence ([^.!?] gaps). Each alternative is one named group.
FLEX = re.compile(r"(?P<frp>fixed rna profiling)"
                  r"|(?P<gemx>gem-?x flex)"
                  r"|(?P<word>\bflex\b)(?=[^.!?]{0,30}?(?:gene expression|assay|kit|probe"
                  r"|sample prep|single[- ]cell|chemistry|v1|v2))"
                  r"|(?:10x|chromium)[^.!?]{0,40}?(?P<after>\bflex\b)", re.I)
INSTRUMENT = re.compile(r"(?:quantstudio|qubit|stepone|lightcycler|viia|cfx)[^.!?]{0,25}$", re.I)
MEDIUM = re.compile(r"^.{0,15}?(?:medium|media)", re.I | re.S)
WORD = re.compile(r"\bflex\b", re.I)
V2 = re.compile(r"gem-x flex|\bflex v2\b|flex gene expression v2", re.I)
WINDOW = 120


def flex_statement(protocol):
    """The first match of FLEX whose "Flex" is not an instrument's or a medium's, or None."""
    for match in FLEX.finditer(protocol):
        if match.group("frp"):
            return match
        word = WORD.search(protocol, match.start(), match.end())
        if word and not INSTRUMENT.search(protocol[max(0, word.start() - 60):word.start()]) \
                and not MEDIUM.search(protocol[word.end():word.end() + 30]):
            return match
    return None


def refusal(protocol):
    """Why `protocol` does not name a Flex kit, or None when it does."""
    if catalogue.names_other_vendor(protocol, "10x"):
        return "other vendor"
    if SPATIAL.search(protocol):
        return "visium"
    if not CONTEXT.search(protocol):
        return "no 10x context"
    if flex_statement(protocol) is None:
        return "no flex statement"
    return None


def window(protocol, match):
    """A verbatim WINDOW-character window centred on `match`."""
    middle = (match.start() + match.end()) // 2
    start = max(0, min(middle - WINDOW // 2, len(protocol) - WINDOW))
    return protocol[start:start + WINDOW]


def classify(protocol):
    if refusal(protocol) is not None:
        return None
    kit = "10x_gemx_flex_v2" if V2.search(protocol) else "10x_flex_v1"
    return kit, window(protocol, flex_statement(protocol))


def shape(r1_length, r2_length):
    return min(r1_length, r2_length) >= 28 and max(r1_length, r2_length) >= 45


def read_old(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def old_rows(old, drawn_studies):
    """Rows of the earlier catalogue `old` whose study the draw does not hold, in OUT_FIELDS."""
    if old and "origin" in old[0]:  # already merged: keep what came from the old catalogue
        return [r for r in old if r["origin"] == "old_catalogue" and r["study"] not in drawn_studies]
    rows = []
    for r in old:
        if r["study"] in drawn_studies:
            continue
        urls = [u for u in (r["r1_url"], r["r2_url"]) if u]
        rows.append({"kit": r["kit"], "study": r["study"], "run": r["run"],
                     "r1_url": r["r1_url"], "r2_url": r["r2_url"], "split": "held_out",
                     "evidence": r["evidence"], "urls": ";".join(urls),
                     "origin": "old_catalogue"})
    return rows


def merge(drawn, old):
    """The draw, then the old catalogue's rows of studies the draw does not hold."""
    for row in drawn:
        row["origin"] = "draw"
    rows = drawn + old_rows(old, {r["study"] for r in drawn})
    runs = [r["run"] for r in rows]
    catalogue.invariant(len(runs) == len(set(runs)), "a run is listed once")
    catalogue.invariant(all(r["split"] == "held_out" for r in rows), "every row is held out")
    return rows


def write(rows, out_path):
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, OUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    for kit in sorted({r["kit"] for r in rows}):
        for origin in ("draw", "old_catalogue"):
            some = [r for r in rows if r["kit"] == kit and r["origin"] == origin]
            print(f"{kit} {origin}: {len(some)} runs from {len({r['study'] for r in some})} "
                  "studies")


def main(out_path, old_path=OLD):
    old = read_old(old_path)  # read before writing: out_path may be old_path
    write(merge(catalogue.draw(QUERY, classify, shape, design=False), old), out_path)


if __name__ == "__main__":
    main(*sys.argv[1:])
