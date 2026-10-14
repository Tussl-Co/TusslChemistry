# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the inDrops (inDrop) catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (9 Oct 2026):
- the query: "indrop" anywhere in the library protocol (ENA wildcards hold no spaces; the search
  ignores case);
- a study qualifies when its library protocol names no vendor's kit in catalogue.VENDORS
  (inDrops is not one of them, so every vendor counts) and names inDrop or inDrops as the method
  (`\\bin-?drops?(?:tm|\\d+)?(?![a-z])`, or "1CellBio"; revision 1), and names no other
  single-cell method as one used (revisions 2, 2b): Drop-seq (not "Drop-seq apparatus"),
  CEL-Seq2 (not "modified CEL-Seq2"), a bare "10x" (not "10x with", a trituration count, nor
  "10X [word] Buffer", a concentration) refuse as "other method named"; Smart-seq refuses as "other
  method named (unclear)", since it may describe a modified inDrop library prep. "CEL-Seq"
  alone is inDrop chemistry (its linear amplification) and does not refuse;
- a protocol stating both v2 and v3 for inDrop refuses as "mixed versions" (revision 3);
- the kit is indrops_v3 when a v3 token ("v3", "V 3", "version 3", 'version "v3') follows an
  inDrop or 1CellBio mention within 60 characters of one clause ("inDrops v3", "inDrop Design
  V3", "inDrop platform with v3 chemistry", "inDrops V3 hydrogel beads"), or "version 3" or
  "v3 beads/hydrogel" lies within 200 characters of an inDrop mention ("Hydrogel beads with
  version 3 oligonucleotide design"); otherwise indrops
  (v1/v2, or no version stated). tchem's inDrops layout (src/layouts/other.rs) is v1/v2: in v3
  the barcodes are split across R1 and an index read, so tchem does not claim v3. Every version
  token seen is recorded in versions.tsv;
- evidence is a verbatim 120-character window of the protocol around the inDrop mention (and the
  version token, when one named the kit);
- a read shape: one read at least 45 bases and both at least 25. The barcode read needs about 47
  (BC1 8-11, the 22-base W1 linker, BC2 8, UMI 6) before polyT; either file may carry it, since
  inDrops runs are often archived swapped;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study. No
  run is design: the layout was designed on runs recorded in no catalogue. Studies with any run
  put through tchem in the wild sweep (examined_wild.tsv) are marked examined, never held out.

Drawn first in scratch on 9 Oct 2026, with refused.tsv (every study the draw leaves out, and
every run the read shape refused) and examined.tsv beside the catalogue.

Revised the same day, after that first draw (kept as catalogue_v1.tsv, refused_v1.tsv) and before
any read was run:
1. the inDrop word rule became `\\bin-?drops?(?:tm|\\d+)?(?![a-z])`: `\\b` failed on "inDropTM
   barcoded beads" (PRJNA996596) and "an inDrops21 setup" (PRJNA679887). It still refuses
   "Findrop" and "indroplets". The version rule's mention uses the same pattern;
2. protocols naming another method used beside inDrop are refused: "followed inDrop, Drop-seq
   and 10X standart protocol respectively" (PRJNA438523), "For CEL-Seq2 sample, cells were
   sorted" (PRJNA473915), "10x cells were prepared following the manufacturer's protocol"
   (PRJNA493298); Smart-seq2 (PRJNA698490, PRJNA860187) is refused as unclear, so it stays
   visible;
3. a study stating both inDrops v2 and v3 is refused as mixed: "Most samples were run using
   'inDrops v2' ... one experiment used 'inDrops v3' beads" (PRJNA471420).
2b. (9 Oct 2026, after the second draw, kept as catalogue_v2.tsv and refused_v2.tsv, and before
   any read was run) revision 2 fired where no other method was used: "10X Danieau Buffer"
   (PRJNA445487), "InDrop library protocol, modified CEL-Seq2" (PRJNA497610), "inDrop(TM)
   barcoded beads in a home-made microfluidic device connected to the Drop-seq apparatus"
   (PRJNA996596). The pattern became `drop-?seq(?! apparatus)|(?<!modified )cel-?seq2|
   \\b10x\\b(?!\\s+with)(?!\\s+(?:\\w+\\s+)?buffer)`.
TruDrop (Southard-Smith et al., BMC Genomics 2020, https://pmc.ncbi.nlm.nih.gov/articles/PMC7331155/) keeps inDrop V2's barcode read,
"[cell barcode 1] GAGTGATTGCTTGTGACGCCTT [cell barcode 2] [UMI] TTTT...", moved to the P5 (TruSeq
Read 1) side; TruDrop studies stay indrops.

Usage: python -I select.py <out dir> <examined_wild.tsv>
"""

import csv
import os
import re
import sys
import time

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw: paired runs, per-study cap, splits)

QUERY = 'library_construction_protocol="*indrop*"'
WINDOW = 120
MENTION = r"\bin-?drops?(?:tm|\d+)?(?![a-z])|1cellbio"
INDROP = re.compile(MENTION, re.I)
OTHER_METHOD = re.compile(r"drop-?seq(?! apparatus)|(?<!modified )cel-?seq2"
                          r"|\b10x\b(?!\s+with)(?!\s+(?:\w+\s+)?buffer)", re.I)
UNCLEAR_METHOD = re.compile(r"smart-?seq", re.I)
# A version token after an inDrop mention, within one clause.
VERSION = re.compile(rf"(?:{MENTION})[^.;]{{0,60}}?"
                     r'\b(?:version\s*"?v?\s?|v\.?\s?)(?P<v>[1-3])\b', re.I)
# A bead or design version within 200 characters of an inDrop mention, across sentences:
# "inDrop protocol on a custom system [...]. Hydrogel beads with version 3 oligonucleotide design".
NEAR = 200
BEADS = re.compile(r'\bversion\s*"?v?\s?(?P<v>[1-3])\b|\bv(?P<w>[1-3])\s(?:beads|hydrogel)', re.I)
SHAPE_REFUSED = "shape"
REFUSED_FIELDS = ["study", "run", "reason", "lengths", "protocol_window"]


def window(protocol, start, end):
    """A verbatim WINDOW-character window centred on protocol[start:end]."""
    middle = (start + end) // 2
    begin = max(0, min(middle - WINDOW // 2, len(protocol) - WINDOW))
    return protocol[begin:begin + WINDOW]


def versions(protocol):
    """Version tokens tied to inDrop: (version, match), in protocol order."""
    found = [(m.group("v"), m) for m in VERSION.finditer(protocol)]
    mentions = [m.start() for m in INDROP.finditer(protocol)]
    for m in BEADS.finditer(protocol):
        if any(m.end() > at - NEAR and m.start() < at + NEAR for at in mentions) \
                and not any(f.start() <= m.start() < f.end() for _, f in found):
            found.append((m.group("v") or m.group("w"), m))
    return sorted(found, key=lambda vm: vm[1].start())


def refusal(protocol):
    """Why `protocol` does not name inDrops, with the span to quote; (None, span) when it does."""
    for vendor in catalogue.VENDORS.values():
        hit = vendor.search(protocol)
        if hit:
            return "other vendor", (hit.start(), hit.end())
    hit = INDROP.search(protocol)
    if hit is None:
        raw = re.search("indrop", protocol, re.I)
        return "no inDrop statement", (raw.start(), raw.end()) if raw else (0, 0)
    other = OTHER_METHOD.search(protocol)
    if other:
        return "other method named", (other.start(), other.end())
    unclear = UNCLEAR_METHOD.search(protocol)
    if unclear:
        return "other method named (unclear)", (unclear.start(), unclear.end())
    found = versions(protocol)
    if {"2", "3"} <= {v for v, _ in found}:
        return "mixed versions", (found[0][1].start(), found[-1][1].end())
    return None, (hit.start(), hit.end())


def classify(protocol):
    reason, span = refusal(protocol)
    if reason is not None:
        return None
    v3 = [m for v, m in versions(protocol) if v == "3"]
    if v3:
        return "indrops_v3", window(protocol, v3[0].start(), v3[0].end())
    return "indrops", window(protocol, *span)


def shape(r1_length, r2_length):
    return max(r1_length, r2_length) >= 45 and min(r1_length, r2_length) >= 25


def ena_rows_once_retried(query):
    """ENA's rows, tried a second time once on failure; a second failure stops the draw."""
    try:
        return catalogue.ena_rows(query)
    except OSError as e:
        print(f"ENA failed ({e}); retrying once", file=sys.stderr)
        time.sleep(30)
        return catalogue.ena_rows(query)


def read_examined(path):
    with open(path) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def refused_rows(ena, drawn, lengths):
    """Every study of the query that gives no catalogue row, and every run the shape refused."""
    studies = {}
    for row in ena:
        studies.setdefault(row["study_accession"], []).append(row)
    drawn_runs = {r["run"] for r in drawn}
    drawn_studies = {r["study"] for r in drawn}
    out = []
    for study, runs in studies.items():
        paired = [r for r in runs if len([f for f in r["fastq_ftp"].split(";") if f]) == 2]
        # as catalogue.draw does: the first paired run's protocol (in ENA's order) names the kit
        first = paired[0] if paired else runs[0]
        paired.sort(key=lambda r: r["run_accession"])
        protocol = first["library_construction_protocol"]
        reason, span = refusal(protocol)
        quote = " ".join(window(protocol, *span).split())
        if reason is None and not paired:
            reason = "no fastq" if all(not r["fastq_ftp"] for r in runs) else "not paired"
        if reason is not None:
            out.append({"study": study, "run": "", "reason": reason, "lengths": "",
                        "protocol_window": quote})
            continue
        for r in paired[:catalogue.SHAPE_CHECKS_PER_STUDY]:
            urls = ["https://" + f for f in sorted(r["fastq_ftp"].split(";"))]
            if r["run_accession"] in drawn_runs or not all(u in lengths for u in urls):
                continue
            out.append({"study": study, "run": r["run_accession"], "reason": SHAPE_REFUSED,
                        "lengths": ",".join(str(lengths[u]) for u in urls),
                        "protocol_window": quote})
        if study not in drawn_studies:
            out.append({"study": study, "run": "", "reason": SHAPE_REFUSED,
                        "lengths": "", "protocol_window": quote})
    return out


def write_tsv(path, fields, rows, header=""):
    with open(path, "w", newline="") as f:
        f.write(header)
        writer = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main(out_dir, examined_path):
    ena = ena_rows_once_retried(QUERY)
    catalogue.ena_rows = lambda query: ena  # one ENA read serves the draw and the refusals
    lengths = {}
    measure = catalogue.read_length

    def recorded(url):
        lengths[url] = measure(url)
        return lengths[url]

    catalogue.read_length = recorded
    wild = read_examined(examined_path)
    examined = frozenset(r["study"] for r in wild)
    rows = catalogue.draw(QUERY, classify, shape, design=False, examined_studies=examined)
    catalogue.write(rows, os.path.join(out_dir, "catalogue.tsv"))

    drawn_studies = {r["study"] for r in rows}
    write_tsv(os.path.join(out_dir, "examined.tsv"), ["run", "study", "what"],
              [r for r in wild if r["study"] in drawn_studies],
              "# Runs put through tchem in the wild sweep before this catalogue was drawn (9 Oct\n"
              "# 2026); their studies are never held out (bench/heldout/scripts/catalogue.py).\n")
    write_tsv(os.path.join(out_dir, "refused.tsv"), REFUSED_FIELDS,
              refused_rows(ena, rows, lengths))
    protocols = {}
    for r in ena:
        protocols.setdefault(r["study_accession"], r["library_construction_protocol"])
    found = []
    for study in sorted(drawn_studies):
        kit = next(r["kit"] for r in rows if r["study"] == study)
        for v, m in versions(protocols[study]):
            found.append({"study": study, "kit": kit, "version": v, "token": m.group(0)[-40:]})
    write_tsv(os.path.join(out_dir, "versions.tsv"), ["study", "kit", "version", "token"], found)


if __name__ == "__main__":
    main(*sys.argv[1:])
