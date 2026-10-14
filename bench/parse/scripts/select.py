# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the Parse Evercode WT catalogue from ENA, by metadata alone, before any run is read.

Rules, fixed before selection (8 Oct 2026):
- a study qualifies when its library protocol names Parse/Evercode and exactly one WT version
  (v1, v2 or v3) next to "WT" or "whole transcriptome", and it has paired runs (two FASTQ
  files). A version beside "fixation" names the fixation kit, not the WT kit, and never counts;
- the shared rules in bench/heldout/scripts/catalogue.py: paired runs, at most 2 per study,
  the earliest study's first run per version for design.

Revision (9 Oct 2026), before any added study was run. It applies only to studies not already in
results/catalogue.tsv. That file is append-only: its rows, splits and design runs do not change.
The rows this revision adds are written by `main` as an appendix. Evidence: the audit of the 189
ENA studies whose protocol names Parse/Evercode. 111 of them were missing from the 8 Oct draw:
58 refused by classify (38 state no version, 19 state a version away from "WT", 1 states two),
28 not returned by QUERY, 17 with no paired run, and 8 that pass the old rules but were not drawn
(their ENA metadata does not show why).
1. QUERY becomes `*parse*` OR `*evercode*`. ENA does not return 28 audit studies for
   `*Parse*Bio*` or `*Evercode*`, although most say "Parse Biosciences" (e.g. PRJNA1094867,
   PRJNA781229, PRJNA820450). The broad query returns all 189. The 41 extra studies use "parse"
   in another sense, and classify refuses every one.
2. The version is no longer required. Held-out Parse runs are scored by family
   (bench/heldout/scripts/run.py, PARSE), so a stated version does not change how a run is scored.
   A protocol that names the WT kit with no single version gets the new kit `parse_wt`
   (e.g. PRJNA1006639 "Evercode Whole Transcriptome Mega kit"; PRJNA1176115 "Version 2 WT Kit").
   `parse_wt` has no design run, because both Parse layouts already have their design runs.
   With exactly one version, the kit is `parse_wt_v{n}` as before.
3. WT must be named. Any of these counts: Evercode WT / Whole Transcriptome, Evercode Mini/Mega,
   WT Mini/Mega/Penta/100k, "Parse WT", "WT kit", "whole transcript..." or "SPLiT-seq" within a
   clause of Parse/Evercode, or a WT catalogue number. A match that contains or spans "fixation"
   does not count. Protocols that name Parse with no product stay out ("Parse Biosciences
   protocol", "Parse Evercode kit", "Evercode V1": e.g. PRJNA1064612, PRJNA1240242, PRJNA1310582).
   Parse's TCR and BCR kits are also "Evercode", so these protocols do not prove WT.
4. Non-WT products are refused, even beside WT, because a study's runs may be either library:
   Evercode TCR/BCR, CRISPR Detect, Gene Select, guide libraries. PRJNA1455443 (TCR only),
   PRJNA1255131 (WT + TCR), PRJNA1427207 (CRISPR Detect) and PRJNA1244860 (CRISPR guide library)
   are refused. Four studies already in the catalogue drew CRISPR guide sublibraries
   (PRJNA1208204, PRJNA1304277, PRJNA1456597, PRJNA1463191).
5. `catalogue.names_other_vendor` now applies, because a study may hold libraries of either
   vendor: PRJNA1284247, PRJNA1216910, PRJNA1033995 (10x libraries beside Parse); PRJNA1458839
   (10x and HIVE); PRJNA1309117 (10x and ScaleBio ATAC). PRJNA975472 is already in the catalogue
   and drew 10x Visium runs. Three studies that only cite a 10x protocol (PRJNA1094867,
   PRJNA1169626, PRJNA993024) are refused too: the rule does not judge use.
6. Shape: R1 >= 25 and R2 >= 58. 58 bases is the structured R2 of the shortest Parse layout (v3).
   This drops runs whose barcode read is missing or too short (e.g. PRJNA1458839 SRR38306962,
   26/50, a 10x run).
7. PRJNA1390207 is refused by name: its first paired runs are bulk QIAseq UPXome RNA-seq
   (MiniSeq, "RNA-seq, Jurkat, CTLA4 RNP, Mock"), and its single-cell runs are Superb-seq,
   a modified Parse protocol.

Usage: python select.py <catalogue.tsv> <appended.tsv>
  Writes the rows of studies not in <catalogue.tsv>, all held out. Read lengths are fetched
  only for those studies.
"""

import csv
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "heldout", "scripts"))
import catalogue  # noqa: E402  (shared ENA draw: paired runs, per-study cap, design split)
from catalogue import invariant  # noqa: E402  (always-on assertion)

QUERY = 'library_construction_protocol="*parse*" OR library_construction_protocol="*evercode*"'
VERSION = re.compile(r"(?:\bwt\b|whole[- ]transcriptome)[^.;]{0,30}?\bv\s?([1-3])(?:\.\d)?\b"
                     r"|\bv([1-3])(?:\.\d)?\b[^.;]{0,20}?(?:\bwt\b|whole[- ]transcriptome)", re.I)
FIXATION = re.compile(r"fixation", re.I)
PARSE = re.compile(r"evercode|parse ?bio", re.I)
# Within one clause ("V3.0" does not end it), never across "fixation".
_GAP = r"(?:(?!fixation)(?:[^.;]|\.\d))"
_TM = r"\W{0,4}(?:tm\W{0,2})?"
WT_NAMED = re.compile(
    r"evercode" + _TM + r"(?:wt|whole[- ]transcriptome)(?![a-z])"
    r"|evercode" + _TM + r"(?:mini|mega)\b"
    r"|\bwt[ _-]?(?:mini|mega|penta|100k)\b"
    r"|parse(?: bio\w*)?\W{0,3}wt\b"
    r"|\bwt (?:library )?kits?\b"
    r"|(?:evercode|parse)" + _GAP + r"{0,60}?whole[- ]transcript\w*"
    r"|whole[- ]transcript\w*(?: \(wt\))?" + _GAP + r"{0,60}?(?:evercode|parse)"
    r"|(?:evercode|parse)" + _GAP + r"{0,60}?split[- ]?seq"
    r"|split[- ]?seq" + _GAP + r"{0,60}?(?:evercode|parse)"
    r"|\b(?:SB2001|(?:EC-?)?W0\d{4}|ECW0\d+|WT100\w*|WT200|UMWT\d+)\b", re.I)
OTHER_PRODUCT = re.compile(
    r"evercode" + _TM + r"(?:tcr|bcr)\b|\b(?:tcr|bcr)\b kit|crispr[- ]detect|gene select"
    r"|\bCRS\d{4}\b|crispr[^.;]{0,30}librar|guide rnas?[^.;]{0,40}librar|immune profiling", re.I)
REFUSED_STUDIES = {
    "PRJNA1390207": "first paired runs are bulk QIAseq UPXome RNA-seq; single-cell runs are "
                    "Superb-seq, a modified Parse protocol",
}
FIELDS = ["kit", "study", "run", "r1_url", "r2_url", "split", "evidence"]


def versions(protocol):
    return [m for m in VERSION.finditer(protocol) if not FIXATION.search(m.group(0))]


def version_of(protocol):
    found = {m.group(1) or m.group(2) for m in versions(protocol)}
    return found.pop() if len(found) == 1 else None


def window(protocol, match):
    start = max(0, match.start() - 40)
    return " ".join(protocol[start:match.end() + 40].split())


def wt_named(protocol):
    return [m for m in WT_NAMED.finditer(protocol) if not FIXATION.search(m.group(0))]


def classify(protocol):
    if not PARSE.search(protocol) or OTHER_PRODUCT.search(protocol):
        return None
    if catalogue.names_other_vendor(protocol, "parse"):
        return None
    named = wt_named(protocol)
    if not named:
        return None
    version = version_of(protocol)
    if version is None:
        return "parse_wt", window(protocol, named[0])
    return f"parse_wt_v{version}", window(protocol, versions(protocol)[0])


def shape(r1_length, r2_length):
    return r1_length >= 25 and r2_length >= 58


def appended(existing_studies):
    """catalogue.draw's rules, for the studies not in `existing_studies`: paired runs, the
    first study protocol classified, the shape on at most SHAPE_CHECKS_PER_STUDY runs, at most
    RUNS_PER_STUDY runs. Every row is held out: the existing catalogue keeps its design runs."""
    studies = {}
    for row in catalogue.ena_rows(QUERY):
        files = [f for f in row["fastq_ftp"].split(";") if f]
        study = row["study_accession"]
        if len(files) == 2 and study not in existing_studies and study not in REFUSED_STUDIES:
            studies.setdefault(study, []).append((row, files))
    rows = []
    for study, runs in studies.items():
        named = classify(runs[0][0]["library_construction_protocol"])
        if named is None:
            continue
        kit, evidence = named
        ordered = sorted(runs, key=lambda rf: rf[0]["run_accession"])
        for row, files in catalogue.shaped(ordered, shape, True):
            urls = ["https://" + f for f in sorted(files)]
            rows.append({"kit": kit, "study": study, "run": row["run_accession"],
                         "first_public": row["first_public"], "r1_url": urls[0],
                         "r2_url": urls[1], "split": "held_out", "evidence": evidence})
    rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
    runs_drawn = [r["run"] for r in rows]
    invariant(len(runs_drawn) == len(set(runs_drawn)), "a run is drawn once")
    invariant(all(r["study"] not in existing_studies for r in rows), "the catalogue is append-only")
    return rows


def main(catalogue_path, out_path):
    with open(catalogue_path) as f:
        existing = {r["study"] for r in csv.DictReader(f, delimiter="\t")}
    rows = appended(existing)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    for kit in sorted({r["kit"] for r in rows}):
        n = sum(r["kit"] == kit for r in rows)
        print(f"{kit}: {n} runs from {len({r['study'] for r in rows if r['kit'] == kit})} studies")


if __name__ == "__main__":
    main(*sys.argv[1:])
