# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws the SeekGene SeekOne catalogue from ENA, by metadata and read lengths alone, before any
held-out run is read.

Rules, fixed before selection (9 Oct 2026):
- query: library protocols naming "seekgene" or "seekone" (ENA wildcards). A protocol must also
  name SeekGene or SeekOne (seek gene, seek one, seeksoul, seekarc, seekspace; any separator).
  Every kit test below reads only the text near a SeekOne, SeekGene or catalogue-number mention
  (40 characters before, 140 after). The method text elsewhere ("3' polyA part", "5' end of the
  primer", "full-length sequencing adapters") is not kit evidence;
- classify, first match wins:
  1. None if the protocol names another vendor's kit (catalogue.VENDORS: 10x, BD, Parse, Scale,
     MGI DNBelab, PIPseq, Singleron, Bio-Rad, ICELL8). None if SeekGene is named only for
     dissociation (SeekMate), counting or service, with no SeekOne kit;
  2. None for SeekSpace, Multiome or ATAC;
  3. MM: the MM kit name ("MM", or "Microwell & Magnetic Beads"), SO01V3 or K00104 ->
     seekgene_mm;
  4. full-length: DD-Q, scFAST, K00801, full-length, full or whole sequence, whole transcriptome,
     FFPE -> seekgene_ddq;
  5. 5': a 5' kit, K00501, immune profiling -> seekgene_dd5;
  6. 3': a 3' kit, K00202 -> seekgene_dd3;
  7. a DD kit (Digital Droplet, DD, Chip S3) with no end stated -> seekgene_dd;
  The evidence is a 120-character window around the match that decided;
- shape: min(r1, r2) >= 28 and max(r1, r2) >= 45 (a barcode read of at least 17 + 12 bases,
  minus one base of slack, and a cDNA read);
- the shared rules in bench/heldout/scripts/catalogue.py (paired runs, at most 2 per study),
  drawn with design=False. The three design runs were read before the draw (9 Oct 2026), so their
  studies are examined: SRR27533148 (PRJNA1064036, DD 3'), SRR23931570 (PRJNA947149, DD 5') and
  SRR22784277 (PRJNA913187, MM). After the draw, exactly these three runs are marked design. Any
  other row from their studies is examined. A design run that the draw did not keep is added as a
  design row with its URLs.

Revisions, 9 Oct 2026, after the first draft draw and before any run of the catalogue was read
or run through the tool (only the three design runs above had been read):
1. mixed: a protocol that names both an MM kit and a DD kit (Digital Droplet or DD beside a 3' or
   5' kit) is refused. The study may hold libraries of either kit. PRJNA1224797 names the MM
   Single Cell 3' kit for its scRNA and the Digital Droplet Single Cell 3' kit for its snRNA, and
   rule 3 had labelled both its runs MM;
2. shape: a run is also accepted when R1 >= 29 (a full B17U12 barcode read) and R2 >= 40. The
   rule above refused PRJEB93881, a DD 3' study ("SeekOne DD Single Cell 3'") whose runs are
   R1 29 and R2 43;
3. evidence: the 120-character window is cut from the protocol itself, around the deciding match,
   not from the joined mention windows. 14 windows of the first draft ran across a " || " join.
   Kit tests still read only the text inside the mention spans.

Only "SeekGene SeekOne DD 3'" (DDV2) and "SeekGene SeekOne DD 5'" (DD5V1) are claimed by the tool.
The other kits are catalogued so their runs are known. Their calls are not scored.

Usage: python select.py <catalogue.tsv> [<refused.tsv>]
(Outside the repository, set TUSSL_BENCH to the repository's bench/ directory.)
"""

import csv
import os
import re
import sys

BENCH = os.environ.get("TUSSL_BENCH") or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       "..", "..")
sys.path.insert(0, os.path.join(BENCH, "heldout", "scripts"))
sys.path.insert(0, BENCH)
import catalogue  # noqa: E402  (shared ENA draw: paired runs, per-study cap)
from invariant import invariant  # noqa: E402

QUERY = ('library_construction_protocol="*seekgene*" OR '
         'library_construction_protocol="*seekone*"')
SEEKGENE = re.compile(r"seek\W*gene|seek\W*one|seek\W*soul|seekarc|seek\W*space", re.I)
MENTION = re.compile(r"seek\W*one|seek\W*gene|seek\W*space|K0\d{4}|SO01V3", re.I)
HYPHENS = re.compile("[‐-―−]")
Q = "['′’ʹ]"
SEEKONE = re.compile(r"seek\W*one", re.I)
EXCLUDED = re.compile(r"seek\s*space|multi-?omics?|multiome|ATAC", re.I)
MM = re.compile(r"\bMM\b|SO01V3|K00104")
MM_NAME = re.compile(r"microwell\s*(?:&|and)\s*magnetic", re.I)
FULL = re.compile(r"DD-Q|scFAST|K00801|full[- ]?length|full[- ]sequence|whole[- ]sequence|"
                  r"whole[- ]tran?s?criptome|FFPE", re.I)
FIVE = re.compile(r"single[- ]?cell\s*5\s*" + Q + r"|5\s*" + Q + r"\s*(?:librar|kit|transcriptome|"
                  r"gene)|K00501|immune profiling", re.I)
THREE = re.compile(r"single[- ]?cell\s*3\s*" + Q + r"|3\s*" + Q + r"\s*(?:librar|kit|transcriptome|"
                   r"gene)|K00202", re.I)
DD = re.compile(r"digital droplet|\bDD\b|chip\s*S3", re.I)
# Revision 1: a DD kit named beside its end, as "Digital Droplet Single Cell 3' ... Kit".
DD_KIT = re.compile(r"(?:digital droplet|\bDD\b)[^.;]{0,40}?(?:single[- ]?cell\s*[35]\s*" + Q +
                    r"|[35]\s*" + Q + r"\s*(?:librar|kit|transcriptome|gene))", re.I)
DESIGN = {"SRR27533148": "PRJNA1064036", "SRR23931570": "PRJNA947149",
          "SRR22784277": "PRJNA913187"}
EXAMINED_STUDIES = frozenset(DESIGN.values())
WHY_EXAMINED = "design run read to confirm the layout, 9 Oct 2026"

refusals = {}  # protocol -> reason, filled by classify for the refusal report


def spans(protocol):
    """The (start, end) spans near each SeekOne, SeekGene or catalogue-number mention."""
    return [(max(0, m.start() - 40), m.end() + 140) for m in MENTION.finditer(protocol)]


def context(protocol):
    """The text near each SeekOne, SeekGene or catalogue-number mention."""
    return " || ".join(protocol[a:b] for a, b in spans(protocol))


def near(pattern, protocol):
    """The first match of `pattern` in the protocol itself that lies inside a mention span."""
    for a, b in spans(protocol):
        m = pattern.search(protocol, a, min(b, len(protocol)))
        if m:
            return m
    return None


def window(text, match):
    """A 120-character evidence window cut from the protocol, centred on `match` (revision 3)."""
    middle = (match.start() + match.end()) // 2
    start = max(0, middle - 60)
    return " ".join(text[start:start + 120].split())


def refuse(protocol, reason):
    refusals[protocol] = reason
    return None


def classify(protocol):
    p = HYPHENS.sub("-", protocol)
    if not SEEKGENE.search(p):
        return refuse(protocol, "no SeekGene or SeekOne name")
    if any(v.search(p) for v in catalogue.VENDORS.values()):
        return refuse(protocol, "names another vendor's kit")
    c = context(p)
    if EXCLUDED.search(c):
        return refuse(protocol, "SeekSpace, Multiome or ATAC")
    if (near(MM, p) or near(MM_NAME, p)) and near(DD_KIT, p):
        return refuse(protocol, "mixed: names both an MM kit and a DD kit")  # revision 1
    for kit, pattern in (("seekgene_mm", MM), ("seekgene_mm", MM_NAME),
                         ("seekgene_ddq", FULL), ("seekgene_dd5", FIVE),
                         ("seekgene_dd3", THREE)):
        m = near(pattern, p)
        if m:
            return kit, window(p, m)
    m = near(DD, p)
    if m and SEEKONE.search(c):
        return "seekgene_dd", window(p, m)
    if SEEKONE.search(c):
        return refuse(protocol, "SeekOne named, no kit")
    return refuse(protocol, "SeekGene named for dissociation, counting or service only")


def gene_expression_shape(r1_length, r2_length):
    """The 9 Oct rule, or (revision 2) a full 29-base barcode read with a cDNA read of 40+."""
    if min(r1_length, r2_length) >= 28 and max(r1_length, r2_length) >= 45:
        return True
    return r1_length >= 29 and r2_length >= 40


def mark_design(rows):
    """The three runs read before the draw are design; other rows of their studies examined."""
    by_run = {r["run"]: r for r in rows}
    for run, study in DESIGN.items():
        if run in by_run:
            invariant(by_run[run]["study"] == study, f"{run} drawn under another study")
            by_run[run]["split"] = "design"
        else:
            rows.append(design_row(run, study))
    for row in rows:
        if row["study"] in EXAMINED_STUDIES and row["run"] not in DESIGN:
            invariant(row["split"] == "examined", f"{row['run']}: {row['split']}")
    designs = [r["run"] for r in rows if r["split"] == "design"]
    invariant(sorted(designs) == sorted(DESIGN), f"design runs {designs}")
    rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
    return rows


def design_row(run, study):
    """A design run the draw did not keep, with its URLs, kit and evidence from its protocol."""
    found = [r for r in catalogue.ena_rows(f"run_accession={run}")]
    invariant(len(found) == 1, f"{run}: {len(found)} ENA rows")
    row = found[0]
    named = classify(row["library_construction_protocol"])
    invariant(named is not None, f"{run}: its protocol names no kit")
    urls = ["https://" + f for f in sorted(f for f in row["fastq_ftp"].split(";") if f)]
    invariant(len(urls) == 2, f"{run}: {len(urls)} files")
    return {"kit": named[0], "study": study, "run": run, "first_public": row["first_public"],
            "r1_url": urls[0], "r2_url": urls[1], "urls": ";".join(urls),
            "evidence": named[1], "split": "design"}


def draw():
    rows = catalogue.draw(QUERY, classify, gene_expression_shape, design=False,
                          examined_studies=EXAMINED_STUDIES)
    return mark_design(rows)


def refused(rows):
    """One row per study left out: by classify, or classified with no run kept."""
    drawn = {r["study"] for r in rows}
    out, seen = [], set()
    ena = catalogue.ena_rows(QUERY)
    paired = {}
    for row in ena:
        two = len([f for f in row["fastq_ftp"].split(";") if f]) == 2
        paired[row["study_accession"]] = paired.get(row["study_accession"], 0) + two
    for row in ena:
        study = row["study_accession"]
        if study in seen or study in drawn:
            continue
        seen.add(study)
        protocol = row["library_construction_protocol"]
        named = classify(protocol)
        if named is None:
            reason = refusals[protocol]
        elif paired[study] == 0:
            reason = f"classified {named[0]}, but no paired run (exactly two FASTQ files)"
        else:
            reason = (f"classified {named[0]}, but none of its first "
                      f"{catalogue.SHAPE_CHECKS_PER_STUDY} paired runs passed the shape rule")
        out.append({"study": study, "first_run": row["run_accession"], "reason": reason,
                    "protocol_start": " ".join(protocol[:200].split())})
    return out


def main(out_path, refused_path=None):
    rows = draw()
    catalogue.write(rows, out_path)
    if refused_path:
        with open(refused_path, "w", newline="") as f:
            writer = csv.DictWriter(f, ["study", "first_run", "reason", "protocol_start"],
                                    delimiter="\t")
            writer.writeheader()
            writer.writerows(refused(rows))


if __name__ == "__main__":
    main(*sys.argv[1:])
