# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Shared drawing of a held-out catalogue from ENA, by metadata alone, before any run is read.
Each kit's select.py supplies the ENA query and a classifier that names a run's kit from its own
library protocol (or None to leave the study out); the rules here are common to every kit:

- only paired runs (exactly two FASTQ files) count;
- at most RUNS_PER_STUDY runs per study, the first by run accession;
- optionally a read shape: `shape(r1_length, r2_length)` must accept the modal read lengths of a
  run's two files, read from their first 64 KB only (lengths, never content), and at most
  SHAPE_CHECKS_PER_STUDY runs per study are checked. This drops archive runs whose barcode read was
  lost on export, or companion libraries filed beside the gene-expression runs;
- per kit, the first run of the earliest-published study is the design run; every other row is
  held out (unless the kit's design data came from elsewhere, such as a vendor's demos);
- `names_other_vendor` lets a classifier refuse a protocol that names another vendor's kit too:
  such a study may hold libraries of either (learnt on Parse, 8 Oct 2026);
- `examined_studies`: studies with any run whose read content was looked at while the kit was
  researched. Their rows stay in the catalogue, marked `examined`, and are never held out: runs
  of one study are not independent of what shaped the layout (decided 8 Oct 2026);
- `files`: how many FASTQ files a run may have (2 by default; 10x ATAC's barcode is an index read,
  so its runs come as 3 or 4). Rows record every file in `urls`; for 2-file runs also r1/r2.
"""

import collections
import sys
import os
import csv
import io
import re
import urllib.parse
import urllib.request
import zlib

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

ENA = "https://www.ebi.ac.uk/ena/portal/api/search?"
FIELDS = "run_accession,study_accession,first_public,fastq_ftp,library_construction_protocol"
RUNS_PER_STUDY = 2
# Kit words by vendor, for refusing protocols that name more than one vendor's kit.
VENDORS = {
    "10x": re.compile(r"10x genomics|chromium|gem-x|cell ranger", re.I),
    "bd": re.compile(r"rhapsody", re.I),
    "parse": re.compile(r"evercode|parse bio", re.I),
    "scale": re.compile(r"scalebio|scale bio|quantumscale", re.I),
    "mgi": re.compile(r"dnbelab", re.I),
    "pipseq": re.compile(r"pipseq|pip-seq|fluent bio", re.I),
    "singleron": re.compile(r"singleron|gexscope", re.I),
    "biorad": re.compile(r"ddseq|surecell", re.I),
    "icell8": re.compile(r"i-?cell ?8|wafergen", re.I),
}
SHAPE_CHECKS_PER_STUDY = 6
SHAPE_BYTES = 65536
OUT_FIELDS = ["kit", "study", "run", "r1_url", "r2_url", "split", "evidence", "urls"]
SPLITS = ("design", "held_out", "examined")
# A range fetch that fails is tried again once: a transient failure must not drop a run, or two
# identical draws would differ.
FETCH_ATTEMPTS = 2


def ena_rows(query):
    params = {"result": "read_run", "query": query, "fields": FIELDS, "limit": "0"}
    text = urllib.request.urlopen(ENA + urllib.parse.urlencode(params), timeout=300).read()
    csv.field_size_limit(10**8)
    return list(csv.DictReader(io.StringIO(text.decode()), delimiter="\t"))


def read_length(url):
    """Modal read length in the first SHAPE_BYTES of a gzipped FASTQ; 0 when unreadable."""
    request = urllib.request.Request(url, headers={"Range": f"bytes=0-{SHAPE_BYTES - 1}"})
    text = None
    for _ in range(FETCH_ATTEMPTS):
        try:
            data = urllib.request.urlopen(request, timeout=60).read()
            text = zlib.decompressobj(31).decompress(data).decode(errors="replace")
            break
        except zlib.error:  # not gzip: no retry will change that
            return 0
        except OSError:  # unreachable this time: try again
            continue
    if text is None:  # unreachable on every attempt: a run that cannot pass the shape rule
        return 0
    lines = text.split("\n")
    lengths = collections.Counter(len(lines[i]) for i in range(1, len(lines) - 3, 4))
    return lengths.most_common(1)[0][0] if lengths else 0


def shaped(runs, shape, two_files):
    """The first RUNS_PER_STUDY runs whose read lengths `shape` accepts: `shape(r1, r2)` for
    two-file runs, `shape(names, lengths)` (file names and lengths, in name order) otherwise."""
    kept = []
    for row, files in runs[:SHAPE_CHECKS_PER_STUDY]:
        names = sorted(files)
        lengths = [read_length("https://" + f) for f in names]
        if shape(*lengths) if two_files else shape(names, lengths):
            kept.append((row, files))
        if len(kept) == RUNS_PER_STUDY:
            break
    invariant(len(kept) <= RUNS_PER_STUDY, f"{len(kept)} runs kept from one study")
    return kept


def names_other_vendor(protocol, own):
    """Whether `protocol` names a kit of a vendor other than `own` (a VENDORS key)."""
    invariant(own in VENDORS, f"unknown vendor {own!r}")
    return any(v.search(protocol) for name, v in VENDORS.items() if name != own)


def draw(query, classify, shape=None, design=True, files=(2,), examined_studies=frozenset()):
    """Rows of the catalogue; `classify(protocol)` returns (kit, evidence) or None; `shape`, when
    given, filters runs by their read lengths; `files` is the accepted file counts per run;
    `examined_studies` are marked, never held out (see the module notes)."""
    invariant(files and all(n >= 1 for n in files), f"file counts {files}")
    two_files = tuple(files) == (2,)
    studies = {}
    for row in ena_rows(query):
        run_files = [f for f in row["fastq_ftp"].split(";") if f]
        if len(run_files) in files:
            studies.setdefault(row["study_accession"], []).append((row, run_files))
    rows = []
    for study, runs in studies.items():
        named = classify(runs[0][0]["library_construction_protocol"])
        if named is None:
            continue
        kit, evidence = named
        ordered = sorted(runs, key=lambda rf: rf[0]["run_accession"])
        chosen = shaped(ordered, shape, two_files) if shape else ordered[:RUNS_PER_STUDY]
        for row, run_files in chosen:
            urls = ["https://" + f for f in sorted(run_files)]
            rows.append({"kit": kit, "study": study, "run": row["run_accession"],
                         "first_public": row["first_public"],
                         "r1_url": urls[0] if two_files else "",
                         "r2_url": urls[1] if two_files else "",
                         "urls": ";".join(urls), "evidence": evidence})
    rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
    seen = set()
    for row in rows:
        row["split"] = "design" if design and row["kit"] not in seen else "held_out"
        seen.add(row["kit"])
        if row["split"] == "held_out" and row["study"] in examined_studies:
            row["split"] = "examined"
    drawn_sound(rows, design)
    return rows


def drawn_sound(rows, design):
    """The draw's rules, checked on its output: each run once, two distinct files, the per-study
    cap, and exactly one design run per kit (none when the design data came from elsewhere)."""
    runs = [r["run"] for r in rows]
    invariant(len(runs) == len(set(runs)), "a run is drawn once")
    invariant(all(len(set(r["urls"].split(";"))) == len(r["urls"].split(";")) for r in rows),
              "a run's files differ")
    invariant(all(r["split"] in SPLITS for r in rows), "every row has a known split")
    per_study = collections.Counter(r["study"] for r in rows)
    invariant(max(per_study.values(), default=0) <= RUNS_PER_STUDY, "the per-study cap holds")
    for kit in {r["kit"] for r in rows}:
        designs = sum(r["kit"] == kit and r["split"] == "design" for r in rows)
        invariant(designs == (1 if design else 0), f"{kit}: {designs} design runs")


def examined_studies(path):
    """The studies listed in a kit's results/examined.tsv (run, study, what); none if absent."""
    if not os.path.exists(path):
        return frozenset()
    with open(path) as f:
        rows = list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))
    studies = frozenset(r["study"] for r in rows)
    invariant(all(r["study"] for r in rows), f"{path}: every examined run names its study")
    return studies


def write(rows, out_path):
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, OUT_FIELDS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    for kit in sorted({r["kit"] for r in rows}):
        n = sum(r["kit"] == kit for r in rows)
        studies = len({r["study"] for r in rows if r["kit"] == kit})
        design = sum(r["kit"] == kit and r["split"] == "design" for r in rows)
        held = [r for r in rows if r["kit"] == kit and r["split"] == "held_out"]
        examined = sum(r["kit"] == kit and r["split"] == "examined" for r in rows)
        print(f"{kit}: {n} runs from {studies} studies ({design} design, {len(held)} held out "
              f"from {len({r['study'] for r in held})} studies, {examined} in examined studies)")
