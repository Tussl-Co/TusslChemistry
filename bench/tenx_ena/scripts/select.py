# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draws a held-out catalogue of 10x Genomics single-cell RNA runs from all of ENA, labelled from
each run's own library protocol text, by metadata and read lengths alone, before any run is read.

Rules, fixed before the draw (9 Oct 2026):
- query: ENA portal, result=read_run, runs whose library_construction_protocol matches any of
  "*10x*", "*chromium*", "*gemx*", "*multiome*", "*single?cell?3*", "*single?cell?5*",
  "*singlecell3*", "*singlecell5*" (ENA wildcards; ENA matches nothing across a space, hyphen or
  apostrophe, so "*gem-x*" and "*single cell 3'*" return 0 rows and are covered by the others).
  Fetched as one query per first_public year (work/fetch.py, cached in raw/; the partition counts
  must add up to the union count). Every run of a study returned by the query is read;
- label: every distinct protocol is split into sentences as the committed
  head_to_head/scripts/lab_labels/build.py does (`sentences`), each sentence passing
  `label10x.is_candidate` (or naming a 10x guide number) is mapped by the committed `label10x.py`,
  imported unchanged. A run's label is its one exact layout; two exact layouts, or any sentence
  `label10x` calls "ambiguous", is "conflict" (stricter than build.py, which skips ambiguous
  sentences: a sentence naming two kits shows the study may hold two assays); only partial
  statements of one end give "end_only (3')" / "end_only (5')"; none give "not stated";
- kept labels and kit keys: "10x 3' v2" tenx_3p_v2, "10x 3' v3/v3.1" tenx_3p_v3,
  "10x 3' v4 (GEM-X)" tenx_3p_v4, "10x 5' v1/v2" tenx_5p_v1v2, "10x 5' v3 (GEM-X)" tenx_5p_v3,
  "10x Multiome GEX" tenx_multiome. Every other label (not stated, end only, conflict, Flex v1,
  Flex v2, 3' v1) refuses the study;
- refusals, per study, over every distinct protocol of its runs, first match wins:
  1. used: the study (PRJ), or any of its runs, experiments or secondary study accessions, is
     named in a used source: every bench/*/results/catalogue*.tsv (all splits), the
     head_to_head result tables (studies.tsv, lab_label_runs.tsv and the others, except
     results.tsv, which is never opened), arc/results/selection.tsv, tenx_lab's catalogue.tsv,
     labels.tsv and excluded.tsv, bench/long_read, every examined*.tsv, and
     exec2/examined_wild.tsv. Run and experiment accessions and GSE series are mapped to their
     ENA study (work/map_used.py);
  2. another vendor: `catalogue.names_other_vendor(protocol, "10x")`;
  3. another assay: a protocol naming Visium or spatial (visium, spatial, space ranger); Flex or
     fixed RNA profiling (label10x FLEX1 or FLEX2); ATAC without Multiome (an ATAC library with
     no Multiome or ARC word: ATAC-only, or ATAC beside a separate RNA kit); V(D)J, TCR or BCR
     libraries with no gene-expression word (label10x GEX_WORD);
  4. labels: the study's runs carry more than one label. "not stated" runs (hashtag, ADT or
     other libraries whose protocol names no kit) do not count here: they are never drawn. Every
     other label counts, end-only and conflict included;
  5. label: the study's one label is not kept;
  6. no paired run: no run with the study's label has exactly two FASTQ files;
  7. shape: none of the first SHAPE_CHECKS_PER_STUDY (6) paired runs passes the shape rule;
- draw: the shared rules of bench/heldout/scripts/catalogue.py: only paired runs whose own
  label is the study's; runs in accession order; at most SHAPE_CHECKS_PER_STUDY checked and the
  first RUNS_PER_STUDY (2) that pass kept; `drawn_sound(rows, design=False)`; every row
  held_out (these layouts have their design runs already);
- shape: modal read lengths (catalogue.read_length, first 64 KB) of the two files in name order:
  R1 >= 26 and R2 >= 50, or R1 >= 50 and R2 >= 26 (the barcode read filed second);
- cap: at most CAP (1,500) studies. The frame is every study that passes rules 1 to 6. If the
  frame is larger than CAP, its accessions are sorted, shuffled with random.Random(SEED)
  (SEED 20261009), and walked in that order; each is shape-checked, and the walk stops when CAP
  studies have a kept run. Studies after the stop are not checked, not refused, and are listed
  in work/not_sampled.tsv;
- pacing: requests to EBI (ENA portal and the FASTQ range reads) about 3 per second; ENA calls
  are retried once; read_length tries twice itself. Every failure is in fetch_failures.tsv;
- evidence: a 120-character window cut from the protocol itself (tab, CR and LF replaced by one
  space each, so lengths are kept), around the match that decided the label (version, GEM-X,
  Multiome, then end word) in the first sentence giving the run's label. Windows are tried
  centred, then shifted; the first whose `label10x.label` is the run's label is kept.

Batch 2 (10 Oct 2026), rules fixed before its walk; batch 1 above is unchanged and is reproduced
by `select.py <dir>` from its cached raw/ and work/:
- the order is batch 1's: the batch-1 frame (from <dir>'s cache) sorted and shuffled with
  random.Random(SEED). Batch 1 walked positions 0 to 1,517 (BATCH1_WALKED = 1,518); batch 2
  continues from position 1,518 and shape-checks in order until BATCH2_STUDIES (2,100) more
  studies have a kept run. The size was 1,100 when the walk began; it was raised to 2,100
  before any run, after disk was freed, 10 Oct 2026. The walk is one batch. Every modal length
  measured is appended at once to <dir2>/work/lengths_cache.tsv; an interrupted pass is resumed
  by re-running the same walk with that file as a cache (`--batch2 <dir2> <cache>`; a cached
  length of 0 is fetched again), so the studies already walked are decided exactly as before.
  (A first pass, begun at size 1,100 and stopped by a session interruption at about 900 drawn
  studies, kept no lengths; the walk was re-run from position 1,518 with this cache.)
- the frame is re-derived with today's ENA (the same eight terms, fetched again into
  <dir2>/raw/) and today's used set: the same sources, which now include
  bench/tenx_ena/results/catalogue.tsv (batch 1) and every other bench catalogue at the
  repository's HEAD (3759832), plus PRJNA1103675 and PRJNA1192233, the studies of 11 tenx_lab
  SRX that the ENA portal could not map and NCBI runinfo did (batch 1's post-draw check). The
  order is not reshuffled. A study of the order that is no
  longer in today's frame (refused today by rules 1 to 6, or not returned), or whose label
  changed, is refused as "changed since batch 1" and recorded in work/changes.tsv; a study whose
  paired runs changed is walked with today's runs and recorded. Studies new to today's frame are
  not in the order: they are listed in work/new_in_frame.tsv and not walked;
- the same rules, shape, 2-run cap and evidence windows; every row held_out; no batch-1 study
  or run again (checked);
- writes <dir2>/appended.tsv (new rows only, catalogue.OUT_FIELDS), refused.tsv (walked studies
  that failed, with their position), not_sampled.tsv (positions left after the batch, with
  today's state), fetch_failures.tsv, and work/.

Usage: python -I -B select.py <dir>            (batch 1)
       python -I -B select.py <dir> --batch2 <dir2> [<lengths cache>]   (batch 2)
(Batch 1: dir holds raw/ena_*.tsv.gz and work/ from work/fetch.py and
work/map_used.py; `--frame` stops before any read-length fetch; batch 1 writes catalogue.tsv, refused.tsv,
fetch_failures.tsv and work/ files.)
Outside the repository, set TUSSL_BENCH to the repository's bench/ directory.
"""

import collections
import csv
import glob
import gzip
import json
import os
import random
import re
import sys
import time

BENCH = os.environ.get("TUSSL_BENCH") or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                       "..", "..")
sys.path.insert(0, os.path.join(BENCH, "heldout", "scripts"))
sys.path.insert(0, os.path.join(BENCH, "head_to_head", "scripts", "lab_labels"))
sys.path.insert(0, BENCH)
import catalogue  # noqa: E402  (shared rules: per-study cap, read_length, write)
import label10x as L  # noqa: E402  (committed mapper, unchanged)
from invariant import invariant  # noqa: E402

KEEP = {"10x 3' v2": "tenx_3p_v2", "10x 3' v3/v3.1": "tenx_3p_v3",
        "10x 3' v4 (GEM-X)": "tenx_3p_v4", "10x 5' v1/v2": "tenx_5p_v1v2",
        "10x 5' v3 (GEM-X)": "tenx_5p_v3", "10x Multiome GEX": "tenx_multiome"}
CAP = 1500
SEED = 20261009
EVID = 120
PACE = 0.34
SPATIAL = re.compile(r"visium|spatial|space\s?ranger", re.I)
ATAC = re.compile(r"atac", re.I)
VDJ = re.compile(r"v\(d\)j|\bvdj\b|\btcr\b|\bbcr\b|t[- ]cell receptor|b[- ]cell receptor|"
                 r"immune repertoire", re.I)
WS = str.maketrans({"\t": " ", "\r": " ", "\n": " "})
csv.field_size_limit(10**9)


def sentences(text):
    """As build.py: lines, then split after '.' or ';'."""
    out = []
    for line in text.splitlines():
        out += [s.strip() for s in re.split(r"(?<=[.;])\s+", line) if s.strip()]
    return out


def run_label(protocol):
    """(label, deciding sentence or None), by build.py's aggregation (ambiguous -> conflict)."""
    exact, partial, first = [], [], {}
    ambiguous = False
    for s in sentences(protocol):
        if not (L.is_candidate(s) or L.guide_numbers(s)):
            continue
        m = L.label(s)
        if m == "ambiguous":
            ambiguous = True
        elif m in L.EXACT:
            exact.append(m)
            first.setdefault(m, s)
        elif m in L.PARTIAL:
            partial.append(m)
    kinds, pk = set(exact), set(partial)
    if ambiguous or len(kinds) > 1:
        return "conflict", None
    if kinds:
        lab = kinds.pop()
        return lab, first[lab]
    if len(pk) == 1:
        return "end_only (" + pk.pop().split()[0] + ")", None
    return "not stated", None


def other_assay(protocol):
    if SPATIAL.search(protocol):
        return "another assay: Visium or spatial"
    if L.FLEX1.search(protocol) or L.FLEX2.search(protocol):
        return "another assay: Flex or fixed RNA profiling"
    if ATAC.search(protocol) and not L.MULTI.search(protocol):
        return "another assay: ATAC without Multiome"
    if VDJ.search(protocol) and not L.GEX_WORD.search(protocol):
        return "another assay: V(D)J only, no gene expression"
    return None


def evidence(protocol, sentence, lab):
    """A 120-character window of the protocol (tab/CR/LF as spaces) around the deciding match."""
    text = protocol.translate(WS)
    base = protocol.find(sentence)
    invariant(base >= 0, "the deciding sentence is in the protocol")
    if len(text) <= EVID:
        return text, True
    hits = []
    for rx in (L.VER, L.GEMX, L.MULTI, L.THREE, L.FIVE, L.CONTEXT):
        hits += [base + (m.start() + m.end()) // 2 for m in rx.finditer(sentence)]
    hits = hits or [base]
    tried = []
    for shift in (60, 40, 80, 20, 100, 0, 110):
        for h in hits:
            a = max(0, min(h - shift, len(text) - EVID))
            if a not in tried:
                tried.append(a)
                if L.label(text[a:a + EVID]) == lab:
                    return text[a:a + EVID], True
    a = max(0, min(hits[0] - 60, len(text) - EVID))
    return text[a:a + EVID], False


class Pacer:
    def __init__(self):
        self.last = 0.0

    def wait(self):
        w = PACE - (time.time() - self.last)
        if w > 0:
            time.sleep(w)
        self.last = time.time()


def load(d):
    """ENA rows of every partition, by study; protocols interned."""
    runs, per_file = {}, {}
    for p in sorted(glob.glob(f"{d}/raw/ena_*.tsv.gz")):
        with gzip.open(p, "rt", encoding="utf-8") as f:
            # ENA's TSV is unquoted: a '"' in a protocol is text, not a quote
            rows = list(csv.DictReader(f, delimiter="\t", quoting=csv.QUOTE_NONE))
        invariant(all(None not in r and None not in r.values() for r in rows), f"{p}: ragged row")
        per_file[os.path.abspath(p)] = len(rows)
        for r in rows:
            runs[r["run_accession"]] = r
    by_study = collections.defaultdict(list)
    for r in runs.values():
        by_study[r["study_accession"]].append(r)
    return runs, by_study, per_file


def used_sets(d):
    studies = json.load(open(f"{d}/work/used_studies.json"))
    tokens = json.load(open(f"{d}/work/used_tokens.json"))
    return studies, tokens


def shape_ok(r1, r2):
    return (r1 >= 26 and r2 >= 50) or (r1 >= 50 and r2 >= 26)


def frame_of(d):
    """Every study of the cached query in `d`, refused by rules 1 to 6 or placed in the frame."""
    queries = json.load(open(f"{d}/work/queries.json"))
    union = next(q["count"] for q in queries if q["kind"] == "union")
    parts = [q for q in queries if q["kind"] == "partition"]
    runs, by_study, per_file = load(d)
    for q in parts:
        q["rows"] = per_file[os.path.abspath(os.path.join(d, "raw", os.path.basename(q["path"])))]
    invariant(all(q["rows"] == q["count"] for q in parts), "each partition holds its count")
    invariant(sum(q["count"] for q in parts) == union, "the partitions add up to the union")
    invariant(len(runs) == union, f"{len(runs)} distinct runs, union {union}")
    used, tokens = used_sets(d)
    labels = {}  # protocol -> (label, sentence)
    refused, frame, info = [], [], {}
    for study in sorted(by_study):
        rr = sorted(by_study[study], key=lambda r: r["run_accession"])
        protocols = list(dict.fromkeys(r["library_construction_protocol"] for r in rr))
        for p in protocols:
            if p not in labels:
                labels[p] = run_label(p)
        run_labs = collections.Counter(labels[r["library_construction_protocol"]][0] for r in rr)
        stated = sorted(k for k in run_labs if k != "not stated")
        names = {study, rr[0]["secondary_study_accession"]} | {r["run_accession"] for r in rr} \
            | {r["experiment_accession"] for r in rr}
        hit = sorted(n for n in names if n in used or n in tokens)
        reason = None
        if hit:
            src = sorted({s for n in hit for s in (used.get(n) or tokens.get(n))})
            reason = "used: " + ",".join(src)
        elif any(catalogue.names_other_vendor(p, "10x") for p in protocols):
            reason = "another vendor named"
        else:
            reason = next((x for x in map(other_assay, protocols) if x), None)
        if reason is None and len(stated) > 1:
            reason = "runs carry more than one label: " + " | ".join(stated)
        if reason is None:
            lab = stated[0] if stated else "not stated"
            if lab not in KEEP:
                reason = f"label: {lab}"
        if reason is None:
            paired = [r for r in rr if labels[r["library_construction_protocol"]][0] == lab
                      and len([f for f in r["fastq_ftp"].split(";") if f]) == 2]
            if not paired:
                reason = "no paired run (exactly two FASTQ files) with the study's label"
        info[study] = {"runs": len(rr), "labels": dict(run_labs),
                       "not_stated_alongside": bool(len(stated) == 1 and "not stated" in run_labs)}
        if reason:
            refused.append({"study": study, "reason": reason, "labels": json.dumps(dict(run_labs)),
                            "first_run": rr[0]["run_accession"],
                            "protocol_start": rr[0]["library_construction_protocol"][:200].translate(WS)})
        else:
            frame.append((study, lab, paired))
    return runs, by_study, labels, refused, frame, info


def walk_order(frame):
    """The batch-1 order: the frame's accessions sorted, then shuffled with SEED."""
    order = sorted(frame)
    if len(order) > CAP:
        random.Random(SEED).shuffle(order)
    return order


def main(d, frame_only=False):
    runs, by_study, labels, refused, frame, info = frame_of(d)
    if frame_only:
        print(len(frame), "frame studies;", collections.Counter(lab for _, lab, _ in frame))
        print(collections.Counter(r["reason"].split(":")[0] for r in refused))
        return
    # cap and sampling
    order = walk_order(frame)
    sampled = len(order) > CAP
    pacer = Pacer()
    rows, failures, shape_log, not_sampled, moved_lost = [], [], [], [], 0
    drawn = 0
    for i, (study, lab, paired) in enumerate(order):
        if drawn == CAP:
            not_sampled.append((study, lab))
            continue
        kept = []
        for r in paired[:catalogue.SHAPE_CHECKS_PER_STUDY]:
            files = sorted(f for f in r["fastq_ftp"].split(";") if f)
            lengths = []
            for f in files:
                pacer.wait()
                n = catalogue.read_length("https://" + f)
                if n == 0:
                    failures.append(("read_length", r["run_accession"], "https://" + f,
                                     "0: unreachable on 2 attempts, or not gzip"))
                lengths.append(n)
            ok = shape_ok(*lengths)
            shape_log.append((r["run_accession"], study, *lengths, "pass" if ok else "fail"))
            if ok:
                kept.append((r, files))
            if len(kept) == catalogue.RUNS_PER_STUDY:
                break
        if not kept:
            refused.append({"study": study, "reason": "shape: no paired run passed (first "
                            f"{catalogue.SHAPE_CHECKS_PER_STUDY} checked)",
                            "labels": json.dumps(info[study]["labels"]),
                            "first_run": paired[0]["run_accession"],
                            "protocol_start": paired[0]["library_construction_protocol"][:200].translate(WS)})
            continue
        drawn += 1
        for r, files in kept:
            p = r["library_construction_protocol"]
            ev, ok = evidence(p, labels[p][1], lab)
            moved_lost += not ok
            urls = ["https://" + f for f in files]
            rows.append({"kit": KEEP[lab], "study": study, "run": r["run_accession"],
                         "first_public": r["first_public"], "r1_url": urls[0], "r2_url": urls[1],
                         "urls": ";".join(urls), "evidence": ev, "split": "held_out"})
        if i % 50 == 0:
            print(f"{i} walked, {drawn} drawn", flush=True)
    rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
    catalogue.drawn_sound(rows, False)
    invariant(all(r["split"] == "held_out" for r in rows), "every row held out")
    invariant(len({r["study"] for r in rows}) <= CAP, "the cap holds")
    catalogue.write(rows, f"{d}/catalogue.tsv")
    refused.sort(key=lambda r: r["study"])
    with open(f"{d}/refused.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, ["study", "reason", "labels", "first_run", "protocol_start"],
                           delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(refused)
    for name in ("fetch_fail_ena.json",):
        for t in json.load(open(f"{d}/work/{name}")):
            failures.append(tuple(t) + ("",) * (4 - len(t)))
    for a in json.load(open(f"{d}/work/fetch_fail_map.json"))["fails"]:
        failures.append(tuple(a) + ("",) * (4 - len(a)))
    with open(f"{d}/fetch_failures.tsv", "w", encoding="utf-8") as f:
        f.write("source\taccession\turl\terror\n")
        f.writelines("\t".join(map(str, t[:4])) + "\n" for t in failures)
    with open(f"{d}/work/shape_log.tsv", "w", encoding="utf-8") as f:
        f.write("run\tstudy\tlen_file1\tlen_file2\tshape\n")
        f.writelines("\t".join(map(str, t)) + "\n" for t in shape_log)
    with open(f"{d}/work/not_sampled.tsv", "w", encoding="utf-8") as f:
        f.write("study\tlabel\n")
        f.writelines(f"{s}\t{lab}\n" for s, lab in not_sampled)
    label_count = collections.Counter(v[0] for v in labels.values())
    summary = {"distinct_runs": len(runs), "studies": len(by_study), "frame_studies": len(frame),
               "sampled": sampled, "seed": SEED if sampled else None,
               "walked": min(len(order), len(order) - len(not_sampled)),
               "not_sampled": len(not_sampled), "drawn_studies": drawn, "rows": len(rows),
               "shape_checked": len(shape_log),
               "shape_failed": sum(t[-1] == "fail" for t in shape_log),
               "evidence_label_not_reproduced": moved_lost,
               "distinct_protocols": len(labels), "protocol_labels": dict(label_count),
               "studies_with_not_stated_alongside_kept_label": sum(
                   1 for s, _, _ in frame if info[s]["not_stated_alongside"]),
               "frame_by_label": dict(collections.Counter(lab for _, lab, _ in frame)),
               "failures": dict(collections.Counter(t[0] for t in failures))}
    json.dump(summary, open(f"{d}/work/summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


# Batch 2 (10 Oct 2026): see the module docstring.
BATCH1_WALKED = 1518   # positions 0..1517 of the order were walked by batch 1
BATCH2_STUDIES = 2100  # raised from 1,100 before any run, after disk was freed, 10 Oct 2026


def batch2(d1, d2, lengths_cache=None):
    """Continue batch 1's order from BATCH1_WALKED with today's ENA (cached in d2) and today's used set."""
    frame1 = frame_of(d1)[4]
    order = walk_order(frame1)
    cat1 = list(csv.DictReader(open(f"{d1}/catalogue.tsv", encoding="utf-8"), delimiter="\t"))
    drawn1 = {r["study"] for r in cat1}
    walked1 = {s for s, _, _ in order[:BATCH1_WALKED]}
    invariant(drawn1 <= walked1 and len(drawn1) == CAP, "batch 1 is the first 1,518 of the order")
    invariant(order[BATCH1_WALKED - 1][0] in drawn1, "batch 1 stopped on a drawn study")
    runs, by_study, labels, refused_today, frame2, info = frame_of(d2)
    today = {s: (lab, paired) for s, lab, paired in frame2}
    why = {r["study"]: r["reason"] for r in refused_today}
    in_order = {s for s, _, _ in order}
    cached = {}  # run -> modal lengths measured by an earlier, interrupted pass of this same walk
    if lengths_cache and os.path.exists(lengths_cache):
        for r in csv.DictReader(open(lengths_cache, encoding="utf-8"), delimiter="\t"):
            if int(r["len_file1"]) > 0 and int(r["len_file2"]) > 0:
                cached[r["run"]] = [int(r["len_file1"]), int(r["len_file2"])]
    reused = 0
    # every fresh measurement is appended at once, so an interrupted pass loses nothing
    cache_out = open(f"{d2}/work/lengths_cache.tsv", "a", encoding="utf-8")
    if cache_out.tell() == 0:
        cache_out.write("run\tstudy\tlen_file1\tlen_file2\n")
    pacer = Pacer()
    rows, refused, changes, failures, shape_log, not_sampled = [], [], [], [], [], []
    drawn, first, last = 0, None, None
    for pos in range(BATCH1_WALKED, len(order)):
        study, lab1, paired1 = order[pos]
        if drawn == BATCH2_STUDIES:
            state = "in today's frame" if study in today else why.get(study, "not returned today")
            not_sampled.append((pos, study, lab1, state))
            continue
        first = pos if first is None else first
        last = pos
        if study not in today:
            reason = why.get(study, "not returned by today's query")
            changes.append((pos, study, "left the frame", reason))
            refused.append({"position": pos, "study": study, "reason": "changed since batch 1: " + reason,
                            "labels": json.dumps(info.get(study, {}).get("labels", {}))})
            continue
        lab, paired = today[study]
        if lab != lab1:
            changes.append((pos, study, "label", f"{lab1} -> {lab}"))
            refused.append({"position": pos, "study": study,
                            "reason": f"changed since batch 1: label {lab1} -> {lab}",
                            "labels": json.dumps(info[study]["labels"])})
            continue
        r1 = [r["run_accession"] for r in paired1]
        r2 = [r["run_accession"] for r in paired]
        if r1 != r2:
            changes.append((pos, study, "paired runs", f"{len(r1)} -> {len(r2)} (today's used)"))
        kept = []
        for r in paired[:catalogue.SHAPE_CHECKS_PER_STUDY]:
            files = sorted(f for f in r["fastq_ftp"].split(";") if f)
            lengths = cached.get(r["run_accession"])
            reused += lengths is not None
            if lengths is None:
                lengths = []
                for f in files:
                    pacer.wait()
                    n = catalogue.read_length("https://" + f)
                    if n == 0:
                        failures.append(("read_length", r["run_accession"], "https://" + f,
                                         "0: unreachable on 2 attempts, or not gzip"))
                    lengths.append(n)
                cache_out.write(f"{r['run_accession']}\t{study}\t{lengths[0]}\t{lengths[1]}\n")
                cache_out.flush()
            ok = shape_ok(*lengths)
            shape_log.append((r["run_accession"], study, *lengths, "pass" if ok else "fail"))
            if ok:
                kept.append((r, files))
            if len(kept) == catalogue.RUNS_PER_STUDY:
                break
        if not kept:
            refused.append({"position": pos, "study": study, "reason": "shape: no paired run passed "
                            f"(first {catalogue.SHAPE_CHECKS_PER_STUDY} checked)",
                            "labels": json.dumps(info[study]["labels"])})
            continue
        drawn += 1
        for r, files in kept:
            p = r["library_construction_protocol"]
            ev, _ = evidence(p, labels[p][1], lab)
            urls = ["https://" + f for f in files]
            rows.append({"kit": KEEP[lab], "study": study, "run": r["run_accession"],
                         "first_public": r["first_public"], "r1_url": urls[0], "r2_url": urls[1],
                         "urls": ";".join(urls), "evidence": ev, "split": "held_out"})
        if drawn % 50 == 0:
            print(f"position {pos}, {drawn} drawn", flush=True)
    rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
    catalogue.drawn_sound(rows, False)
    invariant(not ({r["study"] for r in rows} & drawn1), "no batch-1 study again")
    invariant(not ({r["run"] for r in rows} & {r["run"] for r in cat1}), "no batch-1 run again")
    catalogue.drawn_sound(cat1 + rows, False)
    catalogue.write(rows, f"{d2}/appended.tsv")
    with open(f"{d2}/refused.tsv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, ["position", "study", "reason", "labels"], delimiter="\t",
                           lineterminator="\n")
        w.writeheader()
        w.writerows(refused)
    with open(f"{d2}/not_sampled.tsv", "w", encoding="utf-8") as f:
        f.write("position\tstudy\tlabel\ttoday\n")
        f.writelines("\t".join(map(str, t)) + "\n" for t in not_sampled)
    with open(f"{d2}/work/changes.tsv", "w", encoding="utf-8") as f:
        f.write("position\tstudy\tchange\tdetail\n")
        f.writelines("\t".join(map(str, t)) + "\n" for t in changes)
    new = sorted(s for s in today if s not in in_order)
    with open(f"{d2}/work/new_in_frame.tsv", "w", encoding="utf-8") as f:
        f.write("study\tlabel\n")
        f.writelines(f"{s}\t{today[s][0]}\n" for s in new)
    with open(f"{d2}/work/shape_log.tsv", "w", encoding="utf-8") as f:
        f.write("run\tstudy\tlen_file1\tlen_file2\tshape\n")
        f.writelines("\t".join(map(str, t)) + "\n" for t in shape_log)
    for t in json.load(open(f"{d2}/work/fetch_fail_ena.json")):
        failures.append(tuple(t) + ("",) * (4 - len(t)))
    for a in json.load(open(f"{d2}/work/fetch_fail_map.json"))["fails"]:
        failures.append(tuple(a) + ("",) * (4 - len(a)))
    with open(f"{d2}/fetch_failures.tsv", "w", encoding="utf-8") as f:
        f.write("source\taccession\turl\terror\n")
        f.writelines("\t".join(map(str, t[:4])) + "\n" for t in failures)
    summary = {"order_length": len(order), "first_position": first, "last_position": last,
               "walked": 0 if first is None else last - first + 1, "drawn_studies": drawn,
               "rows": len(rows), "refused": len(refused), "not_sampled": len(not_sampled),
               "changes": dict(collections.Counter(c[2] for c in changes)),
               "today_frame": len(today), "new_in_frame_not_walked": len(new),
               "today_distinct_runs": len(runs), "shape_checked": len(shape_log),
               "lengths_reused_from_cache": reused, "batch2_studies": BATCH2_STUDIES,
               "shape_failed": sum(t[-1] == "fail" for t in shape_log),
               "failures": dict(collections.Counter(t[0] for t in failures))}
    json.dump(summary, open(f"{d2}/work/summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))


if __name__ == "__main__":
    if sys.argv[2:3] == ["--batch2"]:
        batch2(sys.argv[1], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else None)
    else:
        main(sys.argv[1], sys.argv[2:] == ["--frame"])
