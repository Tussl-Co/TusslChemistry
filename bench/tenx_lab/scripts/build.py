# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Copy of bench/head_to_head/scripts/lab_labels/build.py (diff: work/build.diff). Changes: label10x imported
from the committed directory; input = every ENA run of the remaining SRX (srx, srr, study); labels.tsv gains a
run column and loses conflict_reason (no reasons.tsv); the 465 assert is dropped; each run's label, evidence
(first sentence that stated its exact label) and missing records are saved to work/run_labels.json.
Usage: python3 -I -B build.py <W> <outdir>"""
import csv, json, os, re, sys
from collections import Counter

A, OUT = sys.argv[1], sys.argv[2]
sys.path.insert(0, os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."), "head_to_head", "scripts", "lab_labels"))
import label10x as L  # noqa: E402

ENC = re.compile(r"&#|&prime;|&[a-z]+;|â€|�|Â", re.I)
DECISIVE = [L.FIVE, L.THREE, L.VER, L.GEMX, L.MULTI, L.FLEX1, L.FLEX2, L.CONTEXT, L.GUIDE]


def sentences(text):
    out = []
    for line in text.splitlines():
        out += [s.strip() for s in re.split(r"(?<=[.;])\s+", line) if s.strip()]
    return out


def quote(s):
    s = re.sub(r"\t+", " ", s)
    if len(s) <= 300:
        return s
    starts = sorted(m.start() for rx in DECISIVE for m in rx.finditer(s)) or [0]
    wins = [s[a:a + 300] for a in dict.fromkeys(max(0, min(x - 40, len(s) - 300)) for x in starts)]
    # first window (from the earliest decisive match) that keeps the label and every CG number; else the earliest
    for q in wins:
        if L.label(q) == L.label(s) and set(L.guide_numbers(s)) <= set(L.guide_numbers(q)):
            return q
    return wins[0]


def mapped_text(s):
    m, why = L.analyse(s)
    if m == "ambiguous":
        return f"ambiguous ({why})", m
    if m:
        return m, m
    if why.startswith("skipped"):
        return why, None
    g = L.guide_numbers(s)
    if g and not L.is_candidate(s):
        return f"none (guide {','.join(g)})", None
    return ("none" if why == "no end" else f"none ({why})"), None


csv.field_size_limit(10**8)
inp = list(csv.DictReader(open(f"{A}/work/input_runs.tsv", encoding="utf-8"), delimiter="\t"))
ena_idx = json.load(open(f"{A}/raw/ena/index.json", encoding="utf-8"))
ena = {}
for c in ena_idx:
    if not c["path"]:
        continue
    for r in csv.DictReader(open(c["path"], encoding="utf-8"), delimiter="\t"):
        ena[r["run_accession"]] = (r, c["url"])
geo_idx = json.load(open(f"{A}/raw/geo/index.json", encoding="utf-8"))
GEO_KEYS = {"!Sample_title": "title", "!Sample_description": "description", "!Sample_extract_protocol_ch1": "extract_protocol",
            "!Sample_characteristics_ch1": "characteristics", "!Sample_data_processing": "data_processing"}


def geo_fields(g):
    out = []
    for line in open(f"{A}/raw/geo/{g}.txt", encoding="utf-8"):
        k, _, v = line.rstrip("\n").partition(" = ")
        if k in GEO_KEYS and v.strip():
            out.append((GEO_KEYS[k], v))
    return out


stm, labels, missing = [], [], []
info = {}
n_long, guide_lost = 0, 0
enc_sent, quote_mismatch, partial_vs_exact, partial_disagree = 0, 0, [], []
for run in inp:
    srx, srr, study = run["srx"], run["srr"], run["study"]
    exact, partial = [], []
    recs = []
    first_exact = {}
    run_missing = []
    if srr in ena:
        r, url = ena[srr]
        for f, rid in [("experiment_title", r["experiment_accession"]), ("sample_title", r["run_accession"]),
                       ("library_name", r["experiment_accession"]), ("library_construction_protocol", r["experiment_accession"]),
                       ("sample_alias", r["run_accession"]), ("study_title", r["study_accession"]),
                       ("experiment_alias", r["experiment_accession"])]:
            if r[f].strip():
                recs.append(("ena", rid, url, f, r[f]))
        gsms = sorted(set(re.findall(r"GSM\d+", " ".join([r["sample_alias"], r["experiment_alias"], r["experiment_title"]]))))
        for g in gsms:
            st = geo_idx.get(g)
            if not st or st["error"] or not os.path.exists(f"{A}/raw/geo/{g}.txt"):
                missing.append(g); run_missing.append(g); continue
            for fname, v in geo_fields(g):
                recs.append(("geo", g, st["url"], fname, v))
    else:
        missing.append(srr); run_missing.append(srr)
    for source, rid, url, fname, text in recs:
        for s in sentences(text):
            if not (L.is_candidate(s) or L.guide_numbers(s)):
                continue
            if ENC.search(s):
                enc_sent += 1
            mt, m = mapped_text(s)
            q = quote(s)
            if len(s) > 300 and L.label(q) != L.label(s):
                quote_mismatch += 1
            if len(s) > 300:
                n_long += 1
                guide_lost += not set(L.guide_numbers(s)) <= set(L.guide_numbers(q))
            stm.append((srx, study, source, rid, url, fname, q, mt))
            if fname == "study_title" or m is None:
                continue
            if m in L.EXACT:
                exact.append(m)
                first_exact.setdefault(m, (q, s, source, rid, url, fname))
            elif m in L.PARTIAL:
                partial.append(m)
    kinds, pk = set(exact), set(partial)
    if len(kinds) == 1:
        lab = kinds.pop()
        end = "5'" if "5'" in lab else ("3'" if "3'" in lab and "unclaimed" not in lab else None)
        if pk and (end is None or pk != {f"{end} unversioned"}):
            partial_vs_exact.append(srx)
    elif len(kinds) > 1:
        lab = "conflict"
    elif len(pk) == 1:
        lab = "end_only (" + pk.pop().split()[0] + ")"
    else:
        lab = "not stated"
        if len(pk) > 1:
            partial_disagree.append(srx)
    labels.append([srx, srr, study, lab, len(exact) + len(partial)])
    info[srr] = {"srx": srx, "study": study, "label": lab, "evidence": first_exact.get(lab), "missing": run_missing}

def w(name, header, rows):
    with open(f"{OUT}/{name}", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        wr.writerow(header); wr.writerows(rows)


w("statements.tsv", ["srx", "study", "source", "record_id", "url", "field", "quoted_sentence", "mapped"], stm)
w("labels.tsv", ["srx", "run", "study", "label", "n_statements"], labels)
json.dump(info, open(f"{A}/work/run_labels.json", "w", encoding="utf-8"), indent=0)
for row in stm:
    assert len(row[6]) <= 300, row
assert len(labels) == len(inp)

c = Counter(r[3] for r in labels)
studies = {k: len({r[2] for r in labels if r[3] == k}) for k in c}
print(f"{len(stm)} statement rows; {len(labels)} runs; missing records: {len(set(missing))} {sorted(set(missing))[:20] or 'none'}")
print(f"sentences with HTML entities or mojibake: {enc_sent}; quote windows whose label differs from the full sentence: {quote_mismatch}")
print(f"sentences over 300 chars: {n_long}; of these, with a CG number lost from the quote: {guide_lost}")
print(f"runs with an exact label and a partial of another end/kind: {len(partial_vs_exact)} {partial_vs_exact}")
print(f"runs 'not stated' with disagreeing partials: {len(partial_disagree)} {partial_disagree}")
print("mapped values:", dict(Counter(r[7].split(' (')[0] if r[7].startswith(('ambiguous', 'none', 'skipped')) else r[7] for r in stm)))
for k in sorted(c, key=lambda k: -c[k]):
    print(f"  {k}\t{c[k]} runs\t{studies[k]} studies")
