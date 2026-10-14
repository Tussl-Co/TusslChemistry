# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Draw the held-out catalogue from labels.tsv with the shared rules of bench/heldout/scripts/catalogue.py.
Writes catalogue.tsv, excluded.tsv, fetch_failures.tsv, work/shape_log.tsv, work/summary.json.
Usage: python3 -I -B draw.py <W>"""
import csv, json, os, re, sys
from collections import Counter, defaultdict
D = sys.argv[1]
sys.path.insert(0, os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."), "heldout", "scripts"))
sys.path.insert(0, os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."), "head_to_head", "scripts", "lab_labels"))
import catalogue as C  # noqa: E402
import label10x as L  # noqa: E402
csv.field_size_limit(10**8)
KEEP = {"10x 3' v2": "tenx_3p_v2", "10x 3' v3/v3.1": "tenx_3p_v3", "10x 3' v4 (GEM-X)": "tenx_3p_v4",
        "10x 5' v1/v2": "tenx_5p_v1v2", "10x 5' v3 (GEM-X)": "tenx_5p_v3", "10x Multiome GEX": "tenx_multiome"}
EVID = 120
labels = list(csv.DictReader(open(f"{D}/labels.tsv", encoding="utf-8"), delimiter="\t"))
info = json.load(open(f"{D}/work/run_labels.json", encoding="utf-8"))
ena = {}
for c in json.load(open(f"{D}/raw/ena/index.json")):
    for r in csv.DictReader(open(c["path"], encoding="utf-8"), delimiter="\t"):
        ena[r["run_accession"]] = r

def counts(rows):
    c = Counter(r["label"] for r in rows)
    return {k: [c[k], len({r["study"] for r in rows if r["label"] == k})] for k in sorted(c, key=lambda k: -c[k])}

def evidence(lab, q):
    if len(q) <= EVID:
        return q, False
    starts = sorted(m.start() for rx in [L.FIVE, L.THREE, L.VER, L.GEMX, L.MULTI, L.FLEX1, L.FLEX2, L.CONTEXT] for m in rx.finditer(q)) or [0]
    wins = [q[:EVID]] + [q[a:a + EVID] for a in dict.fromkeys(max(0, min(x - 40, len(q) - EVID)) for x in starts)]
    for w in wins:
        if L.label(w) == lab:
            return w, w != q[:EVID]
    return q[:EVID], None

fate = defaultdict(list)  # srx -> reasons for runs not drawn
failures = []
kept = [r for r in labels if r["label"] in KEEP]
for r in labels:
    if r["label"] not in KEEP:
        fate[r["srx"]].append(f"label: {r['label']}")
paired = defaultdict(list)
for r in kept:
    files = [f for f in ena[r["run"]]["fastq_ftp"].split(";") if f]
    if len(files) != 2:
        fate[r["srx"]].append(f"fastq files: {len(files)}"); continue
    paired[r["study"]].append((r, files))
shape = lambda a, b: min(a, b) >= 26 and max(a, b) >= 50
rows, log, fixed, unfixable = [], [], 0, 0
for study in sorted(paired):
    ordered = sorted(paired[study], key=lambda rf: rf[0]["run"])
    n_kept = 0
    for k, (r, files) in enumerate(ordered):
        if n_kept == C.RUNS_PER_STUDY:
            fate[r["srx"]].append(f"study cap: {C.RUNS_PER_STUDY} runs already drawn from {study}"); continue
        if k >= C.SHAPE_CHECKS_PER_STUDY:
            fate[r["srx"]].append(f"not checked: shape-check limit ({C.SHAPE_CHECKS_PER_STUDY} runs per study)"); continue
        names = sorted(files)
        lengths = [C.read_length("https://" + f) for f in names]
        ok = shape(*lengths)
        log.append((r["run"], study, *lengths, "pass" if ok else "fail"))
        for f, n in zip(names, lengths):
            if n == 0:
                failures.append(("read_length", r["run"], "https://" + f, "0: unreachable on 2 attempts, or not gzip"))
        if not ok:
            fate[r["srx"]].append(f"read shape: modal lengths {lengths[0]},{lengths[1]}"); continue
        lab = r["label"]
        q = info[r["run"]]["evidence"][0]
        ev, moved = evidence(lab, q)
        fixed += bool(moved); unfixable += moved is None
        urls = ["https://" + f for f in names]
        rows.append({"kit": KEEP[lab], "study": study, "run": r["run"], "first_public": ena[r["run"]]["first_public"],
                     "r1_url": urls[0], "r2_url": urls[1], "urls": ";".join(urls), "evidence": ev, "split": "held_out"})
        n_kept += 1
rows.sort(key=lambda r: (r["kit"], r["first_public"], r["study"], r["run"]))
C.drawn_sound(rows, False)
C.write(rows, f"{D}/catalogue.tsv")
drawn_srx = {ena[r["run"]]["experiment_accession"] for r in rows}
excl = [tuple(l.rstrip("\n").split("\t")) for l in open(f"{D}/work/excluded_pre.tsv", encoding="utf-8")]
for x in sorted(fate):
    if x not in drawn_srx:
        excl.append((x, "; ".join(dict.fromkeys(fate[x]))))
with open(f"{D}/excluded.tsv", "w", encoding="utf-8") as f:
    f.write("srx\treason\n"); f.writelines(f"{a}\t{b}\n" for a, b in excl)
for c in json.load(open(f"{D}/raw/ena/index.json")):
    if c["error"]:
        failures.append(("ena", ",".join(c["srxs"]), c["url"], c["error"]))
for g, s in sorted(json.load(open(f"{D}/raw/geo/index.json")).items()):
    if s["error"]:
        failures.append(("geo", g, s["url"], s["error"]))
with open(f"{D}/fetch_failures.tsv", "w", encoding="utf-8") as f:
    f.write("source\taccession\turl\terror\n"); f.writelines("\t".join(t) + "\n" for t in failures)
with open(f"{D}/work/shape_log.tsv", "w", encoding="utf-8") as f:
    f.write("run\tstudy\tlen_file1\tlen_file2\tshape\n"); f.writelines("\t".join(map(str, t)) + "\n" for t in log)
rowlab = [{"label": r["kit"], "study": r["study"]} for r in rows]
summ = {"before": counts(labels), "kept_labels": counts(kept), "paired": sum(len(v) for v in paired.values()),
        "after": counts(rowlab), "shape_checked": len(log), "shape_fail": sum(t[-1] == "fail" for t in log),
        "evidence_window_moved": fixed, "evidence_label_lost": unfixable, "failures": Counter(t[0] for t in failures),
        "excluded_reasons": Counter(re.sub(r"\d+", "N", b.split(";")[0].split(":")[0] if b.startswith(("used_", "study_")) else b.split(";")[0]) for a, b in excl)}
json.dump(summ, open(f"{D}/work/summary.json", "w"), indent=1)
print(json.dumps(summ, indent=1))
