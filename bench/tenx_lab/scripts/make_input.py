# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Apply the study exclusions to the ENA-mapped SRX and write work/input_runs.tsv (srx, srr, study).
Writes work/excluded_pre.tsv (srx, reason) and work/unresolved_arc.tsv. Usage: python3 -I -B make_input.py <W>"""
import csv, glob, json, sys
from collections import Counter
D = sys.argv[1]
csv.field_size_limit(10**8)
P = json.load(open(f"{D}/work/prep.json"))
idx = json.load(open(f"{D}/raw/ena/index.json"))
assert not any(c["error"] for c in idx), "an ENA chunk failed"
runs = {}
for c in idx:
    for r in csv.DictReader(open(c["path"], encoding="utf-8"), delimiter="\t"):
        runs.setdefault(r["experiment_accession"], []).append(r)
used_study = dict(P["used_study"])
unresolved = []
for x in dict.fromkeys(P["arc"]):
    if x not in runs:
        unresolved.append(x); continue
    for r in runs[x]:
        used_study.setdefault(r["study_accession"], "used_study:arc/selection")
excl, inp = [], []
pool_srx = list(dict.fromkeys(P["remaining"]))
for x in pool_srx:
    if x not in runs:
        excl.append((x, "no ENA read_run record (chunk fetched without error)")); continue
    studies = sorted({r["study_accession"] for r in runs[x]})
    bad = [used_study[s] for s in studies if s in used_study]
    if bad:
        excl.append((x, bad[0])); continue
    for r in sorted(runs[x], key=lambda r: r["run_accession"]):
        inp.append((x, r["run_accession"], r["study_accession"]))
with open(f"{D}/work/input_runs.tsv", "w", encoding="utf-8") as f:
    f.write("srx\tsrr\tstudy\n"); f.writelines("\t".join(t) + "\n" for t in inp)
with open(f"{D}/work/excluded_pre.tsv", "w", encoding="utf-8") as f:
    f.writelines(f"{a}\t{b}\n" for a, b in excl)
with open(f"{D}/work/unresolved_arc.tsv", "w", encoding="utf-8") as f:
    f.writelines(x + "\n" for x in unresolved)
json.dump(used_study, open(f"{D}/work/used_study.json", "w"), indent=0)
print(len(inp), "runs in", len({t[0] for t in inp}), "SRX,", len({t[2] for t in inp}), "studies;",
      "excluded SRX:", Counter(b.split(":")[0] for a, b in excl), "; arc SRX without ENA record:", len(unresolved))
