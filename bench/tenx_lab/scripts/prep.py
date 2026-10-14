# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Pool SRX (srx column only), SRX/study exclusion sets. Usage: python3 -I -B prep.py <W>
(TUSSL_CHEMISTRY_DATA set: the Arc pool is read from its arc_bench/.)"""
import csv, glob, json, os, sys
W = sys.argv[1]
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")
DATA = os.environ["TUSSL_CHEMISTRY_DATA"]
def rows(p):
    return list(csv.DictReader((l for l in open(p, encoding="utf-8") if not l.startswith("#")), delimiter="\t"))
pool = []
for p in [f"{DATA}/arc_bench/pool.tsv", f"{DATA}/arc_bench/pool_recent.tsv"]:
    with open(p, encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t"); i = hdr.index("srx")
        pool += [l.rstrip("\n").split("\t")[i] for l in f if l.strip()]
uniq = list(dict.fromkeys(pool))
used_srx = {}  # srx -> source file
used_study = {}  # study -> reason
for r in rows(f"{R}/head_to_head/results/lab_label_runs.tsv"):
    used_srx.setdefault(r["srx"], "used_srx:head_to_head/lab_label_runs"); used_study.setdefault(r["study"], "used_study:head_to_head/lab_label_runs")
for r in rows(f"{R}/head_to_head/results/studies.tsv"):
    used_srx.setdefault(r["srx"], "used_srx:head_to_head/studies")
    if r["study"].strip(): used_study.setdefault(r["study"], "used_study:head_to_head/studies")
arc = [r["srx"] for r in rows(f"{R}/arc/results/selection.tsv")]
for x in arc: used_srx.setdefault(x, "used_srx:arc/selection")
for p in sorted(glob.glob(f"{R}/*/results/catalogue.tsv")):
    rr = rows(p)
    if not rr or "split" not in rr[0]: continue
    bench = p.split("/")[-3]
    for r in rr:
        if r["split"] in ("design", "examined"):
            used_study.setdefault(r["study"], f"study_{r['split']}:{bench}/catalogue")
for r in rows(f"{R}/tenx_lab/results/examined_wild.tsv"):
    used_study.setdefault(r["study"], "study_examined:exec2/examined_wild")
remaining = [x for x in uniq if x not in used_srx]
# arc SRX need a study lookup (selection.tsv has no study column)
query = remaining + [x for x in dict.fromkeys(arc) if x not in set(remaining)]
json.dump({"pool_rows": len(pool), "pool_unique": len(uniq), "used_srx": used_srx, "used_study": used_study,
           "remaining": remaining, "arc": arc, "query": query}, open(f"{W}/work/prep.json", "w"), indent=0)
print(len(pool), "pool rows;", len(uniq), "unique;", sum(x in used_srx for x in uniq), "pool SRX already used;",
      len(remaining), "remaining;", len(query), "to query;", len(used_study), "excluded studies")
