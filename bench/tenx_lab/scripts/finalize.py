# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Add the pool SRX excluded as already used to excluded.tsv (once), check the catalogue, print the counts
and a seeded sample. Usage: python3 -I -B finalize.py <W>"""
import csv, json, os, random, sys
from collections import Counter
D = sys.argv[1]
sys.path.insert(0, os.path.join(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."), "head_to_head", "scripts", "lab_labels"))
import label10x as L  # noqa: E402
P = json.load(open(f"{D}/work/prep.json"))
ex = list(csv.reader(open(f"{D}/excluded.tsv", encoding="utf-8"), delimiter="\t"))
have = {r[0] for r in ex[1:]}
pool = list(dict.fromkeys(P["remaining"] + [x for x in P["used_srx"]]))
add = [(x, P["used_srx"][x]) for x in P["used_srx"] if x not in have and x not in set(P["remaining"])]
add = [(x, w) for x, w in add if x in set(json.load(open(f"{D}/work/pool_srx.json")))]
with open(f"{D}/excluded.tsv", "a", encoding="utf-8") as f:
    f.writelines(f"{x}\t{w}\n" for x, w in add)
ex = list(csv.DictReader(open(f"{D}/excluded.tsv", encoding="utf-8"), delimiter="\t"))
assert len({r["srx"] for r in ex}) == len(ex), "an SRX is excluded once"
cat = list(csv.DictReader(open(f"{D}/catalogue.tsv", encoding="utf-8"), delimiter="\t"))
KEY = {"tenx_3p_v2": "10x 3' v2", "tenx_3p_v3": "10x 3' v3/v3.1", "tenx_3p_v4": "10x 3' v4 (GEM-X)",
       "tenx_5p_v1v2": "10x 5' v1/v2", "tenx_5p_v3": "10x 5' v3 (GEM-X)", "tenx_multiome": "10x Multiome GEX"}
lab = {r["run"]: r for r in csv.DictReader(open(f"{D}/labels.tsv", encoding="utf-8"), delimiter="\t")}
for r in cat:
    assert r["split"] == "held_out" and len(r["evidence"]) <= 120 and L.label(r["evidence"]) == KEY[r["kit"]], r
    assert lab[r["run"]]["label"] == KEY[r["kit"]] and lab[r["run"]]["study"] == r["study"]
drawn_srx = {lab[r["run"]]["srx"] for r in cat}
pool_srx = set(json.load(open(f"{D}/work/pool_srx.json")))
assert drawn_srx | {r["srx"] for r in ex} == pool_srx and not drawn_srx & {r["srx"] for r in ex}, "every pool SRX accounted for"
print("added", len(add), "used-SRX rows; excluded SRX", len(ex), "; drawn SRX", len(drawn_srx), "; pool", len(pool_srx))
print("studies drawn:", len({r["study"] for r in cat}), "; studies in >1 kit:",
      sum(1 for s, n in Counter(s for s in {(r["study"], r["kit"]) for r in cat}).items() if False) or
      sum(1 for n in Counter(s for s, k in {(r["study"], r["kit"]) for r in cat}).values() if n > 1))
for r in random.Random(20261009).sample(cat, 5):
    print("\t".join(r[k] for k in ["kit", "study", "run", "split", "evidence"]), "\t", r["r1_url"].split("/")[-1], r["r2_url"].split("/")[-1])
