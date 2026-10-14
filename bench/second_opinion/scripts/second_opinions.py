# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Second opinions on disputed HIGH calls (PLAN.md): Sanger's infer_10x_chemistry.py and a strand
referee built from genes disjoint from tchem's. Run from the repo root with TUSSL_CHEMISTRY_DATA set.

Usage: python second_opinions.py <data> <sanger.py> <sanger_wl> <human.fa> <mouse.fa> <exclude.fa> <species.tsv> <workdir> <out.tsv>
"""
import csv
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, "bench/coverage/scripts")
sys.path.insert(0, "bench/head_to_head/scripts")
sys.path.insert(0, "bench")
import sources  # noqa: E402
import adjudicate  # noqa: E402

WORKERS = int(os.environ.get("WORKERS", "14"))
SANGER = {"10x_3p_v2_or_5p_v1_v2": {"10x 3' v2", "10x 5' v1/v2"}, "10x_3p_v3_family": {"10x 3' v3/v3.1"},
          "10x_3p_v4_gemx": {"10x 3' v4 (GEM-X)"}, "10x_5p_v3_gemx": {"10x 5' v3 (GEM-X)"},
          "10x_multiome_arc_gex": {"10x Multiome GEX"}}


def rows(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def end(layout):
    return "5'" if "5'" in layout else ("3'" if "3'" in layout else "")


def main(data, sanger, sanger_wl, human, mouse, exclude, species_path, work, out_path):
    calls = {r["run"]: r for r in rows("bench/coverage/results/calls.tsv")}
    disputed = {r["run"] for r in rows("bench/coverage/results/adjudication.tsv") if r["outcome"] == "disputed"}
    cases = [c for c in sources.cases(data) if c.run in disputed and calls[c.run]["confidence"] == "high"]
    referees = {"Homo sapiens": adjudicate.kmer_sets(human, exclude), "Mus musculus": adjudicate.kmer_sets(mouse, exclude)}
    species = {r["run"]: r["species"] for r in rows(species_path)}
    os.makedirs(work, exist_ok=True)

    def one(c):
        call = calls[c.run]
        rec = {"source": c.source, "run": c.run, "study": c.study, "label": call["truth"], "call": call["layout"],
               "heads": all(os.path.exists(f) for f in c.files)}
        if not rec["heads"]:
            return rec
        j = os.path.join(work, f"sanger_{c.run}.json")
        subprocess.run(["python3", "-I", sanger, *c.files, "--sample", "s", "--whitelists", sanger_wl,
                        "--outdir", os.path.join(work, "out"), "--json", j, "--manifest", j + ".tsv",
                        "--dry-run", "--allow-missing-whitelists"], capture_output=True, timeout=1800)
        rep = json.load(open(j)) if os.path.exists(j) else {"ok": False, "error": "no json"}
        rec.update(sanger_id=rep.get("chemistry_id") or "", sanger_error=(rep.get("error") or "")[:100])
        rec["species"] = species.get(c.run, "")
        if rec["species"] in referees:
            ss, aa, e = adjudicate.strand(c.files[-1], *referees[rec["species"]])
            rec.update(strand_sense=f"{ss:.5f}", strand_anti=f"{aa:.5f}", strand_end=e)
        return rec

    with ThreadPoolExecutor(WORKERS) as pool:
        out = list(pool.map(one, cases))
    fields = ["source", "run", "study", "species", "label", "call", "heads", "sanger_id", "sanger_error",
              "strand_sense", "strand_anti", "strand_end"]
    with open(out_path, "w") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(out)
    print(len(out), "cases;", sum(r["heads"] for r in out), "with heads")


if __name__ == "__main__":
    main(*sys.argv[1:10])
