# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Re-runs only tussl-chemistry over the head-to-head's cached files (no download, no Sanger), for a
new build: verdict, confidence and library type per sample.

Usage: python rescore.py <h2h results.tsv> <cache> <tool> <wl> <out.tsv>
"""
import csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

def one(r, cache, tool, wl):
    d = f"{cache}/sra_{r['srr']}" if r["source"] == "SRA" else cache
    pair = [os.path.join(d, p) for p in r["ours_pair"].split(",") if p]
    out = subprocess.run([tool, "--whitelists", wl, "--r1", pair[0]] + (["--r2", pair[1]] if len(pair) > 1 else []),
                         capture_output=True, text=True, timeout=900).stdout
    v = next((l[9:].strip() for l in out.splitlines() if l.startswith("VERDICT")), "")
    lib = next((l.split("-> ")[-1].strip() for l in out.splitlines() if l.strip().startswith("library: constant")), "")
    strand = next((l.split("-> ")[-1].strip() for l in out.splitlines() if "cDNA read strand" in l), "")
    return {"srx": r["srx"], "verdict": v, "library": lib, "strand": strand}

def main(results, cache, tool, wl, out_path):
    rows = [r for r in csv.DictReader(open(results), delimiter="\t") if r["condition"] == "ok"]
    with ThreadPoolExecutor(8) as ex, open(out_path, "w") as out:
        w = csv.DictWriter(out, ["srx", "verdict", "library", "strand"], delimiter="\t")
        w.writeheader()
        for rec in ex.map(lambda r: one(r, cache, tool, wl), rows):
            w.writerow(rec); out.flush()

if __name__ == "__main__":
    main(*sys.argv[1:6])
