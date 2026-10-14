# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Clean CPU timing for the head-to-head: a stratified subset, one sample at a time (no other
benchmark processes), each tool's own CPU via wait4, on the files the head-to-head already cached.

Usage: python timing.py <results.tsv> <cache> <tool> <wl> <sanger.py> <sanger_wl> <out.tsv> [per_stratum]
"""
import csv, glob, os, random, sys
sys.path.insert(0, os.path.dirname(__file__)); sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "arc", "scripts"))
from run import cpu_of  # noqa: E402

def files_for(r, cache):
    if r["source"] == "SRA":
        return sorted(glob.glob(f"{cache}/sra_{r['srr']}/{r['srr']}_*.fastq.gz"))
    return sorted(glob.glob(f"{cache}/{r['srr']}*.fastq.gz"))

def main(results, cache, tool, wl, sanger, sanger_wl, out_path, per=10):
    rows = [r for r in csv.DictReader(open(results), delimiter="\t") if r["condition"] == "ok"]
    rng, chosen = random.Random(7), []
    for key in sorted({(r["arc_whitelist"], r["source"]) for r in rows}):
        group = [r for r in rows if (r["arc_whitelist"], r["source"]) == key]
        rng.shuffle(group); chosen += group[:int(per)]
    with open(out_path, "w") as out:
        out.write("srx\tarc_whitelist\tsource\tours_cpu\tsanger_cpu\n")
        for r in chosen:
            files = files_for(r, cache)
            pair = [os.path.join(os.path.dirname(files[0]), p) for p in r["ours_pair"].split(",") if p]
            _, _, ours = cpu_of([tool, "--whitelists", wl, "--r1", pair[0]] + (["--r2", pair[1]] if len(pair) > 1 else []))
            j = os.path.join(cache, f"timing_{r['srx']}.json")
            _, _, sang = cpu_of(["python3", "-I", sanger, *files, "--sample", "s", "--whitelists", sanger_wl, "--outdir",
                                 os.path.join(cache, "timing_out"), "--json", j, "--manifest", j + ".tsv", "--dry-run",
                                 "--allow-missing-whitelists"])
            out.write(f"{r['srx']}\t{r['arc_whitelist']}\t{r['source']}\t{ours:.2f}\t{sang:.2f}\n"); out.flush()

if __name__ == "__main__":
    main(*sys.argv[1:8], *sys.argv[8:9])
