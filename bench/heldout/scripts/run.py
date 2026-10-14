# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Runs a frozen held-out catalogue once: fetches the first 200k records of each file by byte
range (the arc benchmark's fetch_head), runs tchem on each pair, and records the call beside the
layout the run's own metadata names (EXPECTED, by the catalogue's kit column). Design rows (a
`split` column of "design") are skipped: they settled the layout, so they prove nothing. Rows
whose files are not both on ENA are recorded as not run.

Usage: python run.py <catalogue.tsv> <heads_dir> <tchem> <whitelist_dir> <out.tsv>
"""

import collections
import csv
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import zlib

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "arc", "scripts"))
from run_bench import fetch_head  # noqa: E402  (shared range-safe head fetcher)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # -I leaves the script dir out
from catalogue import invariant  # noqa: E402  (always-on assertion)

FIELDS = ["kit", "study", "gse", "run", "status", "verdict", "confidence", "r1_length", "detail"]
# The layouts each catalogue kit may rightly be called. A kit whose version metadata does not fix
# its read layout accepts its family: Parse's WT v3 kits shipped in both layouts, and runs labelled
# v2 carry the v3 layout. Version agreement is reported on its own line.
PARSE = ("Parse Evercode WT", "Parse Evercode WT v3")
SEEKGENE_DD = ("SeekGene SeekOne DD 3'", "SeekGene SeekOne DD 5'")
SINGLERON = ("Singleron GEXSCOPE V1", "Singleron GEXSCOPE V2", "Singleron GEXSCOPE V3")
EXPECTED = {
    "10x_flex_v1": ("10x Flex v1",),
    # Parked: no layout yet, and its rows cannot run (R1 not on ENA).
    "10x_gemx_flex_v2": (),
    "parse_wt": PARSE,
    "parse_wt_v1": PARSE,
    "parse_wt_v2": PARSE,
    "parse_wt_v3": PARSE,
    "mgi_c4": ("MGI DNBelab C4 / C-series HT",),
    "mgi_c_ht": ("MGI DNBelab C4 / C-series HT",),
    # A protocol naming the C-series kit without C4 or HT (9 Oct 2026): one layout either way.
    "mgi_c_series": ("MGI DNBelab C4 / C-series HT",),
    # Archive metadata almost never states BD's bead version: scored as the BD family.
    "bd_rhapsody": ("BD Rhapsody (original V1 beads)", "BD Rhapsody (Enhanced Beads)",
                    "BD Rhapsody (Enhanced V2/V3 beads)"),
    "pipseq_v3": ("Illumina Single Cell 3' / PIPseq v3-v4",),
    "pipseq_v4": ("Illumina Single Cell 3' / PIPseq v3-v4",),
    # v5 (PIPseq V, Illumina Single Cell 3' RNA Prep) has its own layout from 9 Oct 2026.
    "pipseq_v5": ("Illumina Single Cell 3' RNA Prep / PIPseq V",),
    # A kit version stated beside GEXSCOPE (Library Kit V2) is taken to name the layout of that
    # version, as Singleron's own chemistry names do (an assumption, 9 Oct 2026: a disagreement
    # shows as a wrong call and goes to adjudication). Without one, the study is the family.
    "singleron": SINGLERON,
    # inDrops v1/v2 (the Klein lab gel-bead layout, TruDrop included); v3 splits the barcode
    # across R1 and an index read and is not claimed: no layout is right for it.
    "indrops": ("inDrops",),
    "indrops_v3": (),
    # SeekGene (9 Oct 2026, before any run): a DD kit with no end stated, and DD-Q (scFAST-seq,
    # FFPE), whose Read 1 is the same 17-base barcode and 12-base UMI (SeekSoul Tools coded it
    # B17U12), are right at the family level only. MM is not claimed: no layout is right for it.
    "seekgene_dd3": ("SeekGene SeekOne DD 3'",),
    "seekgene_dd5": ("SeekGene SeekOne DD 5'",),
    "seekgene_dd": SEEKGENE_DD,
    "seekgene_ddq": SEEKGENE_DD,
    "seekgene_mm": (),
    "singleron_v1": ("Singleron GEXSCOPE V1",),
    # 10x runs labelled from the lab's own record (bench/tenx_lab, label10x.py's labels).
    "tenx_3p_v2": ("10x 3' v2",),
    "tenx_3p_v3": ("10x 3' v3/v3.1",),
    "tenx_3p_v4": ("10x 3' v4 (GEM-X)",),
    "tenx_5p_v1v2": ("10x 5' v1/v2",),
    "tenx_5p_v3": ("10x 5' v3 (GEM-X)",),
    "tenx_multiome": ("10x Multiome GEX",),
    "singleron_v2": ("Singleron GEXSCOPE V2",),
    "singleron_v3": ("Singleron GEXSCOPE V3",),
}


def run_one(row, heads, tool, whitelists):
    out = {k: row.get(k, "") for k in ("kit", "study", "gse", "run")}
    urls = [row["r1_url"], row["r2_url"]]
    if not all(u.startswith("https://") for u in urls):
        return {**out, "status": "not run (R1 not on ENA)"}
    paths = [os.path.join(heads, f"{row['run']}_R{i}.fq.gz") for i in (1, 2)]
    try:
        for url, path in zip(urls, paths):
            fetch_head(url.removeprefix("https://"), path)
    except (OSError, EOFError, zlib.error) as error:  # recorded, never retried silently
        return {**out, "status": f"fetch failed: {error}"}
    run = subprocess.run([tool, "--whitelists", whitelists, "--r1", paths[0], "--r2", paths[1]],
                         capture_output=True, text=True, timeout=600)
    verdict = re.search(r"^VERDICT\s+(.*?)(?:\s+\(confidence (\w+)\))?$", run.stdout, re.M)
    length = re.search(r"^R1\s+\d+ reads, length (\d+)", run.stdout, re.M)
    candidates = re.findall(r"^  (\S.*?)\s+on R[12]\s+support\s+([\d.]+)%", run.stdout, re.M)
    return {**out, "status": "ok" if verdict else f"exit {run.returncode}",
            "verdict": verdict.group(1) if verdict else run.stderr.strip()[:200],
            "confidence": (verdict.group(2) or "") if verdict else "",
            "r1_length": length.group(1) if length else "",
            "detail": "; ".join(f"{n} {s}%" for n, s in candidates[:3])}


def main(catalogue, heads, tool, whitelists, out_path):
    os.makedirs(heads, exist_ok=True)
    rows = [r for r in csv.DictReader(open(catalogue), delimiter="\t")
            if r.get("split", "held_out") != "design"]
    unknown = {r["kit"] for r in rows} - EXPECTED.keys()
    if unknown:
        sys.exit(f"no expected layout for kit(s) {sorted(unknown)}: add them to EXPECTED")
    with ThreadPoolExecutor(6) as pool:
        results = list(pool.map(lambda r: run_one(r, heads, tool, whitelists), rows))
    invariant(len(results) == len(rows), "one result per catalogue row")
    invariant([r["run"] for r in results] == [r["run"] for r in rows], "results keep row order")
    invariant(all(r.get("status") for r in results), "every row has a status")
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(results)
    for kit in sorted({r["kit"] for r in results}):
        want = EXPECTED[kit]
        calls = [r for r in results if r["kit"] == kit and r["status"] == "ok"]
        right = sum(r["verdict"] in want for r in calls)
        high = sum(r["verdict"] in want and r["confidence"] == "HIGH" for r in calls)
        wrong = sum(r["verdict"] not in want and r["confidence"] in ("HIGH", "MEDIUM")
                    for r in calls)
        print(f"{kit}: {sum(r['kit'] == kit for r in results)} rows, {len(calls)} run: {right} "
              f"called {' or '.join(want)} ({high} HIGH), {wrong} confidently called something "
              f"else")
        if len(want) > 1:
            print(f"  by layout: {dict(collections.Counter(r['verdict'] for r in calls))}")


if __name__ == "__main__":
    main(*sys.argv[1:])
