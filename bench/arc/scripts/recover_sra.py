# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Recover benchmark samples whose ENA export was unusable (barcode read dropped, or no FASTQ) by
pulling the first 200k spots from NCBI SRA with every read kept (fastq-dump includes technical
reads by default), then scoring them exactly as run_bench.py does.

Usage: python recover_sra.py <results.tsv> <cache_dir> <tool> <whitelist_dir> <fastq-dump> <out.tsv>
"""

import csv
import glob
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(__file__))
import run_bench as rb  # noqa: E402

SPOTS = 200_000


def runs(srx):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + srx
           + "&result=read_run&fields=run_accession&format=json")
    return [r["run_accession"] for r in json.loads(rb.get(url) or b"[]")]


def recover(row, cache, tool, wl, fastq_dump):
    rec = {k: row[k] for k in ("srx", "tech_10x", "arc_whitelist")}
    try:
        srr = row.get("srr") or (runs(row["srx"]) or [""])[0]
        if not srr:
            return {**rec, "condition": "sra: no run"}
        out = os.path.join(cache, srr)
        os.makedirs(out, exist_ok=True)
        if not glob.glob(f"{out}/{srr}_*.fastq.gz"):
            subprocess.run([fastq_dump, "-X", str(SPOTS), "--split-files", "--gzip", "-O", out, srr],
                           capture_output=True, timeout=900, check=True)
        paths = sorted(glob.glob(f"{out}/{srr}_*.fastq.gz"), key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))
        lens = [rb.mode_len(p) for p in paths]
        usable = [p for p, n in zip(paths, lens) if n >= rb.MIN_READ_LEN]
        rec.update(srr=srr, n_files=len(paths), file_lens=",".join(map(str, lens)),
                   pair_used=",".join(os.path.basename(p) for p in usable[:2]))
        if not usable:
            return {**rec, "condition": "sra: no usable reads"}
        args = [tool, "--whitelists", wl, "--r1", usable[0]] + (["--r2", usable[1]] if len(usable) > 1 else [])
        text = subprocess.run(args, capture_output=True, text=True, timeout=600).stdout
        name, conf, swapped, (pa, pt) = rb.parse(text)
        our_wl, our_end = rb.VERDICT_MAP.get(name, ("", ""))
        if name.startswith("AMBIGUOUS between 10x 3' v2"):
            our_wl, our_end = "737K-august-2016", "ambiguous"
        arc_end = rb.ARC_END.get(row["tech_10x"], "out of scope")
        rec.update(condition="sra ok" if len(usable) >= 2 else "sra: single read file", verdict=name,
                   confidence=conf, swapped=swapped, polyA=pa, polyT=pt, our_whitelist=our_wl, our_end=our_end,
                   arc_end=arc_end, wl_agree=our_wl == row["arc_whitelist"], end_agree=our_end == arc_end)
        return rec
    except Exception as e:
        return {**rec, "condition": "sra error: " + repr(e)[:100]}


def main(results, cache, tool, wl, fastq_dump, out_path):
    rows = [r for r in csv.DictReader(open(results), delimiter="\t")
            if r["condition"] in ("single read file", "no ENA fastq") or r["condition"].startswith("error")]
    with ThreadPoolExecutor(4) as ex, open(out_path, "w") as out:
        w = csv.DictWriter(out, rb.FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for i, rec in enumerate(ex.map(lambda r: recover(r, cache, tool, wl, fastq_dump), rows), 1):
            w.writerow(rec)
            out.flush()
            print(f"[{i}/{len(rows)}] {rec['srx']} {rec.get('condition')} {rec.get('verdict', '')}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:7])
