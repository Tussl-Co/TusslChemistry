# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Head-to-head: tussl-chemistry and Sanger cellgeni's infer_10x_chemistry.py on the same files.

Per SRX: ENA's first run; the first 200k records of each file (re-packed). If ENA lacks a usable
barcode read (fewer than two files of 20+ bases, or no FASTQ at all), the first 200k spots come
from SRA instead, every read kept. Both tools get the same files: ours the first two files of
20+ bases (as run_bench.py), Sanger's every file (it expects the full set, index reads included).
Per tool: verdict, exit status, CPU seconds.

Usage: python run.py <selection.tsv> <cache> <tool> <wl> <sanger.py> <sanger_wl> <fastq-dump> <out.tsv>
"""

import csv
import glob
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "arc", "scripts"))
import run_bench as rb  # noqa: E402

SPOTS = 200_000
FIELDS = ["srx", "tech_10x", "arc_whitelist", "source", "srr", "n_files", "file_lens", "ours_pair",
          "ours_verdict", "ours_conf", "ours_cpu", "sanger_ok", "sanger_id", "sanger_conf", "sanger_error",
          "sanger_cpu", "cdna", "condition"]


def cpu_of(args, timeout=1800):
    """Runs a command; returns (stdout, rc, CPU seconds of that one process tree).

    os.wait4 reports the child's own resource usage. RUSAGE_CHILDREN deltas would also count
    other tools finishing in parallel threads, which inflated the first 176 rows' timings.
    """
    with subprocess.Popen(args, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True) as p:
        out = p.stdout.read()
        _, status, usage = os.wait4(p.pid, 0)
        p.returncode = os.waitstatus_to_exitcode(status)
    return out, p.returncode, usage.ru_utime + usage.ru_stime


def ena_files(srx, cache):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + srx
           + "&result=read_run&fields=run_accession,fastq_ftp&format=json")
    runs = json.loads(rb.get(url) or b"[]")
    if not runs:
        return "", []
    run = runs[0]
    paths = []
    for u in [u for u in (run.get("fastq_ftp") or "").split(";") if u]:
        p = os.path.join(cache, os.path.basename(u))
        rb.fetch_head(u, p)
        paths.append(p)
    return run["run_accession"], paths


def sra_files(srr, cache, fastq_dump):
    out = os.path.join(cache, "sra_" + srr)
    os.makedirs(out, exist_ok=True)
    if not glob.glob(f"{out}/{srr}_*.fastq.gz"):
        subprocess.run([fastq_dump, "-X", str(SPOTS), "--split-files", "--gzip", "-O", out, srr],
                       capture_output=True, timeout=1800, check=True)
    return sorted(glob.glob(f"{out}/{srr}_*.fastq.gz"), key=lambda p: int(p.rsplit("_", 1)[1].split(".")[0]))


def one(row, a):
    rec = {k: row[k] for k in ("srx", "tech_10x", "arc_whitelist")}
    try:
        srr, paths = ena_files(row["srx"], a["cache"])
        lens = [rb.mode_len(p) for p in paths]
        source = "ENA"
        if sum(n >= rb.MIN_READ_LEN for n in lens) < 2 and srr:
            paths, source = sra_files(srr, a["cache"], a["fastq_dump"]), "SRA"
            lens = [rb.mode_len(p) for p in paths]
        if not srr:
            return {**rec, "condition": "no run"}
        usable = [p for p, n in zip(paths, lens) if n >= rb.MIN_READ_LEN]
        rec.update(source=source, srr=srr, n_files=len(paths), file_lens=",".join(map(str, lens)),
                   ours_pair=",".join(os.path.basename(p) for p in usable[:2]))
        if not usable:
            return {**rec, "condition": "no usable reads"}

        args = [a["tool"], "--whitelists", a["wl"], "--r1", usable[0]] + (["--r2", usable[1]] if len(usable) > 1 else [])
        out, _, cpu = cpu_of(args)
        name, conf, _, _ = rb.parse(out)
        rec.update(ours_verdict=name, ours_conf=conf, ours_cpu=f"{cpu:.2f}", cdna=usable[-1] if len(usable) > 1 else "")

        jpath = os.path.join(a["cache"], f"sanger_{row['srx']}.json")
        sargs = ["python3", "-I", a["sanger"], *paths, "--sample", "s", "--whitelists", a["sanger_wl"],
                 "--outdir", os.path.join(a["cache"], "sanger_out"), "--json", jpath,
                 "--manifest", jpath + ".tsv", "--dry-run", "--allow-missing-whitelists"]
        _, rc, scpu = cpu_of(sargs)
        rep = json.load(open(jpath)) if os.path.exists(jpath) else {"ok": False, "error": f"no json, rc {rc}"}
        rec.update(sanger_ok=rep.get("ok"), sanger_id=rep.get("chemistry_id") or "", sanger_conf=rep.get("confidence") or "",
                   sanger_error=(rep.get("error") or "")[:120], sanger_cpu=f"{scpu:.2f}", condition="ok")
        return rec
    except Exception as e:
        return {**rec, "condition": "error: " + repr(e)[:120]}


def main(sel, cache, tool, wl, sanger, sanger_wl, fastq_dump, out_path):
    os.makedirs(cache, exist_ok=True)
    a = dict(cache=cache, tool=tool, wl=wl, sanger=sanger, sanger_wl=sanger_wl, fastq_dump=fastq_dump)
    rows = list(csv.DictReader((l for l in open(sel) if not l.startswith("#")), delimiter="\t"))
    done = set()
    if os.path.exists(out_path):  # resumable: keep finished rows
        done = {r["srx"] for r in csv.DictReader(open(out_path), delimiter="\t") if r.get("condition") == "ok"}
    todo = [r for r in rows if r["srx"] not in done]
    mode = "a" if done else "w"
    with ThreadPoolExecutor(8) as ex, open(out_path, mode) as out:
        w = csv.DictWriter(out, FIELDS, delimiter="\t", extrasaction="ignore")
        if not done:
            w.writeheader()
        t0 = time.time()
        for i, rec in enumerate(ex.map(lambda r: one(r, a), todo), 1):
            w.writerow(rec)
            out.flush()
            print(f"[{i}/{len(todo)} {time.time() - t0:.0f}s] {rec['srx']} {rec.get('source', '')} "
                  f"ours={rec.get('ours_verdict', '')[:28]} sanger={rec.get('sanger_id') or rec.get('sanger_error', '')[:28]} "
                  f"{rec['condition'][:40]}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:9])
