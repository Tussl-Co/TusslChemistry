# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Cell Ranger count --chemistry auto on a seeded sample of human 10x-vs-10x disputed runs (PLAN.md, 3).
Run from the repo root with TUSSL_CHEMISTRY_DATA set; one run at a time (the aligner index needs ~32 GB).

Usage: python cellranger_sample.py <data> <cellranger> <transcriptome> <species.tsv> <n> <workdir> <out.tsv>
"""
import csv
import glob
import os
import random
import re
import subprocess
import sys
import time

sys.path.insert(0, "bench/coverage/scripts")
sys.path.insert(0, "bench")
import sources  # noqa: E402

SEED = 20261010


def rows(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def detected(run_dir, log):
    """The chemistry Cell Ranger settled on, from its pipeline outputs or its log."""
    for p in glob.glob(f"{run_dir}/**/chemistry_defs.json", recursive=True) + \
            glob.glob(f"{run_dir}/**/chemistry_def.json", recursive=True):
        found = re.findall(r'"name"\s*:\s*"([^"]+)"', open(p).read())
        if found:
            return found[0]
    m = re.search(r"[Cc]hemistry[^\n]*?:\s*([^\n]+)", log)
    return m.group(1).strip()[:120] if m else ""


def main(data, cellranger, transcriptome, species_path, n, work, out_path):
    calls = {r["run"]: r for r in rows("bench/coverage/results/calls.tsv")}
    disputed = {r["run"] for r in rows("bench/coverage/results/adjudication.tsv") if r["outcome"] == "disputed"}
    human = {r["run"] for r in rows(species_path) if r["species"] == "Homo sapiens"}
    pool = sorted((c for c in sources.cases(data) if c.run in disputed & human
                   and calls[c.run]["confidence"] == "high"
                   and calls[c.run]["layout"].startswith("10x") and calls[c.run]["truth"].startswith("10x")
                   and all(os.path.exists(f) for f in c.files)), key=lambda c: c.run)
    sample = random.Random(SEED).sample(pool, min(int(n), len(pool)))
    print(len(pool), "eligible;", len(sample), "sampled", flush=True)
    os.makedirs(work, exist_ok=True)
    fields = ["run", "study", "label", "call", "cellranger", "status", "seconds", "message"]
    done = {r["run"] for r in rows(out_path)} if os.path.exists(out_path) else set()
    with open(out_path, "a") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        if not done:
            w.writeheader()
        for c in sample:
            if c.run in done:
                continue
            fq = os.path.join(work, "fastq", c.run)
            os.makedirs(fq, exist_ok=True)
            for i, src in enumerate(c.files[:2], 1):
                dst = os.path.join(fq, f"{c.run}_S1_L001_R{i}_001.fastq.gz")
                if not os.path.exists(dst):
                    os.symlink(src, dst)
            t = time.time()
            p = subprocess.run([cellranger, "count", f"--id={c.run}", f"--transcriptome={transcriptome}",
                                f"--fastqs={fq}", f"--sample={c.run}", "--create-bam=false", "--nosecondary",
                                "--chemistry=auto", "--localcores=16", "--localmem=48", "--disable-ui"],
                               cwd=work, capture_output=True, text=True)
            log = p.stdout + p.stderr
            err = re.findall(r"(?:error|Error|ERROR)[^\n]*", log)
            w.writerow({"run": c.run, "study": c.study, "label": calls[c.run]["truth"], "call": calls[c.run]["layout"],
                        "cellranger": detected(os.path.join(work, c.run), log),
                        "status": "ok" if p.returncode == 0 else f"exit {p.returncode}",
                        "seconds": f"{time.time() - t:.0f}", "message": (err[0] if err else "")[:200]})
            f.flush()
            subprocess.run(["rm", "-rf", os.path.join(work, c.run, "SC_RNA_COUNTER_CS")])


if __name__ == "__main__":
    main(*sys.argv[1:8])
