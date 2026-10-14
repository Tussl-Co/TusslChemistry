# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Parity oracle: every input the benchmarks have used, run through one binary, each report kept
byte for byte. A refactor passes only when `compare` finds every report identical (stdout, exit
code and error line; timing on stderr is excluded).

Usage: python oracle.py record|compare|verdicts <tool> <whitelist_dir> <oracle_dir>
  verdicts: for a deliberate behaviour change; lists every case whose VERDICT block changed, so
  each change is reviewed before the oracle is re-recorded.
Standard library only. Inputs live under $TUSSL_CHEMISTRY_DATA (see cases()).
"""

import csv
import glob
import hashlib
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "chaos", "scripts"))
import chaos  # noqa: E402

RNA = os.environ["TUSSL_CHEMISTRY_DATA"]  # root of the benchmark inputs
EXTRA = [  # validation inputs not already among the chaos bases
    ("v3_5k", f"{RNA}/samples/5k_pbmc_v3/R1.fastq.gz", f"{RNA}/samples/5k_pbmc_v3/R2.fastq.gz"),
    ("v3_10k", f"{RNA}/samples/pbmc_10k_v3/R1.fastq.gz", f"{RNA}/samples/pbmc_10k_v3/R2.fastq.gz"),
    ("v31_nextgem", f"{RNA}/samples/10k_PBMC_3p_nextgem_Chromium_Controller/R1.fastq.gz",
     f"{RNA}/samples/10k_PBMC_3p_nextgem_Chromium_Controller/R2.fastq.gz"),
    ("v3_trim26", f"{RNA}/samples/preflight_mangled/v3_trim26_R1.fastq.gz", f"{RNA}/samples/preflight_mangled/v3_R2_head.fastq.gz"),
    ("wgs", f"{RNA}/giab-chr20-35x/HG002_chr20_R1.fq", f"{RNA}/giab-chr20-35x/HG002_chr20_R2.fq"),
    ("v2_no_r2", f"{RNA}/samples/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R1_001.fastq.gz", None),
    ("visium", f"{RNA}/samples/ena_SRR24710534/SRR24710534_1.fastq.gz", f"{RNA}/samples/ena_SRR24710534/SRR24710534_2.fastq.gz"),
]


def bench_rows(results, cache_of):
    for r in csv.DictReader(open(results), delimiter="\t"):
        pair = [p for p in (r.get("pair_used") or "").split(",") if p]
        if pair:
            yield r["srx"], *[os.path.join(cache_of(r), p) for p in pair[:2]], *([None] if len(pair) == 1 else [])


def cases():
    out = [(f"base_{n}", r1, r2) for n, r1, r2, *_ in chaos.BASES] + [(f"extra_{n}", r1, r2) for n, r1, r2 in EXTRA]
    for r1 in sorted(glob.glob(f"{RNA}/chaos_work/*.R1.*")):
        name = os.path.basename(r1).split(".R1.")[0]
        r2 = r1.replace(".R1.", ".R2.")
        out.append((f"chaos_{name}", r1, r2 if os.path.exists(r2) else None))
    out += [(f"wild_{s}", a, b) for s, a, b in bench_rows(f"{RNA}/wild/results.tsv", lambda r: f"{RNA}/wild/cache")]
    out += [(f"arc_{s}", a, b) for s, a, b in bench_rows(f"{RNA}/arc_bench/results.tsv", lambda r: f"{RNA}/arc_bench/cache")]
    out += [(f"sra_{s}", a, b) for s, a, b in bench_rows(f"{RNA}/arc_bench/results_sra.tsv",
                                                         lambda r: f"{RNA}/arc_bench/sra/{r['srr']}")]
    return out


# Extra arguments for every run, e.g. ORACLE_ARGS="--format json" with its own oracle_dir.
EXTRA_ARGS = os.environ.get("ORACLE_ARGS", "").split()


def run(tool, wl, r1, r2):
    args = [tool, "--whitelists", wl, "--r1", r1] + (["--r2", r2] if r2 else []) + EXTRA_ARGS
    p = subprocess.run(args, capture_output=True, text=True, timeout=900)
    errors = [l for l in p.stderr.splitlines() if l.startswith("Error") or "panicked" in l]
    return f"rc={p.returncode}\n{p.stdout}\n--- stderr errors\n" + "\n".join(errors) + "\n"


def main(mode, tool, wl, odir):
    os.makedirs(odir, exist_ok=True)
    todo = cases()
    with ThreadPoolExecutor(8) as ex:
        results = list(ex.map(lambda c: (c[0], run(tool, wl, c[1], c[2])), todo))
    diffs = 0
    if mode == "verdicts":
        def verdict(t):
            lines = t.splitlines()
            start = next((i for i, l in enumerate(lines) if l.startswith("VERDICT")), None)
            block = lines[start:start + 3] if start is not None else [l for l in lines if l.startswith("rc=") or "Error" in l]
            return " / ".join(l.strip() for l in block if l.strip() and not l.strip().startswith("hint"))
        for name, text in sorted(results):
            path = os.path.join(odir, name + ".out")
            old = verdict(open(path).read()) if os.path.exists(path) else "(new case)"
            new = verdict(text)
            if old != new:
                diffs += 1
                print(f"CHANGED {name}\n    was: {old[:160]}\n    now: {new[:160]}")
        print(f"verdicts: {len(results)} cases, {diffs} changed")
        sys.exit(0)
    for name, text in results:
        path = os.path.join(odir, name + ".out")
        if mode == "record":
            open(path, "w").write(text)
        elif not os.path.exists(path) or open(path).read() != text:
            diffs += 1
            print(f"DIFF {name}")
    digest = hashlib.sha256("".join(t for _, t in sorted(results)).encode()).hexdigest()[:16]
    print(f"{mode}: {len(results)} cases, digest {digest}" + ("" if mode == "record" else f", {diffs} differ"))
    sys.exit(1 if diffs else 0)


if __name__ == "__main__":
    main(*sys.argv[1:5])
