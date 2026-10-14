# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Reads public runs by accession (`tchem <accession>`, streamed from ENA) and checks each call
against the same run's call from its cached head (bench/coverage). The cached heads are the first
reads of the same ENA files, so the two calls must be the same: layout, verdict and confidence.

Runs are picked without looking at their calls: the first PER_SOURCE runs of each held-out
source, by accession. Needs the network and curl; dated by the run, like any download.

Usage: python run.py <tchem> <whitelist_dir> <calls.tsv> <out.tsv>
"""

import os
import csv
import datetime
import json
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

SOURCES = ("flex", "bd", "mgi", "parse", "pipseq")
PER_SOURCE = 4
TIMEOUT_SECONDS = 1200
FIELDS = ["source", "run", "local", "remote", "agree", "seconds", "megabytes"]


def call(report):
    layout = report["layout"]["name"] if report["layout"] else ""
    return f"{report['verdict']} {layout} {report['confidence'] or ''}".strip()


def main(tchem, whitelists, calls_path, out_path):
    with open(calls_path) as f:
        rows = list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))
    picked = []
    for source in SOURCES:
        ok = sorted((r for r in rows if r["source"] == source and r["status"] == "ok"),
                    key=lambda r: r["run"])
        picked += ok[:PER_SOURCE]
    invariant(len(picked) == PER_SOURCE * len(SOURCES), "every source has enough runs")
    results = []
    for r in picked:
        local = f"{r['verdict']} {r['layout']} {r['confidence']}".strip()
        start = datetime.datetime.now()
        p = subprocess.run([tchem, "--whitelists", whitelists, "--format", "json", r["run"]],
                           capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
        seconds = (datetime.datetime.now() - start).total_seconds()
        if p.returncode != 0:
            remote = next((l for l in p.stderr.splitlines() if l.startswith("Error")), "failed")
        else:
            remote = call(json.loads(p.stdout))
        streamed = re.search(r"streamed from ENA \(([0-9.]+) MB\)", p.stderr)
        megabytes = streamed.group(1) if streamed else ""
        results.append({"source": r["source"], "run": r["run"], "local": local, "remote": remote,
                        "agree": remote == local, "seconds": f"{seconds:.1f}",
                        "megabytes": megabytes})
        print(f"{r['source']:7} {r['run']:12} {remote == local!s:5} {seconds:6.1f}s "
              f"{megabytes:>6} MB  {remote}")
    with open(out_path, "w", newline="") as f:
        version = subprocess.run([tchem, "--version"], capture_output=True, text=True).stdout
        f.write(f"# Runs read by accession from ENA vs their cached heads ({version.strip()}, "
                f"{datetime.date.today().isoformat()}).\n")
        w = csv.DictWriter(f, FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(results)
    agree = sum(r["agree"] for r in results)
    sizes = sorted(float(r["megabytes"]) for r in results if r["megabytes"])
    if sizes:
        print(f"streamed {sizes[0]:.1f} to {sizes[-1]:.1f} MB per run, median "
              f"{sizes[len(sizes) // 2]:.1f} MB")
    print(f"{agree} of {len(results)} agree")
    if agree != len(results):
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
