# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fetches the head of each held-out run's two FASTQs in a catalogue, by byte range, with the
benches' shared fetcher (bench/arc/scripts/run_bench.py fetch_head), into the heads directory
bench/coverage reads (<data>/<bench>_heldout/heads/<run>_R1.fq.gz and _R2). A file already there is
kept. Every run that fails is written to the failures table, never dropped.

Usage: python fetch_heads.py <catalogue.tsv> <heads_dir> <failures.tsv>
"""

import csv
import os
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "arc",
                                "scripts"))
from run_bench import fetch_head  # noqa: E402  (shared range-safe head fetcher)

WORKERS = 6  # ENA's FTP over HTTPS serves a handful of range requests at once without refusing


def held_out(catalogue):
    with open(catalogue, newline="") as f:
        rows = list(csv.DictReader((line for line in f if not line.startswith("#")),
                                   delimiter="\t"))
    invariant(rows, f"{catalogue}: no rows")
    return [r for r in rows if r.get("split", "held_out") == "held_out"]


def fetch_run(row, heads):
    """None when both heads are in place; otherwise (run, url, why)."""
    for mate, url in ((1, row["r1_url"]), (2, row["r2_url"])):
        if not url.startswith("https://"):
            return row["run"], url, "not an https URL"
        out = os.path.join(heads, f"{row['run']}_R{mate}.fq.gz")
        try:
            fetch_head(url.removeprefix("https://"), out)
        except OSError as e:  # unreachable or refused: recorded, never dropped
            return row["run"], url, repr(e)[:200]
    return None


def main(catalogue, heads, failures_path):
    os.makedirs(heads, exist_ok=True)
    rows = held_out(catalogue)
    with ThreadPoolExecutor(WORKERS) as pool:
        failures = [f for f in pool.map(lambda r: fetch_run(r, heads), rows) if f]
    invariant(len(failures) <= len(rows), "at most one failure per run")
    with open(failures_path, "w") as f:
        f.write("run\turl\terror\n")
        for failure in failures:
            f.write("\t".join(failure) + "\n")
    print(f"{len(rows)} held-out runs, {len(failures)} failed")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
