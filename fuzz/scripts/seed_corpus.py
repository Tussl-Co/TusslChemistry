# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Seeds the fuzz corpora from held-out read heads: up to two runs per layout tchem called HIGH.
  detect/<n>   the first 300 records of R1, a 0xFF byte, the first 300 of R2 (detect.rs splits there)
  reads/<n>_*  the first 8 records of R1 as plain FASTQ, gzip, BGZF and an unaligned BAM
Run from the repository root with TUSSL_CHEMISTRY_DATA set (heads as bench/coverage reads them).

Usage: python fuzz/scripts/seed_corpus.py <data> <corpus dir>
"""

import collections
import csv
import gzip
import os
import sys

sys.path[:0] = ["bench/coverage/scripts", "bench", "bench/chaos/scripts"]
import probes  # noqa: E402  (bgzf and bam writers)
import sources  # noqa: E402

PER_LAYOUT = 2
PAIR_RECORDS = 300
SMALL_RECORDS = 8


def head(path, records):
    with gzip.open(path, "rb") as f:
        return b"".join(f.readline() for _ in range(records * 4))


def main(data, out):
    os.makedirs(f"{out}/reads", exist_ok=True)
    os.makedirs(f"{out}/detect", exist_ok=True)
    calls = {r["run"]: r for r in csv.DictReader(
        (l for l in open("bench/coverage/results/calls.tsv") if not l.startswith("#")), delimiter="\t")}
    seen, n = collections.Counter(), 0
    for case in sources.cases(data):
        call = calls.get(case.run)
        if not call or call["confidence"] != "high" or seen[call["layout"]] >= PER_LAYOUT:
            continue
        if not all(os.path.exists(f) for f in case.files):
            continue
        seen[call["layout"]] += 1
        n += 1
        pair = head(case.files[0], PAIR_RECORDS) + b"\xff" + head(case.files[1], PAIR_RECORDS)
        open(f"{out}/detect/{n:03d}", "wb").write(pair)
        small = head(case.files[0], SMALL_RECORDS)
        open(f"{out}/reads/{n:03d}_plain", "wb").write(small)
        open(f"{out}/reads/{n:03d}_gz", "wb").write(gzip.compress(small))
        open(f"{out}/reads/{n:03d}_bgzf", "wb").write(probes.bgzf(small))
        lines = small.split(b"\n")
        records = [(lines[i][1:].split()[0], lines[i + 1]) for i in range(0, SMALL_RECORDS * 4, 4)]
        probes.bam(f"{out}/reads/{n:03d}_bam", records)
    print(f"{n} runs from {len(seen)} layouts")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
