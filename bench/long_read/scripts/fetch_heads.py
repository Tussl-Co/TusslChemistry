# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fetches the first 64 MiB of every directly downloadable run in the long-read catalogue, by HTTP
byte range, so a 150 GB BAM costs the same as a small FASTQ. Resumable: existing heads are kept.

Usage: python fetch_heads.py <catalogue.tsv> <heads_dir>
"""

import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

HEAD_BYTES = 64 << 20


def url_of(row):
    urls = [u.strip() for u in row["file_urls"].replace(";", ",").split(",") if u.strip()]
    if len(urls) != 1 or urls[0].endswith(".tar.gz") or "://" not in urls[0] and "ftp.sra.ebi.ac.uk" not in urls[0]:
        return None
    url = urls[0]
    return url if "://" in url else "https://" + url


def head_name(row, url):
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in row["run"])
    return safe + (".bam" if url.endswith(".bam") else ".fq.gz" if url.endswith(".gz") else ".fq")


def fetch(row, out_dir):
    url = url_of(row)
    if url is None:
        return row["run"], "skipped (not one direct file)"
    path = os.path.join(out_dir, head_name(row, url))
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return row["run"], "kept"
    tmp = path + ".part"
    cmd = ["curl", "-sS", "--fail", "-L", "--retry", "3", "-r", f"0-{HEAD_BYTES - 1}", "-o", tmp, url]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        return row["run"], f"failed: {result.stderr.strip()[:120]}"
    os.replace(tmp, path)
    return row["run"], f"fetched {os.path.getsize(path)} bytes"


def main(catalogue, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rows = list(csv.DictReader(open(catalogue), delimiter="\t"))
    with ThreadPoolExecutor(6) as pool:
        for run, status in pool.map(lambda r: fetch(r, out_dir), rows):
            print(f"{run}\t{status}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
