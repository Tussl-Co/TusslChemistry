# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Run tussl-chemistry over the frozen Arc benchmark selection and record, per sample, the input
condition (ENA files, read lengths, pair used) separately from the verdict.

For each SRX: ENA's first run with FASTQs; the first ~25 MB of each file, re-packed as 200k whole
records (cached); files under 20 bases (index reads) dropped; the first two remaining files passed
as R1/R2 in ENA's order, so mislabelled archives are tested as found.

Usage: python run_bench.py <selection.tsv> <cache_dir> <tool> <whitelist_dir> <out.tsv>
Standard library only.
"""

import csv
import gzip
import json
import os
import re
import subprocess
import sys
import urllib.request
import zlib
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

PREFIX_BYTES = 25_000_000
READS = 200_000
MIN_READ_LEN = 20

# Our verdict -> (whitelist, end); the frozen mapping (see selection header).
VERDICT_MAP = {
    "10x 3' v2": ("737K-august-2016", "3'"),
    "10x 5' v1/v2": ("737K-august-2016", "5'"),
    "10x 3' v3/v3.1": ("3M-february-2018", "3'"),
    "10x 3' v4 (GEM-X)": ("3M-3pgex-may-2023", "3'"),
    "10x 5' v3 (GEM-X)": ("3M-5pgex-jan-2023", "5'"),
}
ARC_END = {"3_prime_gex": "3'", "5_prime_gex": "5'"}


def get(url, rng=None):
    headers = {"User-Agent": "curl/8.5.0"}
    if rng:
        headers["Range"] = f"bytes=0-{rng - 1}"
    for attempt in range(4):
        try:
            return urllib.request.urlopen(urllib.request.Request(url, headers=headers), timeout=60).read()
        except Exception:
            if attempt == 3:
                raise


def runs_of(srx):
    url = ("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + srx
           + "&result=read_run&fields=run_accession,fastq_ftp,fastq_bytes&format=json")
    return json.loads(get(url) or b"[]")


def damaged_record(lines):
    """Index of the first record that is not a FASTQ record, or None. A prefix of a single-member gzip
    is decoded without its CRC (that sits at the file's end), so bytes lost in transfer can decode
    into plausible text; on 10 Oct 2026, 6 of 10,732 held-out heads had such records."""
    for i in range(0, len(lines) - 3, 4):
        h, s, plus, q = lines[i:i + 4]
        if not h.startswith(b"@") or not plus.startswith(b"+") or len(s) != len(q):
            return i // 4
    return None


def fetch_head(url, out, attempts=3):
    if os.path.exists(out):
        return
    for attempt in range(attempts):
        parts, buf = [], get("https://" + url, PREFIX_BYTES)
        while buf:
            d = zlib.decompressobj(31)
            try:
                parts.append(d.decompress(buf))
            except zlib.error:
                break
            if not d.eof:
                break
            buf = d.unused_data
        lines = b"".join(parts).split(b"\n")
        n = min((len(lines) - 1) // 4, READS) * 4
        # A body that is not gzip (a redirect, an error page) decodes to no records, none damaged.
        if n == 0:
            if attempt == attempts - 1:
                raise OSError(f"no reads on {attempts} fetches")
            continue
        bad = damaged_record(lines[:n])
        if bad is None:
            break
        if attempt == attempts - 1:
            raise OSError(f"record {bad + 1} damaged on {attempts} fetches")
    with gzip.open(out + ".tmp", "wb") as f:
        f.write(b"\n".join(lines[:n]) + b"\n")
    os.rename(out + ".tmp", out)


def mode_len(path):
    lens = Counter()
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                lens[len(line.rstrip("\n"))] += 1
            if i > 40_000:
                break
    return lens.most_common(1)[0][0] if lens else 0


def parse(out):
    v = re.search(r"^VERDICT\s+(.*)$", out, re.M)
    verdict = v.group(1).strip() if v else "NO OUTPUT"
    m = re.match(r"(.+?)\s+\(confidence (\w+)\)", verdict)
    name, conf = (m.group(1), m.group(2)) if m else (verdict, "")
    strand = re.search(r"polyA runs ([\d.]+)% vs polyT runs ([\d.]+)%", out)
    return name, conf, "SWAPPED" in out, strand.groups() if strand else ("", "")


def bench(row, cache, tool, wl):
    rec = {"srx": row["srx"], "tech_10x": row["tech_10x"], "arc_whitelist": row["arc_whitelist"]}
    try:
        runs = [r for r in runs_of(row["srx"]) if r.get("fastq_ftp")]
        if not runs:
            return {**rec, "condition": "no ENA fastq"}
        run = runs[0]
        urls = run["fastq_ftp"].split(";")
        paths, lens = [], []
        for u in urls:
            p = os.path.join(cache, os.path.basename(u))
            fetch_head(u, p)
            paths.append(p)
            lens.append(mode_len(p))
        usable = [p for p, n in zip(paths, lens) if n >= MIN_READ_LEN]
        rec.update(srr=run["run_accession"], n_files=len(urls), file_lens=",".join(map(str, lens)),
                   pair_used=",".join(os.path.basename(p) for p in usable[:2]))
        if not usable:
            return {**rec, "condition": "no usable reads"}
        args = [tool, "--whitelists", wl, "--r1", usable[0]] + (["--r2", usable[1]] if len(usable) > 1 else [])
        out = subprocess.run(args, capture_output=True, text=True, timeout=600).stdout
        name, conf, swapped, (pa, pt) = parse(out)
        our_wl, our_end = VERDICT_MAP.get(name, ("", ""))
        if name.startswith("AMBIGUOUS between 10x 3' v2"):
            our_wl, our_end = "737K-august-2016", "ambiguous"
        arc_end = ARC_END.get(row["tech_10x"], "out of scope")
        rec.update(condition="ok" if len(usable) >= 2 else "single read file", verdict=name, confidence=conf,
                   swapped=swapped, polyA=pa, polyT=pt, our_whitelist=our_wl, our_end=our_end, arc_end=arc_end,
                   wl_agree=our_wl == row["arc_whitelist"], end_agree=our_end == arc_end)
        return rec
    except Exception as e:
        return {**rec, "condition": "error: " + repr(e)[:100]}


FIELDS = ["srx", "srr", "tech_10x", "arc_whitelist", "condition", "n_files", "file_lens", "pair_used",
          "verdict", "confidence", "swapped", "polyA", "polyT", "our_whitelist", "our_end", "arc_end",
          "wl_agree", "end_agree"]


def main(sel, cache, tool, wl, out_path):
    os.makedirs(cache, exist_ok=True)
    rows = list(csv.DictReader((l for l in open(sel) if not l.startswith("#")), delimiter="\t"))
    with ThreadPoolExecutor(4) as ex, open(out_path, "w") as out:
        w = csv.DictWriter(out, FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for i, rec in enumerate(ex.map(lambda r: bench(r, cache, tool, wl), rows), 1):
            w.writerow(rec)
            out.flush()
            print(f"[{i}/{len(rows)}] {rec['srx']} {rec.get('condition')} {rec.get('verdict', '')}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:6])
