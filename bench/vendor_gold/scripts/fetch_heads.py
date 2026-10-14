# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fetches the head of R1 and R2 for every dataset in the vendor gold catalogue, by HTTP byte range,
so a 100 GB FASTQ tarball costs a few tens of megabytes. Inside a tar, member headers are walked one
range request at a time (GNU long names included) to find a lane's R1 and R2; each head keeps whole
records only, re-gzipped. Resumable: existing heads are kept.

Usage: python fetch_heads.py <catalogue.tsv> <heads_dir>
"""

import csv
import gzip
import os
import re
import sys
import urllib.request
import zlib
from concurrent.futures import ThreadPoolExecutor

MEMBER_BYTES = 32 << 20  # compressed bytes read from the start of each FASTQ
READS_MAX = 200_000
TAR_MEMBERS_MAX = 4_000  # a tarball with more members than this is not walked
AGENT = {"User-Agent": "curl/8.5.0"}


def ranged(url, start, length, attempts=3):
    """`ranged_once`, retried: a large file's first request can time out while the CDN fetches it."""
    for attempt in range(attempts):
        try:
            return ranged_once(url, start, length)
        except (TimeoutError, OSError):
            if attempt == attempts - 1:
                raise
    return b""


def ranged_once(url, start, length):
    """Exactly the bytes asked for. A server may answer a range request with the whole file
    (HTTP 200): then only the wanted bytes are read, never the rest of a 100 GB body."""
    req = urllib.request.Request(url, headers={**AGENT, "Range": f"bytes={start}-{start + length - 1}"})
    with urllib.request.urlopen(req, timeout=120) as r:
        if r.status == 206:
            return r.read(length)
        # Whole file sent: skip to `start` within the stream, then take `length`.
        skipped = 0
        while skipped < start:
            chunk = r.read(min(1 << 20, start - skipped))
            if not chunk:
                return b""
            skipped += len(chunk)
        return r.read(length)


def tar_size(field):
    """A tar size field: octal text, or for members of 8 GiB and more (GNU), a base-256 number
    flagged by the high bit of its first byte."""
    if field[0] & 0x80:
        return int.from_bytes(bytes([field[0] & 0x7F]) + field[1:], "big")
    return int(field.split(b"\0")[0].strip() or b"0", 8)


def tar_members(url):
    """(name, data offset, size) of every member, from the headers alone."""
    offset, long_name = 0, None
    for _ in range(TAR_MEMBERS_MAX):
        header = ranged(url, offset, 512)
        if len(header) < 512 or not header.strip(b"\0"):
            return
        name = header[:100].split(b"\0")[0].decode(errors="replace")
        prefix = header[345:500].split(b"\0")[0].decode(errors="replace")
        size = tar_size(header[124:136])
        kind = header[156:157]
        data = offset + 512
        offset = data + (size + 511) // 512 * 512
        if kind == b"L":  # GNU long name: the next member's name is this member's data
            long_name = ranged(url, data, size).split(b"\0")[0].decode(errors="replace")
            continue
        if kind in (b"x", b"g"):  # pax headers carry no FASTQ
            continue
        full = long_name or (prefix + "/" + name if prefix else name)
        long_name = None
        yield full, data, size


def records(raw):
    """Whole FASTQ records from the start of a (multi-member) gzip stream cut short."""
    parts, buffer = [], raw
    while buffer:
        d = zlib.decompressobj(31)
        try:
            parts.append(d.decompress(buffer))
        except zlib.error:
            break
        if not d.eof:
            break
        buffer = d.unused_data
    lines = b"".join(parts).split(b"\n")
    n = min((len(lines) - 1) // 4, READS_MAX) * 4
    return b"\n".join(lines[:n]) + b"\n", n // 4


# Libraries that share a tarball with gene expression but are not it.
OTHER_LIBRARIES = re.compile(
    r"atac|vdj|tcr|bcr|_t_|_b_|antibody|_ab_|crispr|feature|_fb_|cmo|hto|adt|multiplex|i7|index",
    re.I,
)


def pick(names):
    """The first lane's gene-expression R1 and R2 among a tarball's FASTQ members: files named
    for gene expression first, then any not named for another library type."""
    fastqs = sorted(n for n in names if re.search(r"\.f(ast)?q(\.gz)?$", n))
    gex = [n for n in fastqs if re.search(r"gex|gene_?exp|_ge_", n, re.I)]
    fastqs = gex or [n for n in fastqs if not OTHER_LIBRARIES.search(n.rsplit("/", 1)[-1].lower())]
    for lane in ("_L001_", "_L002_", ""):
        r1 = [n for n in fastqs if lane in n and re.search(r"_R1[_.]", n)]
        r2 = [n for n in fastqs if lane in n and re.search(r"_R2[_.]", n)]
        if r1 and r2:
            return r1[0], r2[0]
    return None, None


def fetch(row, out_dir):
    name = re.sub(r"[^A-Za-z0-9._-]", "_", row["dataset"])
    r1_out, r2_out = (os.path.join(out_dir, f"{name}_{r}.fq.gz") for r in ("R1", "R2"))
    if os.path.exists(r1_out) and os.path.exists(r2_out):
        return row["dataset"], "kept"
    try:
        if row["container"] == "tar":
            members = {n: (start, size) for n, start, size in tar_members(row["files_url"])}
            r1, r2 = pick(members)
            if r1 is None:
                return row["dataset"], f"no R1/R2 pair among {len(members)} members"
            sources = [(row["files_url"], *members[r1]), (row["files_url"], *members[r2])]
        else:
            if not row.get("r2_url"):
                return row["dataset"], "no R2 URL"
            sources = [(row["files_url"], 0, MEMBER_BYTES), (row["r2_url"], 0, MEMBER_BYTES)]
        counts = []
        for (url, start, size), out in zip(sources, (r1_out, r2_out)):
            text, n = records(ranged(url, start, min(size, MEMBER_BYTES)))
            with gzip.open(out + ".part", "wb", compresslevel=1) as f:
                f.write(text)
            os.replace(out + ".part", out)
            counts.append(n)
        return row["dataset"], f"fetched {counts[0]} / {counts[1]} reads"
    except Exception as error:  # a dataset that cannot be fetched is reported, not fatal
        return row["dataset"], f"failed: {str(error)[:120]}"


def main(catalogue, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    rows = list(csv.DictReader(open(catalogue), delimiter="\t"))
    with ThreadPoolExecutor(6) as pool:
        for dataset, status in pool.map(lambda r: fetch(r, out_dir), rows):
            print(f"{dataset}\t{status}", flush=True)


if __name__ == "__main__":
    main(*sys.argv[1:3])
