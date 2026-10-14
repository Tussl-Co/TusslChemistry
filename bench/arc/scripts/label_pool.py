# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Label a random pool of Arc scBaseCount human samples with the 10x whitelist Arc's own read-based
pipeline (scRecounter's STARsolo grid) selected, recovered from the cell barcodes in each h5ad.

Only obs/_index is read, via HTTP range requests; the count matrix is never downloaded.

Usage: python label_pool.py <sample_metadata.parquet> <whitelist_dir> <n> <seed> <out.tsv> [min_srx [exclude.tsv]]
  min_srx: only SRX accessions numbered at least this (recent submissions, to enrich newer kits);
  exclude.tsv: skip samples already in an earlier pool.
Needs pyarrow and h5py.
"""

import io
import multiprocessing
import os
import random
import sys
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed

import h5py
import pyarrow.parquet as pq

BLOCK = 1 << 20


class HttpFile(io.RawIOBase):
    """Read-only, seekable view of a remote file through HTTP range requests, cached per block."""

    def __init__(self, url):
        self.url, self.pos, self.blocks = url, 0, {}
        req = urllib.request.Request(url, method="HEAD")
        self.size = int(urllib.request.urlopen(req, timeout=60).headers["Content-Length"])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=0):
        self.pos = {0: offset, 1: self.pos + offset, 2: self.size + offset}[whence]
        return self.pos

    def _block(self, i):
        if i not in self.blocks:
            lo, hi = i * BLOCK, min((i + 1) * BLOCK, self.size) - 1
            req = urllib.request.Request(self.url, headers={"Range": f"bytes={lo}-{hi}"})
            self.blocks[i] = urllib.request.urlopen(req, timeout=60).read()
        return self.blocks[i]

    def readinto(self, buf):
        n = min(len(buf), self.size - self.pos)
        got = 0
        while got < n:
            i, off = divmod(self.pos + got, BLOCK)
            chunk = self._block(i)[off:off + n - got]
            buf[got:got + len(chunk)] = chunk
            got += len(chunk)
        self.pos += got
        return got

LISTS = {
    "737K-august-2016": "10x_v2_permit.txt",
    "3M-february-2018": "10x_v3_permit.txt",
    "3M-3pgex-may-2023": "10x_v4_3p_permit.txt",
    "3M-5pgex-jan-2023": "10x_v3_5p_permit.txt",
    # Fixed 8 Oct 2026: this named 10x_arc_v1_gex_label_only.txt, which holds the ARC *ATAC* list
    # (no barcode in common with the GEX list, either orientation), so no Multiome GEX sample could
    # be recovered as 737K-arc-v1. The recorded pools (results/pool.tsv, pool_recent.tsv) predate
    # the fix: of the 51 samples Arc labels multiome there, 40 are "unassigned" and 11 carry a 3'
    # list, and their hit_737K-arc-v1 column was counted against the ATAC list.
    "737K-arc-v1": "10x_arc_v1_gex_permit.txt",
}
BARCODES_PER_SAMPLE = 200
ASSIGN_AT = 0.8


def fetch_barcodes(file_path):
    try:
        url = file_path.replace("gs://", "https://storage.googleapis.com/")
        with h5py.File(io.BufferedReader(HttpFile(url), BLOCK), "r") as h:
            barcodes = [b.decode() for b in h["obs"]["_index"][:BARCODES_PER_SAMPLE]]
        return [b.split("-")[0] for b in barcodes], ""
    except Exception as e:  # recorded, not fatal: a pool row we could not read
        return [], repr(e)[:120]


def main(meta_path, wl_dir, n, seed, out_path, min_srx=0, exclude=None):
    rows = pq.read_table(meta_path).to_pylist()
    if min_srx:
        rows = [r for r in rows if r["srx_accession"].startswith("SRX") and int(r["srx_accession"][3:]) >= min_srx]
    if exclude:
        seen = {l.split("\t", 1)[0] for l in open(exclude)}
        rows = [r for r in rows if r["srx_accession"] not in seen]
    random.Random(seed).shuffle(rows)
    pool = rows[:n]
    lists = {name: set(open(f"{wl_dir}/{f}").read().split()) for name, f in LISTS.items()}

    def label(row, barcodes, err):
        if err:
            return "error", {}, 0
        fracs = {name: sum(b in s for b in barcodes) / max(len(barcodes), 1) for name, s in lists.items()}
        best = max(fracs, key=fracs.get)
        return (best if fracs[best] >= ASSIGN_AT else "unassigned"), fracs, len(barcodes)

    # Processes, not threads: h5py serialises every call behind one global lock, and our file
    # object does network I/O inside those calls. Workers return barcodes; lookups stay here.
    workers = int(os.environ.get("WORKERS", "32"))
    spawn = multiprocessing.get_context("spawn")
    with ProcessPoolExecutor(workers, mp_context=spawn) as ex, open(out_path, "w") as out:
        out.write("srx\ttech_10x\tcell_prep\tobs_count\tarc_whitelist\tn_barcodes\t"
                  + "\t".join(f"hit_{k}" for k in LISTS) + "\terror\n")
        futures = {ex.submit(fetch_barcodes, r["file_path"]): r for r in pool}
        for i, fut in enumerate(as_completed(futures), 1):
            row = futures[fut]
            barcodes, err = fut.result()
            arc_wl, fracs, nb = label(row, barcodes, err)
            out.write(f"{row['srx_accession']}\t{row['tech_10x']}\t{row['cell_prep']}\t{row['obs_count']}\t"
                      f"{arc_wl}\t{nb}\t" + "\t".join(f"{fracs.get(k, 0):.3f}" for k in LISTS) + f"\t{err}\n")
            out.flush()
            if i % 100 == 0:
                print(i, "labelled", flush=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), sys.argv[5],
         int(sys.argv[6]) if len(sys.argv) > 6 else 0, sys.argv[7] if len(sys.argv) > 7 else None)
