# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Checks every head file's FASTQ structure: 4-line records, '@' header, '+' line,
sequence and quality of equal length, sequence in ACGTN. Prints one line per bad file.
Usage: python scan_heads.py <out.tsv> <file>..."""
import gzip, sys, re
from concurrent.futures import ProcessPoolExecutor
SEQ = re.compile(rb"^[ACGTNacgtn.]*$")
def check(p):
    try:
        with gzip.open(p, "rb") as f:
            n = 0
            while True:
                h = f.readline()
                if not h: return p, n, ""
                s, plus, q = f.readline(), f.readline(), f.readline()
                n += 1
                if not h.startswith(b"@"): return p, n, "header"
                if not plus.startswith(b"+"): return p, n, "plus line"
                s, q = s.rstrip(b"\n"), q.rstrip(b"\n")
                if len(s) != len(q): return p, n, "seq/qual length"
                if not SEQ.match(s): return p, n, "sequence alphabet"
    except Exception as e:
        return p, -1, "unreadable: " + repr(e)[:80]
if __name__ == "__main__":
    files = sys.argv[2:]
    with ProcessPoolExecutor(8) as pool, open(sys.argv[1], "w") as out:
        out.write("file\trecord\tproblem\n")
        bad = 0
        for p, n, why in pool.map(check, files, chunksize=8):
            if why:
                bad += 1; out.write(f"{p}\t{n}\t{why}\n")
    print(len(files), "files;", bad, "bad")
