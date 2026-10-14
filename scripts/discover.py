# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Reference-free structure discovery for a barcode read: finds the commonest fixed k-mer (a linker),
aligns reads on it, and labels each position relative to it as fixed (one base dominates), random
(UMI-like, near 2 bits) or concentrated (barcode-like: few distinct k-mers cover most reads).

Adapter read-through (Nextera, TruSeq) is never an anchor. After one structure is described, the
reads that carry its anchor are set aside and the search repeats, so a library mixing two bead
designs shows both (up to three structures, each in at least 10% of reads).

Usage: python discover.py <fastq.gz> [reads] [anchor_k]
"""

import collections
import gzip
import math
import sys


def reads(path, n):
    out = []
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 == 1:
                out.append(line.strip())
                if len(out) >= n:
                    break
    return out


def entropy(counter):
    total = sum(v for k, v in counter.items() if k in "ACGT")
    return -sum(v / total * math.log2(v / total) for k, v in counter.items() if k in "ACGT" and v) if total else 0.0


ADAPTERS = ("CTGTCTCTTATACACATCT", "AGATGTGTATAAGAGACAG", "AGATCGGAAGAGC", "GCTCTTCCGATCT")


def is_adapter(kmer):
    return any(kmer in a or a in kmer for a in ADAPTERS) or any(kmer[:8] in a for a in ADAPTERS)


def main(path, n=50_000, k=12):
    rs = reads(path, n)
    total = len(rs)
    for structure in range(1, 4):
        per_read = collections.Counter()
        for r in rs:
            # Anchor candidates exclude homopolymer runs and adapter read-through.
            per_read.update({r[i:i + k] for i in range(len(r) - k + 1) if len(set(r[i:i + k])) > 2 and not is_adapter(r[i:i + k])})
        if not per_read:
            return
        anchor, hits = per_read.most_common(1)[0]
        if hits < 0.10 * total:
            return
        print(f"\n=== structure {structure}: anchor {anchor} in {100 * hits / total:.1f}% of all reads")
        describe(rs, anchor)
        rs = [r for r in rs if anchor not in r]


def describe(rs, anchor):
    pos = collections.Counter(r.find(anchor) for r in rs if anchor in r)
    print(f"start positions {pos.most_common(5)}")
    aligned = [(r, r.find(anchor)) for r in rs if anchor in r]
    lo = -max(p for _, p in aligned)
    hi = max(len(r) - p for r, p in aligned)
    rows = []
    for rel in range(lo, min(hi, 70)):
        col = collections.Counter(r[p + rel] for r, p in aligned if 0 <= p + rel < len(r))
        if sum(col.values()) < len(aligned) * 0.5:
            continue
        base, top = col.most_common(1)[0]
        frac, h = top / sum(col.values()), entropy(col)
        kind = "F" if frac > 0.85 else "." if h > 1.85 else "b" if h < 1.85 and frac < 0.6 else "?"
        rows.append((rel, base, kind, frac, h))
    print("relative to anchor start (F fixed, . random/UMI, b barcode-like, ? mixed):")
    print("  " + "".join(base if kind == "F" else kind for _, base, kind, _, _ in rows))
    # Concentration of 8-mers in each non-fixed run: barcodes repeat, UMIs do not.
    run, runs = [], []
    for rel, _, kind, _, _ in rows + [(None, None, "F", 0, 0)]:
        if kind != "F":
            run.append(rel)
        elif run:
            runs.append((run[0], run[-1]))
            run = []
    for a, b in runs:
        ln = b - a + 1
        seqs = collections.Counter(r[p + a:p + b + 1] for r, p in aligned if 0 <= p + a and p + b < len(r))
        tot = sum(seqs.values())
        cover = sum(c for _, c in seqs.most_common(384)) / tot if tot else 0
        print(f"  run {a:+d}..{b:+d} ({ln} bp): {len(seqs)} distinct, top-384 cover {100 * cover:.1f}% "
              f"-> {'barcode-like' if cover > 0.5 else 'random (UMI/cDNA)'}")


if __name__ == "__main__":
    main(sys.argv[1], *(int(x) for x in sys.argv[2:4]))
