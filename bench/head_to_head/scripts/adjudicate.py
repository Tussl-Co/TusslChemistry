# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Adjudication evidence independent of both tools, per the frozen plan in selection.tsv:
  T: what ENA's titles say (3'/5'; library type),
  N: strand of the cDNA read on ~40 nuclear genes' spliced transcripts (build_nuclear_ref.py).

Usage: python adjudicate.py <h2h results.tsv> <nuclear_genes.fa> <out.tsv> [exclude.fa]
       python adjudicate.py --calibrate <nuclear_genes.fa> [--exclude exclude.fa] <label:path>...
  exclude.fa: k-mers of this set (either strand) are dropped, so a referee built from a disjoint gene set
  shares no sequence with a tool that uses exclude.fa.
"""

import csv
import gzip
import json
import os
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor

K = 25
STEP = 5
READS = 100_000
MIN_SHARE = 0.002
MIN_RATIO = float(__import__("os").environ.get("REF_RATIO", "3.0"))  # referee bar; set B uses 2.0 (calibrated)
COMP = str.maketrans("ACGTN", "TGCAN")
FIVE = re.compile(r"5'|5’|5 prime|5prime|5-prime|\b5p\b|five prime|5'GEX|5 GEX|5' ?GEX", re.I)
THREE = re.compile(r"3'|3’|3 prime|3prime|3-prime|\b3p\b|three prime", re.I)
VDJ = re.compile(r"\b(TCR|BCR|VDJ|V\(D\)J|TRB|IGH)\b", re.I)
FEATURE = re.compile(r"\b(ADT|HTO|CITE|hashtag|antibody|feature|CRISPR|sgRNA|guide|FB)\b", re.I)


def low_complexity(k, run=8):
    return any(c * run in k for c in "ACGT")


def kmer_sets(fasta, exclude=None):
    sense, anti = kmer_sets_raw(fasta)
    # Low-complexity k-mers (A-rich transcript ends) match polyT stretches on the other strand.
    drop_lc = {k for k in sense if low_complexity(k)}
    sense, anti = sense - drop_lc, anti - {k.translate(COMP)[::-1] for k in drop_lc}
    if exclude:
        ex_s, ex_a = kmer_sets_raw(exclude)
        drop = ex_s | ex_a
        sense, anti = sense - drop, anti - drop
    return sense, anti


def kmer_sets_raw(fasta):
    sense = set()
    seq = []
    for line in open(fasta):
        if line.startswith(">"):
            if seq:
                s = "".join(seq)
                sense.update(s[i:i + K] for i in range(len(s) - K + 1))
            seq = []
        else:
            seq.append(line.strip())
    if seq:
        s = "".join(seq)
        sense.update(s[i:i + K] for i in range(len(s) - K + 1))
    anti = {k.translate(COMP)[::-1] for k in sense}
    both = sense & anti  # overlapping antisense genes and palindromes cannot vote
    return sense - both, anti - both


def strand(path, sense, anti):
    s = a = n = 0
    with gzip.open(path, "rt") as f:
        for i, line in enumerate(f):
            if i % 4 != 1:
                continue
            r = line.strip()
            n += 1
            kms = [r[j:j + K] for j in range(0, len(r) - K + 1, STEP)]
            s += sum(k in sense for k in kms) >= 2
            a += sum(k in anti for k in kms) >= 2
            if n >= READS:
                break
    n = max(n, 1)
    ss, aa = s / n, a / n
    hi, lo = max(ss, aa), min(ss, aa)
    end = ("3'" if ss > aa else "5'") if hi >= MIN_SHARE and hi >= MIN_RATIO * lo else ""
    return ss, aa, end


def titles(srx):
    url = (f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={srx}&result=read_run"
           "&fields=experiment_title,library_name,sample_title,study_title&format=json")
    try:
        r = json.loads(urllib.request.urlopen(url, timeout=60).read())[0]
    except Exception:
        return ""
    return " | ".join(r.get(k) or "" for k in ("experiment_title", "library_name", "sample_title"))


def title_end(t):
    five, three = bool(FIVE.search(t)), bool(THREE.search(t))
    return "5'" if five and not three else "3'" if three and not five else ""


def title_type(t):
    return "vdj" if VDJ.search(t) else "feature" if FEATURE.search(t) else "gex"


def main(results, fasta, out_path, exclude=None):
    sense, anti = kmer_sets(fasta, exclude)
    rows = [r for r in csv.DictReader(open(results), delimiter="\t") if r.get("condition") == "ok"]

    def one(r):
        t = titles(r["srx"])
        ss, aa, n_end = strand(os.path.expandvars(r["cdna"]), sense, anti) if r.get("cdna") else (0.0, 0.0, "")
        t_end = title_end(t)
        truth = t_end if t_end and (not n_end or n_end == t_end) else n_end if not t_end else "undetermined"
        return {**{k: r[k] for k in ("srx", "tech_10x", "arc_whitelist")}, "title": t[:200], "title_end": t_end,
                "nuc_sense": f"{100 * ss:.2f}", "nuc_anti": f"{100 * aa:.2f}", "nuc_end": n_end,
                "truth_end": truth or "undetermined", "lib_type": title_type(t)}

    with ThreadPoolExecutor(8) as ex, open(out_path, "w") as out:
        fields = ["srx", "tech_10x", "arc_whitelist", "title", "title_end", "nuc_sense", "nuc_anti", "nuc_end",
                  "truth_end", "lib_type"]
        w = csv.DictWriter(out, fields, delimiter="\t")
        w.writeheader()
        for rec in ex.map(one, rows):
            w.writerow(rec)


def calibrate(fasta, cases):
    exclude = None
    if cases and cases[0] == "--exclude":
        exclude, cases = cases[1], cases[2:]
    sense, anti = kmer_sets(fasta, exclude)
    print(f"k-mers: sense {len(sense)}, antisense {len(anti)}")
    for c in cases:
        label, path = c.split(":", 1)
        ss, aa, end = strand(path, sense, anti)
        print(f"{label:24} sense {100 * ss:5.2f}%  antisense {100 * aa:5.2f}%  -> {end or 'undetermined'}")


if __name__ == "__main__":
    if sys.argv[1] == "--calibrate":
        calibrate(sys.argv[2], sys.argv[3:])
    else:
        main(*sys.argv[1:5])
