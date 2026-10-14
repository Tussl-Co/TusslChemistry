# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Chaos suite: break correctly-called samples the ways archives and people break FASTQs, and
check three invariants on every run:
  1. no panic: either a verdict, or a clean "Error:" exit,
  2. never HIGH-confidence wrong: a HIGH verdict must name the base sample's true layout,
  3. deterministic: two runs print the same report.
Each mutation also states the behaviour expected of it (checked and reported separately).

Usage: python chaos.py <tool> <whitelist_dir> <work_dir> <out.tsv>
Standard library only. Base samples live under $TUSSL_CHEMISTRY_DATA (see BASES).
"""

import gzip
import os
import random
import re
import subprocess
import sys

D = os.environ["TUSSL_CHEMISTRY_DATA"]  # root of the benchmark inputs
S = f"{D}/samples"
P = f"{D}/pbmc_1k_v3_fastqs"
V = f"{D}/vendor_gold/heads"
W = f"{D}/wild/cache"
# name, R1, R2, true verdict, which file carries the barcode (0 = R1, 1 = R2), and the length the
# barcode read is trimmed to (short enough to cut into its barcode or UMI). One base per layout.
# Bases are design runs, vendor data or runs outside every held-out set: a check run on every
# commit is one the code is tuned against (bench/VALIDATION.md). The five marked "held out" are
# counted by bench/coverage today; whether they move to design is open (VALIDATION.md, Open 4).
BASES = [
    ("10x_3p_v3", f"{P}/pbmc_1k_v3_S1_L001_R1_001.fastq.gz", f"{P}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz", "10x 3' v3/v3.1", 0, 22),
    ("10x_3p_v2", f"{S}/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R1_001.fastq.gz", f"{S}/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R2_001.fastq.gz", "10x 3' v2", 0, 22),
    ("10x_5p_v2", f"{S}/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001_R1_001.fastq.gz", f"{S}/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001_R2_001.fastq.gz", "10x 5' v1/v2", 0, 22),
    ("bd_v1", f"{S}/bd_v1_smk_head/R1.fastq.gz", f"{S}/bd_v1_smk_head/R2.fastq.gz", "BD Rhapsody (original V1 beads)", 0, 40),
    ("bd_eb", f"{S}/bd_eb_wta_head/R1.fastq.gz", f"{S}/bd_eb_wta_head/R2.fastq.gz", "BD Rhapsody (Enhanced Beads)", 0, 40),
    ("parse", f"{S}/ena_SRR21740391/SRR21740391_1.fastq.gz", f"{S}/ena_SRR21740391/SRR21740391_2.fastq.gz", "Parse Evercode WT", 1, 60),
    ("splitseq", f"{S}/ena_SRR6750044/SRR6750044_1.fastq.gz", f"{S}/ena_SRR6750044/SRR6750044_2.fastq.gz", "SPLiT-seq (Rosenberg 2018)", 1, 60),
    ("dropseq", f"{S}/ena_SRR1853178/SRR1853178_1.fastq.gz", f"{S}/ena_SRR1853178/SRR1853178_2.fastq.gz", "Drop-seq", 0, 22),
    # Design runs (sources.DESIGN, or split=design in their catalogue).
    ("10x_multiome", f"{V}/pbmc_granulocyte_sorted_10k_R1.fq.gz", f"{V}/pbmc_granulocyte_sorted_10k_R2.fq.gz", "10x Multiome GEX", 0, 22),
    ("10x_flex_v1", f"{V}/10k_K562_singleplex_Multiplex_R1.fq.gz", f"{V}/10k_K562_singleplex_Multiplex_R2.fq.gz", "10x Flex v1", 0, 22),
    ("bd_v2v3", f"{V}/BD_WTANext_3070_R1.fq.gz", f"{V}/BD_WTANext_3070_R2.fq.gz", "BD Rhapsody (Enhanced V2/V3 beads)", 0, 40),
    ("quantumscale", f"{V}/ScaleBio_test_202506_quantumV1_R1.fq.gz", f"{V}/ScaleBio_test_202506_quantumV1_R2.fq.gz", "ScaleBio QuantumScale RNA", 1, 10),
    ("parse_v3", f"{D}/parse_heldout/heads/SRR34855738_R1.fq.gz", f"{D}/parse_heldout/heads/SRR34855738_R2.fq.gz", "Parse Evercode WT v3", 1, 30),
    ("mgi", f"{D}/mgi_heldout/heads/SRR25532325_R1.fq.gz", f"{D}/mgi_heldout/heads/SRR25532325_R2.fq.gz", "MGI DNBelab C4 / C-series HT", 0, 15),
    ("pipseq", f"{D}/pipseq_heldout/heads/SRR25991440_R1.fq.gz", f"{D}/pipseq_heldout/heads/SRR25991440_R2.fq.gz", "Illumina Single Cell 3' / PIPseq v3-v4", 0, 30),
    # Outside every held-out set (scripts/validate.sh runs them already).
    ("scalebio", f"{W}/SRR28437723_1.fastq.gz", f"{W}/SRR28437723_2.fastq.gz", "ScaleBio RNA", 0, 20),
    ("indrops", f"{W}/ERR10146537_2.fastq.gz", f"{W}/ERR10146537_1.fastq.gz", "inDrops", 0, 30),  # archived swapped
    # SeekGene design runs (bench/seekgene/RESEARCH.md): read to confirm the layouts.
    ("seekgene_dd3", f"{D}/seekgene_design/SRR27533148_R1.fq.gz", f"{D}/seekgene_design/SRR27533148_R2.fq.gz", "SeekGene SeekOne DD 3'", 0, 22),
    ("seekgene_dd5", f"{D}/seekgene_design/SRR23931570_R1.fq.gz", f"{D}/seekgene_design/SRR23931570_R2.fq.gz", "SeekGene SeekOne DD 5'", 0, 22),
    ("pipseq_v", f"{D}/pipseq_heldout/heads/SRR33161101_R1.fq.gz", f"{D}/pipseq_heldout/heads/SRR33161101_R2.fq.gz", "Illumina Single Cell 3' RNA Prep / PIPseq V", 0, 30),  # its design run
    # Held out: no design run is recorded for these layouts.
    ("10x_3p_v4", f"{V}/10k_hgmm_3p_gemx_Multiplex_R1.fq.gz", f"{V}/10k_hgmm_3p_gemx_Multiplex_R2.fq.gz", "10x 3' v4 (GEM-X)", 0, 22),
    ("10x_5p_v3", f"{V}/10k_hgmm_5p_gemx_Multiplex_R1.fq.gz", f"{V}/10k_hgmm_5p_gemx_Multiplex_R2.fq.gz", "10x 5' v3 (GEM-X)", 0, 22),
    ("singleron_v1", f"{D}/singleron/heads/SRR19645295_R1.fq.gz", f"{D}/singleron/heads/SRR19645295_R2.fq.gz", "Singleron GEXSCOPE V1", 0, 30),
    ("singleron_v2", f"{D}/singleron/heads/SRR36633509_R1.fq.gz", f"{D}/singleron/heads/SRR36633509_R2.fq.gz", "Singleron GEXSCOPE V2", 0, 30),
    ("singleron_v3", f"{D}/singleron/heads/SRR38275575_R1.fq.gz", f"{D}/singleron/heads/SRR38275575_R2.fq.gz", "Singleron GEXSCOPE V3", 0, 30),
]
N = 100_000
COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


def load(path, n=N):
    recs = []
    with gzip.open(path, "rb") as f:
        while len(recs) < n:
            lines = [f.readline() for _ in range(4)]
            if not lines[3]:
                break
            recs.append([l.rstrip(b"\r\n") for l in lines])
    return recs


def write(path, recs, eol=b"\n", fasta=False):
    with gzip.open(path, "wb", compresslevel=1) as f:
        for r in recs:
            if fasta:
                f.write(b">" + r[0][1:] + eol + r[1] + eol)
            else:
                f.write(eol.join(r) + eol)
    return path


def seq_map(recs, fn):
    out = []
    for h, s, p, q in recs:
        s2 = fn(s)
        out.append([h, s2, p, q[: len(s2)] if len(q) >= len(s2) else q + b"F" * (len(s2) - len(q))])
    return out


def n_storm(rate, rng):
    return lambda s: bytes(b"N"[0] if rng.random() < rate else c for c in s)


def substitutions(rate, rng):
    return lambda s: bytes(rng.choice(b"ACGT") if rng.random() < rate else c for c in s)


CONSTANT = []  # set by main: the constant read's first 17 bases, from the SeekGene list


def constant_read(s):
    """One sequence for every read, as long as the read: a 7-mer, a SeekGene list 10-mer, then a
    fixed amplicon-like filler."""
    filler = b"TTCCAGTGGCACCCGAGAATTCCAG" * 12
    return (CONSTANT[0] + filler)[: len(s)]


# Mutation -> (builder, expectation). A builder returns (r1, r2_or_None) paths, or raises for
# inputs the tool should refuse. Expectations: "same" (true verdict), "not_high", "error",
# "no_match" (no layout claimed), "swapped_same", "unswapped_same_or_lower" (files as given).
def mutations(base, a, b, other, work):
    name, _, _, truth, bc, cut = base
    rng = random.Random(7)
    j = lambda tag: os.path.join(work, f"{name}.{tag}")
    bc_recs, cd_recs = (a, b) if bc == 0 else (b, a)

    def bc_mut(fn):  # mutate the barcode read, keep the cDNA read
        m = seq_map(bc_recs, fn)
        return (m, cd_recs) if bc == 0 else (cd_recs, m)

    def pair(tag, ra, rb):
        return write(j(tag + ".R1.fq.gz"), ra), (write(j(tag + ".R2.fq.gz"), rb) if rb is not None else None)

    yield "swap", lambda: pair("swap", b, a), "swapped_same"
    # The cDNA read replaced by one constant sequence that sits on a SeekGene list entry, as a tag
    # library's amplicon read did on 9 BD runs (coverage at c157518): the list must not vouch for
    # it, so the files are not swapped and the barcode read keeps its verdict. A constant cDNA read
# carries no strand, so kits that only the strand separates (10x 3' v2 / 5' v1/v2 on 26-base
# reads) may fall to an honest ambiguity.
    def cd_mut(fn):  # mutate the cDNA read, keep the barcode read
        m = seq_map(cd_recs, fn)
        return (bc_recs, m) if bc == 0 else (m, bc_recs)

    yield "constant_mate_on_a_list", lambda: pair("const", *cd_mut(constant_read)), "unswapped_same_or_lower"
    yield "revcomp_barcode_read", lambda: pair("rc", *bc_mut(lambda s: s.translate(COMP)[::-1])), "not_high"
    yield "trim_barcode_read", lambda: pair("trim", *bc_mut(lambda s: s[:cut])), "not_high"
    yield "mix_with_other_chemistry", lambda: pair(
        "mix", [x for p in zip(a[::2], other[0][1::2]) for x in p], [x for p in zip(b[::2], other[1][1::2]) for x in p]), "not_high"
    yield "tiny_500_reads", lambda: pair("tiny", a[:500], b[:500]), "same_or_lower"
    yield "crlf_line_endings", lambda: (write(j("crlf.R1.fq.gz"), a, eol=b"\r\n"), write(j("crlf.R2.fq.gz"), b, eol=b"\r\n")), "same"
    yield "lowercase_bases", lambda: pair("lower", seq_map(a, bytes.lower), seq_map(b, bytes.lower)), "same"
    yield "n_storm_30pct", lambda: pair("nstorm", *bc_mut(n_storm(0.3, rng))), "not_high_wrong"
    yield "substitutions_3pct", lambda: pair("subs", seq_map(a, substitutions(0.03, rng)), seq_map(b, substitutions(0.03, rng))), "same_or_lower"
    yield "r2_is_copy_of_r1", lambda: pair("copy", a, a), "not_high_wrong"
    yield "barcode_read_missing", lambda: pair("cdnaonly", cd_recs, None), "no_match"
    yield "r2_half_the_reads", lambda: pair("half", a, b[: len(b) // 2]), "same"
    yield "fasta_not_fastq", lambda: (write(j("fasta.R1.fa.gz"), a, fasta=True), write(j("fasta.R2.fa.gz"), b, fasta=True)), "error"
    yield "empty_file", lambda: (write(j("empty.R1.fq.gz"), []), write(j("empty.R2.fq.gz"), [])), "error"

    def truncated_gzip():
        r1 = write(j("trunc.R1.fq.gz"), a[:5000])
        data = open(r1, "rb").read()
        open(r1, "wb").write(data[: len(data) // 2])
        return r1, write(j("trunc.R2.fq.gz"), b[:5000])
    yield "truncated_gzip", truncated_gzip, "error_or_same"

    def bad_quality():
        broken = [[h, s, p, q[:-5]] for h, s, p, q in a[:2000]]
        return write(j("badq.R1.fq.gz"), broken), write(j("badq.R2.fq.gz"), b[:2000])
    yield "quality_length_mismatch", bad_quality, "error"


def run(tool, wl, r1, r2):
    args = [tool, "--whitelists", wl, "--r1", r1] + (["--r2", r2] if r2 else [])
    p = subprocess.run(args, capture_output=True, text=True, timeout=600)
    return p.returncode, p.stdout, p.stderr


def verdict(out):
    m = re.search(r"^VERDICT\s+(.*)$", out, re.M)
    if not m:
        return None, None
    v = m.group(1).strip()
    c = re.match(r"(.+?)\s+\(confidence (\w+)\)", v)
    return (c.group(1), c.group(2)) if c else (v, None)


def judge(expect, truth, rc, out, err, swapped):
    name, conf = verdict(out)
    panicked = "panicked" in err
    clean_error = rc != 0 and "Error" in err and not panicked
    high_wrong = conf == "HIGH" and name != truth
    ok = {
        "same": name == truth,
        "same_or_lower": name == truth or (conf != "HIGH" and not panicked and rc == 0),
        "swapped_same": name == truth and swapped,
        "unswapped_same_or_lower": not swapped and rc == 0 and (name == truth or conf != "HIGH"),
        "not_high": conf != "HIGH" and rc == 0,
        "not_high_wrong": not high_wrong and rc == 0,
        "no_match": rc == 0 and (name or "").startswith("no known layout"),
        "error": clean_error,
        "error_or_same": clean_error or name == truth,
    }[expect]
    return ok, panicked, high_wrong, name, conf


def main(tool, wl, work, out_path):
    os.makedirs(work, exist_ok=True)
    with open(os.path.join(wl, "seekgene_p3cb_part2.txt")) as f:
        CONSTANT.append(b"CGAGAAT" + f.readline().strip().encode())
    loaded = {b[0]: (load(b[1]), load(b[2])) for b in BASES}
    rows, fails = [], 0
    for i, base in enumerate(BASES):
        a, b = loaded[base[0]]
        other = loaded[BASES[(i + 3) % len(BASES)][0]]
        for mname, build, expect in mutations(base, a, b, other, work):
            r1, r2 = build()
            rc, out, err = run(tool, wl, r1, r2)
            rc2, out2, _ = run(tool, wl, r1, r2)
            ok, panicked, high_wrong, name, conf = judge(expect, base[3], rc, out, err, "SWAPPED" in out)
            deterministic = out == out2 and rc == rc2
            invariant_ok = not panicked and not high_wrong and deterministic
            status = "PASS" if ok and invariant_ok else "FAIL"
            fails += status == "FAIL"
            # The scratch directory is this machine's: the record names it $WORK.
            message = (err.strip().splitlines() or [""])[-1].replace(os.path.abspath(work), "$WORK")
            detail = message[:90] if rc else f"{name} ({conf})"
            rows.append((base[0], mname, expect, status, str(panicked), str(high_wrong), str(deterministic), detail))
            print(f"{status}  {base[0]:10} {mname:26} expect {expect:15} -> {detail}", flush=True)
    with open(out_path, "w") as f:
        f.write("base\tmutation\texpect\tstatus\tpanicked\thigh_wrong\tdeterministic\tdetail\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"\n{len(rows) - fails}/{len(rows)} pass", flush=True)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main(*sys.argv[1:5])
