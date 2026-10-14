# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Chaos probes for the input paths chaos.py does not mutate: BAM read names, --emit paths, the
no-match fingerprint, IUPAC bases and named pipes. Each probe is hostile input that once aborted,
hung, leaked or misled (engineering audit and chaos run, 8 Oct 2026), and states what it must do
now. Every probe also checks the invariants: no abort (exit 134 or a signal), no panic message,
no hang (TIMEOUT seconds).

Usage: python probes.py <tool> <whitelist_dir> <work_dir> <out.tsv>
Standard library only. The base sample is pbmc_1k_v3 under $TUSSL_CHEMISTRY_DATA.
"""

import gzip
import json
import os
import random
import shutil
import struct
import subprocess
import sys
import zlib

P = f"{os.environ['TUSSL_CHEMISTRY_DATA']}/pbmc_1k_v3_fastqs"
BASE_R1 = f"{P}/pbmc_1k_v3_S1_L001_R1_001.fastq.gz"
BASE_R2 = f"{P}/pbmc_1k_v3_S1_L001_R2_001.fastq.gz"
READS = 3000  # enough for a HIGH 10x call, few enough to write fast
TIMEOUT = 30  # seconds; a probe that runs longer has hung
KMER = 12  # the fingerprint's anchor length: no 12-mer of a read may be shared
LINKER = b"GTCACTGCTAGCATCG"  # the made-up kit's fixed sequence, which may be shared
PIPELINES = ("cellranger", "starsolo", "simpleaf", "kb")


def fq(path, reads):
    with gzip.open(path, "wb") as f:
        for name, seq in reads:
            f.write(b"@" + name + b"\n" + seq + b"\n+\n" + b"I" * len(seq) + b"\n")


def bgzf(raw):
    out = b""
    for i in range(0, len(raw), 60000):
        chunk = raw[i:i + 60000]
        c = zlib.compressobj(6, zlib.DEFLATED, -15)
        data = c.compress(chunk) + c.flush()
        out += (b"\x1f\x8b\x08\x04" + b"\0" * 4 + b"\0\xff" + struct.pack("<H", 6) + b"BC"
                + struct.pack("<HH", 2, len(data) + 25) + data
                + struct.pack("<II", zlib.crc32(chunk) & 0xFFFFFFFF, len(chunk)))
    return out + bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")


def bam(path, reads):
    """An unaligned BAM: no references, one record per read, names as given."""
    text = b"@HD\tVN:1.6\tSO:unsorted\n@RG\tID:x\tPL:ILLUMINA\n"
    raw = b"BAM\1" + struct.pack("<i", len(text)) + text + struct.pack("<i", 0)
    codes = {c: i for i, c in enumerate(b"=ACMGRSVTWYHKDBN")}
    for name, seq in reads:
        packed = bytearray((len(seq) + 1) // 2)
        for i, b in enumerate(seq):
            packed[i // 2] |= codes[b] << (4 if i % 2 == 0 else 0)
        body = struct.pack("<iiBBHHHiiii", -1, -1, len(name) + 1, 0, 4680, 0, 4, len(seq), -1, -1, 0)
        body += name + b"\0" + bytes(packed) + b"\x28" * len(seq)
        raw += struct.pack("<i", len(body)) + body
    with open(path, "wb") as f:
        f.write(bgzf(raw))


def head(path, n):
    with gzip.open(path) as f:
        return [line.strip() for i, line in zip(range(4 * n), f) if i % 4 == 1]


def run(tool, wl, args, cwd=None):
    try:
        p = subprocess.run([tool, "--whitelists", wl, *args], capture_output=True, timeout=TIMEOUT,
                           cwd=cwd)
        return p.returncode, p.stdout, p.stderr
    except subprocess.TimeoutExpired:
        return "TIMEOUT", b"", b""


def shared(text, reads):
    """The first 12-mer of any read found in text, the linker's own excepted."""
    linker = {LINKER[i:i + KMER] for i in range(len(LINKER) - KMER + 1)}
    for read in reads:
        for i in range(len(read) - KMER + 1):
            piece = read[i:i + KMER]
            if piece not in linker and piece in text:
                return piece.decode()
    return None


def link_and_fingerprint(out):
    """What a no-match report shares: the JSON fingerprint, or the text report's issue link."""
    if out.lstrip().startswith(b"{"):
        return json.dumps(json.loads(out).get("fingerprint")).encode()
    return b"\n".join(line for line in out.splitlines() if b"issues/new" in line)


def probes(work, r1, r2):
    """(name, args, expect, check) for every probe; check(rc, out) -> problem or None."""
    rng = random.Random(7)
    rand = lambda n: bytes(rng.choice(b"ACGT") for _ in range(n))
    j = lambda name: os.path.join(work, name)
    good1, good2 = j("g_S1_L001_R1_001.fastq.gz"), j("g_S1_L001_R2_001.fastq.gz")
    fq(good1, [(b"r%d" % i, s) for i, s in enumerate(r1)])
    fq(good2, [(b"r%d" % i, s) for i, s in enumerate(r2)])
    exit_code = lambda want: lambda rc, out: None if rc == want else f"exit {rc}, wanted {want}"
    out = []
    # The fingerprint goes into a public issue link: no read may be in it, at any read count.
    nm1 = [rand(12) + LINKER + rand(30) for _ in range(2000)]
    nm2 = [rand(90) for _ in range(2000)]
    fq(j("nm_R1.fq.gz"), [(b"q%d" % i, s) for i, s in enumerate(nm1)])
    fq(j("nm_R2.fq.gz"), [(b"q%d" % i, s) for i, s in enumerate(nm2)])
    pair = ["--r1", j("nm_R1.fq.gz"), "--r2", j("nm_R2.fq.gz")]
    for k in (1, 2, 3, 6, 50, 2000):
        def no_leak(rc, o, k=k):
            leak = shared(link_and_fingerprint(o), nm1[:k] + nm2[:k])
            return f"read 12-mer {leak} shared" if leak else (None if rc == 0 else f"exit {rc}")
        out.append((f"fingerprint_json_reads_{k}", ["--format", "json", "--reads", str(k), *pair],
                    "no read shared", no_leak))
        out.append((f"fingerprint_text_reads_{k}", ["--reads", str(k), *pair], "no read shared", no_leak))
    # Read names are file bytes: control characters are kept visible, never asserted on.
    for tag, name in (("newline", b"a\nb"), ("nul", b"a\0b"), ("escape", b"a\x1b[31mb"), ("star", b"*")):
        bam(j(f"bam_{tag}.bam"), [(name + b"%d" % i, s) for i, s in enumerate(r2)])
        out.append((f"bam_name_{tag}", [j(f"bam_{tag}.bam")], "a verdict", exit_code(0)))
    esc = b"a\x1b]0;title\x07b"
    fq(j("esc_R1.fq.gz"), [(esc + b"%d" % i, s) for i, s in enumerate(r1)])
    fq(j("esc_R2.fq.gz"), [(esc + b"%d" % i, s) for i, s in enumerate(r2)])
    out.append(("fastq_name_escape", ["--r1", j("esc_R1.fq.gz"), "--r2", j("esc_R2.fq.gz")],
                "no raw ESC in the report",
                lambda rc, o: "raw ESC in the report" if b"\x1b" in o else exit_code(0)(rc, o)))
    # Paths a pipeline would misread are refused (exit 3) or rewritten, never asserted on.
    for d in ("nl\ndir", "a,b", "dash"):
        os.makedirs(j(d), exist_ok=True)
    for mate, src in (("R1", good1), ("R2", good2)):
        for d, stem in (("nl\ndir", "x"), ("a,b", "x"), ("dash", "-o")):
            dst = j(f"{d}/{stem}_S1_L001_{mate}_001.fastq.gz")
            if not os.path.exists(dst):
                shutil.copy(src, dst)
    files = lambda d, stem: ["--r1", j(f"{d}/{stem}_S1_L001_R1_001.fastq.gz"),
                             "--r2", j(f"{d}/{stem}_S1_L001_R2_001.fastq.gz")]
    for pipe in PIPELINES:
        for fmt in ("text", "json"):
            emit = ["--emit", pipe, "--format", fmt]
            out.append((f"emit_{pipe}_{fmt}_newline_path", emit + files("nl\ndir", "x"), "exit 3",
                        exit_code(3)))
            out.append((f"emit_{pipe}_{fmt}_same_file", emit + ["--r1", good1, "--r2", good1],
                        "exit 3", exit_code(3)))
        comma = 0 if pipe == "kb" else 3  # kb takes files as separate arguments
        out.append((f"emit_{pipe}_comma_path", ["--emit", pipe] + files("a,b", "x"), f"exit {comma}",
                    exit_code(comma)))
        def no_option(rc, o):
            bad = [w for w in o.split() if w.startswith(b"-o_")]
            return f"path read as an option: {bad[0].decode()}" if bad else exit_code(0)(rc, o)
        out.append((f"emit_{pipe}_dash_path", ["--emit", pipe, "--r1", "-o_S1_L001_R1_001.fastq.gz",
                                                "--r2", "-o_S1_L001_R2_001.fastq.gz"],
                    "no path as an option", no_option, j("dash")))
    fifo = j("fifo_R1.fq")
    if not os.path.exists(fifo):
        os.mkfifo(fifo)
    out.append(("emit_named_pipe", ["--emit", "starsolo", "--r1", fifo, "--r2", good2], "exit 3, no hang",
                exit_code(3)))
    # IUPAC codes in a prefix every read shares: a verdict, never an anchor assertion.
    for code in "RYKN":
        reads = [b"ACGT" + code.encode() + b"ACGTACGTAC" + rand(40) for _ in range(2000)]
        fq(j(f"iu_{code}_R1.fq.gz"), [(b"i%d" % i, s) for i, s in enumerate(reads)])
        for fmt in ("text", "json"):
            out.append((f"iupac_prefix_{code}_{fmt}", ["--format", fmt, "--r1", j(f"iu_{code}_R1.fq.gz"),
                                                       "--r2", j("nm_R2.fq.gz")], "a verdict", exit_code(0)))
    return [p if len(p) == 5 else (*p, None) for p in out]


def main(tool, wl, work, out_path):
    # Absolute: the dash probes run from their own directory.
    tool, wl, work = os.path.abspath(tool), os.path.abspath(wl), os.path.abspath(work)
    os.makedirs(work, exist_ok=True)
    r1, r2 = head(BASE_R1, READS), head(BASE_R2, READS)
    rows, fails = [], 0
    for name, args, expect, check, cwd in probes(work, r1, r2):
        if name == "emit_named_pipe":  # a writer that never comes: the reader must not wait
            subprocess.Popen(["sh", "-c", f"cat {BASE_R1} > '{args[3]}'"])
        rc, out, err = run(tool, wl, args, cwd)
        problem = check(rc, out) if rc != "TIMEOUT" else f"hung past {TIMEOUT} s"
        if rc in (134, -6) or b"panicked" in err:
            problem = "aborted: " + (err.decode(errors="replace").strip().splitlines() or [""])[0][:90]
        status = "FAIL" if problem else "PASS"
        fails += status == "FAIL"
        rows.append((name, expect, status, str(rc), problem or ""))
        print(f"{status}  {name:34} expect {expect:22} exit {rc} {problem or ''}", flush=True)
    with open(out_path, "w") as f:
        f.write("probe\texpect\tstatus\texit\tproblem\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"\n{len(rows) - fails}/{len(rows)} pass", flush=True)
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main(*sys.argv[1:5])
