# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Runs every command `tchem --emit` writes, on real reads, with the pipelines themselves.

A tiny reference stands in for a genome: the 607 housekeeping transcripts in data/reference, each
its own contig, with a GTF naming its gene. STAR, kallisto (kb) and simpleaf indices are built
from it once. For each case (a sample and its labels) and pipeline, tchem emits a command, the
placeholders are filled with those indices, and the command runs on the first READS reads. A
command passes when the pipeline exits 0 and assigns reads to cells and genes; for STARsolo the
same run with the opposite strand is recorded beside it, so a wrong strand would show. Cell
Ranger cannot run without its full install and a reference, so its command is checked against
its argument parser only (an invalid chemistry is refused there before anything else).

Usage: python run.py <tchem> <whitelist_dir> <tools_dir> <data_dir> <work_dir> <out.tsv>
  tools_dir holds STAR_2.7.11b/, simpleaf-x86_64-unknown-linux-gnu/, kb venv (kb/venv),
  cellranger-10.1.0/ (bin and lib/bin only), and the alevin-fry and piscem binaries.
"""

import csv
import gzip
import json
import os
import random
import re
import shutil
import subprocess
import sys

import jsonschema

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

READS = 200_000  # the cached heads' size; fewer leaves 10k-cell samples with no cell at alevin-fry's 10-read floor
INTRON = 300
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
TRANSCRIPTS = os.path.join(REPO, "data", "reference", "housekeeping_human.fa")
# Every emitted object is checked against the published contract (needs jsonschema).
EMIT_SCHEMA = jsonschema.Draft202012Validator(
    json.load(open(os.path.join(REPO, "schema", "emit.v1.json"))))
FIELDS = ["case", "pipeline", "status", "detail", "command"]


def cases(data):
    """(name, files as given, labels, extra tchem args): 10x 3' v3, swapped, 5' v2 and v3,
    3' v4, Multiome GEX, and Drop-seq (MEDIUM, so only with --emit-anyway)."""
    p = f"{data}/pbmc_1k_v3_fastqs/pbmc_1k_v3_S1_L001"
    f5 = f"{data}/samples/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001"
    vg = f"{data}/vendor_gold/heads"
    ds = f"{data}/samples/ena_SRR1853178/SRR1853178"
    return [
        ("10x 3' v3", [f"{p}_R1_001.fastq.gz", f"{p}_R2_001.fastq.gz"], ["--r1", "--r2"], []),
        ("10x 3' v3, swapped", [f"{p}_R2_001.fastq.gz", f"{p}_R1_001.fastq.gz"], ["--r1", "--r2"], []),
        ("10x 5' v1/v2", [f"{f5}_R1_001.fastq.gz", f"{f5}_R2_001.fastq.gz"], ["--r1", "--r2"], []),
        ("10x 5' v3", [f"{vg}/10k_hgmm_5p_gemx_Multiplex_R1.fq.gz", f"{vg}/10k_hgmm_5p_gemx_Multiplex_R2.fq.gz"],
         ["--r1", "--r2"], []),
        ("10x 3' v4", [f"{vg}/10k_hgmm_3p_gemx_Multiplex_R1.fq.gz", f"{vg}/10k_hgmm_3p_gemx_Multiplex_R2.fq.gz"],
         ["--r1", "--r2"], []),
        ("10x Multiome GEX", [f"{vg}/human_brain_3k_R1.fq.gz", f"{vg}/human_brain_3k_R2.fq.gz"],
         ["--r1", "--r2"], []),
        ("Drop-seq", [f"{ds}_1.fastq.gz", f"{ds}_2.fastq.gz"], ["--r1", "--r2"], ["--emit-anyway"]),
    ]


def subset(src, dst):
    """The first READS records of a gzipped FASTQ, under the same file name, so Illumina names
    (which Cell Ranger reads) survive."""
    with gzip.open(src, "rt") as fin, gzip.open(dst, "wt") as fout:
        for _ in range(READS * 4):
            line = fin.readline()
            if not line:
                break
            fout.write(line)


def reference(work):
    """The tiny genome and GTF: one contig per transcript, one exon spanning it."""
    fa, gtf, t2g = (os.path.join(work, n) for n in ("genome.fa", "genes.gtf", "t2g.tsv"))
    if os.path.exists(gtf):
        return fa, gtf, t2g
    records, name, seq = [], None, []
    for line in open(TRANSCRIPTS):
        if line.startswith(">"):
            if name:
                records.append((name, "".join(seq)))
            name, seq = line[1:].split(), []
        else:
            seq.append(line.strip())
    records.append((name, "".join(seq)))
    rng = random.Random(7)  # the introns are fixed, so every build is the same reference
    with open(fa, "w") as f, open(gtf, "w") as g, open(t2g, "w") as m:
        for (tx, gene), seq in records:
            # Two exons around a 300-base intron (GT...AG), as simpleaf's default index needs.
            half = len(seq) // 2
            intron = "GT" + "".join(rng.choice("ACGT") for _ in range(INTRON - 4)) + "AG"
            f.write(f">{tx}\n{seq[:half]}{intron}{seq[half:]}\n")
            attrs = f'gene_id "{gene}"; transcript_id "{tx}"; gene_name "{gene}";'
            end = len(seq) + INTRON
            for kind, start, stop in [("gene", 1, end), ("transcript", 1, end), ("exon", 1, half),
                                      ("exon", half + INTRON + 1, end)]:
                g.write(f"{tx}\tbench\t{kind}\t{start}\t{stop}\t.\t+\t.\t{attrs}\n")
            m.write(f"{tx}\t{gene}\t{gene}\n")
    return fa, gtf, t2g


def run(argv, cwd, env=None, timeout=1800):
    p = subprocess.run(argv, cwd=cwd, capture_output=True, text=True, timeout=timeout,
                       env={**os.environ, **(env or {})})
    return p.returncode, p.stdout + p.stderr


def indices(tools, work):
    """STAR, kallisto and simpleaf indices of the tiny reference, built once."""
    fa, gtf, t2g = reference(work)
    star = os.path.join(work, "star_index")
    if not os.path.exists(os.path.join(star, "SA")):
        os.makedirs(star, exist_ok=True)
        rc, out = run([f"{tools}/STAR_2.7.11b/Linux_x86_64_static/STAR", "--runMode", "genomeGenerate",
                       "--genomeDir", star, "--genomeFastaFiles", fa, "--sjdbGTFfile", gtf,
                       "--genomeSAindexNbases", "8", "--runThreadN", "8"], work)
        invariant(rc == 0, out[-2000:])
    kallisto = os.path.join(work, "kallisto.idx")
    kb_bin = f"{tools}/kb/venv/lib/python3.14/site-packages/kb_python/bins/linux/kallisto/kallisto"
    if not os.path.exists(kallisto):
        rc, out = run([kb_bin, "index", "-i", kallisto, fa], work)
        invariant(rc == 0, out[-2000:])
    home = os.path.join(work, "af_home")
    simpleaf = f"{tools}/simpleaf-x86_64-unknown-linux-gnu/simpleaf"
    af_index = os.path.join(work, "simpleaf_index")
    if not os.path.exists(os.path.join(af_index, "index")):
        rc, out = run([simpleaf, "set-paths", "--piscem", f"{tools}/piscem", "--alevin-fry",
                       f"{tools}/alevin-fry"], work, {"ALEVIN_FRY_HOME": home})
        invariant(rc == 0, out[-2000:])
        rc, out = run([simpleaf, "index", "--fasta", fa, "--gtf", gtf, "--rlen", "91",
                       "--output", af_index, "--threads", "8"], work, {"ALEVIN_FRY_HOME": home})
        invariant(rc == 0, out[-2000:])
    # Cell Ranger checks the reference exists while parsing; an empty directory gets past that.
    os.makedirs(os.path.join(work, "no_reference"), exist_ok=True)
    return {"STAR_INDEX": star, "KB_INDEX": kallisto, "KB_T2G": t2g,
            "SIMPLEAF_INDEX": os.path.join(af_index, "index"),
            "CELLRANGER_REF": os.path.join(work, "no_reference")}, home


def tools_path(tools):
    """PATH with every pipeline's program first."""
    dirs = [f"{tools}/STAR_2.7.11b/Linux_x86_64_static", f"{tools}/kb/venv/bin",
            f"{tools}/simpleaf-x86_64-unknown-linux-gnu", f"{tools}/cellranger-10.1.0/bin"]
    return ":".join(dirs + [os.environ["PATH"]])


def outcome(pipeline, cwd, rc, out):
    """pass/fail and the figures that show the pipeline counted reads into cells and genes."""
    if pipeline == "cellranger":
        parsed = "invalid value" not in out and "unexpected argument" not in out
        return ("parsed" if parsed else "fail"), out.strip().splitlines()[-1][:150] if out else ""
    if rc != 0:
        return "fail", out.strip()[-300:]
    if pipeline == "starsolo":
        summary = dict(l.split(",", 1) for l in open(f"{cwd}/Solo.out/Gene/Summary.csv"))
        return "pass", (f"valid barcodes {summary['Reads With Valid Barcodes'].strip()}, "
                        f"reads in genes {summary['Reads Mapped to Gene: Unique Gene'].strip()}")
    if pipeline == "kb":
        info = json.load(open(f"{cwd}/kb_count/run_info.json"))
        inspect = json.load(open(f"{cwd}/kb_count/inspect.json"))
        return "pass", (f"pseudoaligned {info['p_pseudoaligned']}%, "
                        f"barcodes on list {inspect.get('percentageReadsOnOnlist', 'n/a')}%")
    info = json.load(open(f"{cwd}/simpleaf_quant/af_map/map_info.json"))
    return "pass", f"mapped {info.get('percent_mapped', 'n/a')}%"


def cellranger_names(tools, work):
    """Every Cell Ranger chemistry the layout table holds, against Cell Ranger's own parser: an
    unknown name is refused there ("invalid value for '--chemistry'") before anything runs."""
    table = open(os.path.join(REPO, "src", "layouts", "tenx.rs")).read()
    names = re.findall(r'cellranger_chemistry: Some\("([^"]+)"\)', table)
    invariant(names, "the 10x layouts name their Cell Ranger chemistries")
    cr = f"{tools}/cellranger-10.1.0/bin/cellranger"
    rows = []
    for name in names:
        rc, out = run([cr, "count", "--id", "check", "--create-bam", "true", "--chemistry", name],
                      work)
        refused = "for '--chemistry" in out and "invalid value" in out
        rows.append({"case": f"layout table: {name}", "pipeline": "cellranger",
                     "status": "fail" if refused else "parsed",
                     "detail": "refused by the parser" if refused else "accepted by the parser"})
    return rows


def main(tchem, whitelists, tools, data, work, out_path):
    os.makedirs(work, exist_ok=True)
    values, home = indices(tools, work)
    env = {**values, "ALEVIN_FRY_HOME": home, "PATH": tools_path(tools)}
    rows = []
    for name, files, labels, extra in cases(data):
        slug = "".join(c if c.isalnum() else "_" for c in name)
        reads = os.path.join(work, "reads", slug)
        os.makedirs(reads, exist_ok=True)
        local = []
        for f in files:
            dst = os.path.join(reads, os.path.basename(f))
            if not os.path.exists(dst):
                subset(f, dst)
            local.append(dst)
        given = [x for pair in zip(labels, local) for x in pair]
        for pipeline in ("cellranger", "starsolo", "simpleaf", "kb"):
            cwd = os.path.join(work, "runs", slug, pipeline)
            shutil.rmtree(cwd, ignore_errors=True)
            os.makedirs(cwd)
            rc, emitted = run([tchem, "--whitelists", whitelists, "--emit", pipeline, "--format",
                               "json", *extra, *given], cwd)
            if rc != 0:
                reason = next((l for l in emitted.splitlines() if l.startswith("Error")), emitted)
                rows.append({"case": name, "pipeline": pipeline, "status": f"refused (exit {rc})",
                             "detail": reason.removeprefix("Error: ")})
                continue
            command = json.loads(emitted[:emitted.rindex("}") + 1])
            schema_errors = list(EMIT_SCHEMA.iter_errors(command))
            invariant(not schema_errors,  # the message is built either way: no [0] on an empty list
                      f"{name} {pipeline}: {[e.message for e in schema_errors[:1]]}")
            argv = [values.get(a[2:-1], a) if a.startswith("${") else a for a in command["argv"]]
            rc, out = run(argv, cwd, env)
            status, detail = outcome(pipeline, cwd, rc, out)
            if pipeline == "simpleaf" and "contains no output cells" in out:
                # The tiny reference leaves a 10k-cell sample with no barcode at alevin-fry's
                # 10-read floor; the bench (never the emitted command) lowers it to 1 and says so.
                shutil.rmtree(os.path.join(cwd, "simpleaf_quant"), ignore_errors=True)
                rc, out = run([*argv, "--min-reads", "1"], cwd, env)
                status, detail = outcome(pipeline, cwd, rc, out)
                status += " (bench: --min-reads 1)"
            rows.append({"case": name, "pipeline": pipeline, "status": status, "detail": detail,
                         "command": " ".join(command["argv"])})
            if pipeline == "starsolo" and status == "pass":
                flip = {"Forward": "Reverse", "Reverse": "Forward"}
                i = argv.index("--soloStrand")
                opposite = [*argv[:i + 1], flip[argv[i + 1]], *argv[i + 2:]]
                ocwd = cwd + "_opposite_strand"
                shutil.rmtree(ocwd, ignore_errors=True)
                os.makedirs(ocwd)
                rc, out = run(opposite, ocwd, env)
                status, _ = outcome(pipeline, ocwd, rc, out)
                genes = lambda d: float(dict(l.split(",", 1) for l in open(
                    f"{d}/Solo.out/Gene/Summary.csv"))["Reads Mapped to Gene: Unique Gene"])
                emitted_share, opposite_share = genes(cwd), genes(ocwd)
                # The emitted strand should put at least twice as many reads in genes.
                verdict = "emitted strand wins" if emitted_share >= 2 * opposite_share else "inconclusive"
                rows.append({"case": name, "pipeline": "starsolo strand control",
                             "status": verdict if status == "pass" else status,
                             "detail": f"reads in genes: emitted strand {emitted_share}, "
                                       f"opposite {opposite_share}"})
    rows += cellranger_names(tools, work)
    with open(out_path, "w", newline="") as f:
        version = subprocess.run([tchem, "--version"], capture_output=True, text=True).stdout.strip()
        f.write(f"# Emitted commands run by the pipelines themselves ({version}; first {READS} reads;\n"
                f"# STAR 2.7.11b, kb 0.30.2, simpleaf 0.30.0, Cell Ranger 10.1.0 argument parser).\n")
        w = csv.DictWriter(f, FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        # Paths on this machine say nothing to a reader: name the roots instead.
        portable = lambda v: str(v).replace(work, "$WORK").replace(tools, "$TOOLS").replace(data, "$TUSSL_CHEMISTRY_DATA")
        w.writerows({k: portable(v) for k, v in r.items()} for r in rows)
    for r in rows:
        print(f"{r['case']:20} {r['pipeline']:38} {r['status']:18} {r['detail'][:110]}")


if __name__ == "__main__":
    if len(sys.argv) != 7:
        sys.exit(__doc__)
    main(*sys.argv[1:])
