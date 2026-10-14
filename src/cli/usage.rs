// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The help text `tchem --help` prints.

pub const USAGE: &str = "\
Names the layout of a single-cell RNA library from its first reads. It knows 25 layouts from
10x, BD Rhapsody, Parse, Singleron, SeekGene, ScaleBio, MGI, Illumina (PIPseq), SPLiT-seq,
inDrops and Drop-seq. The README lists them.

Usage: tchem --whitelists <DIR> --r1 <FILE> [--r2 <FILE>] [options]
       tchem --whitelists <DIR> [options] <FILE>...   (a lane's R1, R2, I1, I2: index reads are
                                                       recognised and set aside)
       tchem --whitelists <DIR> [options] <ACCESSION>   (a public run, e.g. SRR12345678 or an
                                                       experiment SRX...: its first reads are
                                                       streamed from ENA; needs curl)
       tchem --whitelists <DIR> --pack

Options:
  --r1 <FILE>          The file labelled R1. Swapped files are detected.
  --r2 <FILE>          The file labelled R2.
  --whitelists <DIR>   Directory holding the barcode lists (scripts/fetch_whitelists.sh).
                       [env: TUSSL_CHEMISTRY_WHITELISTS]
  --reads <N>          Reads sampled from the start of each file, 1 to 1000000. [default: 100000]
  --all                Show every candidate's evidence, not only the verdict and the top three.
  --format <FORMAT>    text (a report for people) or json (one object for pipelines).
                       [default: text]
  --emit <PIPELINE>    Print a ready command for cellranger, starsolo, simpleaf or kb instead
                       of the report (shell; JSON with --format json). Refused (exit 3) below
                       HIGH, for libraries other than gene expression, for layouts the
                       pipeline cannot take, and where it would read the files wrongly.
  --emit-anyway        With --emit: emit below HIGH too.
  --pack               Write each list's packed form (<list>.packed) beside it for faster
                       loads, then exit.
  --help               Print this help.
  --version            Print the version.
";
