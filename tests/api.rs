// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The library as a mapper would use it: barcode lists from a directory, reads pushed in from
//! memory, a structured verdict back. Synthetic data: 10x 3' v3 reads whose barcodes come from a
//! small made-up v3 list, with every other list made of different barcodes.

// Test code: a panic is how a test reports failure.
#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)]

use std::path::{Path, PathBuf};
use tussl_chemistry::{
    Arg, Confidence, DetectRefusal, EmitRefusal, Input, Pipeline, ReadKind, Reads, SCHEMA_VERSION,
    Verdict, Whitelists, detect, packed_path,
};

const READS: usize = 5_000;
const CDNA_LENGTH: usize = 90;

/// A deterministic stream of bases.
struct Bases(u64);

impl Bases {
    fn next(&mut self, length: usize) -> Vec<u8> {
        (0..length)
            .map(|_| {
                self.0 = self.0.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                b"ACGT"[usize::try_from(self.0 >> 62).unwrap()]
            })
            .collect()
    }
}

/// A test's own directory under the system temporary directory, removed when the test ends, so
/// that `cargo test` leaves nothing behind.
struct Scratch(PathBuf);

impl Scratch {
    fn new(name: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("tchem-api-{name}-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        Self(dir)
    }
}

impl std::ops::Deref for Scratch {
    type Target = Path;
    fn deref(&self) -> &Path {
        &self.0
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        // A failure to remove is not the test's failure: the directory is only temporary.
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A whitelist directory: every list the layouts name, of the right barcode lengths, with the
/// 10x v3 list holding `v3`.
fn whitelist_dir(name: &str, v3: &[Vec<u8>], bases: &mut Bases) -> Scratch {
    let dir = Scratch::new(name);
    let write = |file: &str, barcodes: &[Vec<u8>]| {
        let text: Vec<String> =
            barcodes.iter().map(|b| String::from_utf8(b.clone()).unwrap()).collect();
        std::fs::write(dir.join(file), text.join("\n") + "\n").unwrap();
    };
    write("10x_v3_permit.txt", v3);
    for (file, length) in [
        ("10x_v2_permit.txt", 16),
        ("10x_v4_3p_permit.txt", 16),
        ("10x_v3_5p_permit.txt", 16),
        ("bd_cls1.txt", 9),
        ("bd_cls2.txt", 9),
        ("bd_cls3.txt", 9),
        ("splitseq_round1.txt", 8),
        ("splitseq_round2.txt", 8),
        ("splitseq_round3.txt", 8),
        ("singleron_v1_bc.txt", 8),
        ("singleron_bc1.txt", 9),
        ("singleron_bc2.txt", 9),
        ("singleron_bc3.txt", 9),
        ("10x_arc_v1_gex_permit.txt", 16),
        ("parse_v3_rounds23.txt", 8),
        ("parse_v3_round1.txt", 8),
        ("mgi_c_series_bc.txt", 10),
        ("bd_v2_cls1.txt", 9),
        ("bd_v2_cls2.txt", 9),
        ("bd_v2_cls3.txt", 9),
        ("scale_quantum_rt.txt", 8),
        ("seekgene_p3cb_part2.txt", 10),
    ] {
        let other: Vec<Vec<u8>> = (0..96).map(|_| bases.next(length)).collect();
        write(file, &other);
    }
    dir
}

/// 10x 3' v3 reads: R1 is a listed 16-base barcode and a 12-base random UMI, R2 random cDNA.
fn tenx_v3_reads(v3: &[Vec<u8>], bases: &mut Bases) -> (Reads, Reads) {
    let mut r1 = Reads::with_capacity(READS, READS * 28);
    let mut r2 = Reads::with_capacity(READS, READS * CDNA_LENGTH);
    for index in 0..READS {
        let mut barcode_read = v3[index % v3.len()].clone();
        barcode_read.extend(bases.next(12));
        r1.push(&barcode_read);
        r2.push(&bases.next(CDNA_LENGTH));
    }
    (r1, r2)
}

#[test]
fn detects_a_layout_from_reads_in_memory() {
    let mut bases = Bases(11);
    let v3: Vec<Vec<u8>> = (0..500).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("detect", &v3, &mut bases);
    let whitelists = Whitelists::load(&dir).unwrap();
    assert_eq!(whitelists.dir(), &*dir);
    assert_eq!(whitelists.lists().len(), 23);
    let (r1, r2) = tenx_v3_reads(&v3, &mut bases);

    let report = detect(&r1, Some(&r2), &whitelists).unwrap();
    let Verdict::Match { layout, confidence } = report.verdict() else {
        panic!("expected a match, got {:?}", report.verdict())
    };
    assert_eq!(layout.name, "10x 3' v3/v3.1");
    assert_eq!(layout.geometry, "1{b[16]u[12]x:}2{r:}");
    // The layout carries its held-out record (data/validation.tsv): 10x 3' v3 has one.
    assert!(layout.validation.held_out_high_calls > 0, "{:?}", layout.validation);
    assert_eq!(layout.validation.wrong, 0);
    assert_eq!(layout.barcode_input, Input::First);
    assert_eq!(layout.barcode_lists, vec![dir.join("10x_v3_permit.txt")]);
    assert_ne!(confidence, Confidence::Low, "a clean listed match is not LOW");
    assert!(!report.swapped());
    assert_eq!(report.layout(), Some(layout));
    // The directory has no Flex list (it ships only with Cell Ranger): loading goes on without
    // it, and the report names the layout it could not check.
    assert_eq!(report.not_checked().collect::<Vec<_>>(), vec!["10x Flex v1"]);

    // The same reads given the wrong way round: detected, and the barcode read is the second.
    let swapped = detect(&r2, Some(&r1), &whitelists).unwrap();
    assert!(swapped.swapped());
    assert_eq!(swapped.layout().unwrap().barcode_input, Input::Second);

    // The command line's two outputs come from the same report.
    let text = report.to_text(false);
    assert!(text.contains("VERDICT  10x 3' v3/v3.1"), "{text}");
    let json = report.to_json("R1.fq.gz".as_ref(), Some("R2.fq.gz".as_ref())).unwrap();
    assert!(json.contains(&format!("\"schema_version\": \"{SCHEMA_VERSION}\"")), "{json}");
    assert!(json.contains("\"verdict\": \"match\""), "{json}");
    assert!(json.contains("\"barcode_read\": \"R1.fq.gz\""), "{json}");
    assert!(json.contains("\"not_checked\": [\n    \"10x Flex v1\"\n  ]"), "{json}");
    assert!(text.contains("not checked, list not fetched: 10x Flex v1"), "{text}");

    // An archive's model is the caller's text: a newline in it, or a second note, aborted.
    let mut report = report;
    report.note_archive_model("Illumina\nNovaSeq 6000");
    report.note_archive_model("Illumina NovaSeq 6000");
    let json = report.to_json("R1.fq.gz".as_ref(), Some("R2.fq.gz".as_ref())).unwrap();
    assert!(json.contains("Illumina NovaSeq 6000"), "{json}");
}

#[test]
fn refuses_empty_samples_and_packs_lists() {
    let mut bases = Bases(23);
    let v3: Vec<Vec<u8>> = (0..50).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("pack", &v3, &mut bases);
    let packed = Whitelists::pack(&dir).unwrap();
    assert!(packed_path(&dir.join("10x_v3_permit.txt")).exists());
    let loaded = Whitelists::load(&dir).unwrap();
    assert_eq!(loaded.barcode_count(), packed.barcode_count());

    let empty = Reads::with_capacity(0, 0);
    let (r1, _) = tenx_v3_reads(&v3, &mut bases);
    assert_eq!(detect(&empty, None, &loaded).err(), Some(DetectRefusal::NoReads(Input::First)));
    assert_eq!(
        detect(&r1, Some(&empty), &loaded).err(),
        Some(DetectRefusal::NoReads(Input::Second))
    );
}

/// A long-read run in two files (ENA serves some that way) goes to the long-read detector, and a
/// missing list directory names the command that fetches the lists (user test, 12 Oct 2026).
#[test]
fn long_reads_in_two_files_and_a_missing_list_directory() {
    let mut bases = Bases(31);
    let v3: Vec<Vec<u8>> = (0..50).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("long", &v3, &mut bases);
    let whitelists = Whitelists::load(&dir).unwrap();
    let mut r1 = Reads::with_capacity(READS, READS * 600);
    let mut r2 = Reads::with_capacity(READS, READS * 600);
    for _ in 0..READS {
        r1.push(&bases.next(600));
        r2.push(&bases.next(600));
    }
    assert_eq!(detect(&r1, None, &whitelists).unwrap().read_kind(), ReadKind::Long);
    assert_eq!(detect(&r1, Some(&r2), &whitelists).unwrap().read_kind(), ReadKind::Long);

    let missing = Whitelists::load(&dir.join("not-there")).err().unwrap().to_string();
    assert!(missing.starts_with("cannot open barcode list"), "{missing}");
    assert!(missing.contains("scripts/fetch_whitelists.sh"), "{missing}");
}

/// `--emit` from the library: the command carries the layout's positions, list and strand, files
/// in the pipeline's order (also when swapped), and refuses what the pipeline would read wrongly.
#[test]
fn emits_commands_and_refuses_what_would_be_read_wrongly() {
    let mut bases = Bases(13);
    let v3: Vec<Vec<u8>> = (0..500).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("emit", &v3, &mut bases);
    let whitelists = Whitelists::load(&dir).unwrap();
    let (r1, r2) = tenx_v3_reads(&v3, &mut bases);
    let report = detect(&r1, Some(&r2), &whitelists).unwrap();
    let words = |pipeline, a: &str, b: &str| -> Result<Vec<String>, EmitRefusal> {
        let emitted = report.emit(pipeline, a.as_ref(), Some(b.as_ref()), true)?;
        Ok(emitted
            .argv
            .iter()
            .map(|arg| match arg {
                Arg::Text(t) => t.clone(),
                Arg::Var(p) => format!("${{{}}}", p.name),
            })
            .collect())
    };
    let star = words(Pipeline::Starsolo, "x_R1.fq.gz", "x_R2.fq.gz").unwrap();
    let after = |flag: &str| star.iter().position(|w| w == flag).map(|i| star[i + 1].clone());
    assert_eq!(after("--readFilesIn").as_deref(), Some("x_R2.fq.gz"), "cDNA first: {star:?}");
    assert_eq!(star.last().map(String::as_str), Some("Gene"));
    assert_eq!(after("--soloUMIstart").as_deref(), Some("17"));
    assert_eq!(after("--soloUMIlen").as_deref(), Some("12"));
    assert_eq!(after("--readFilesCommand").as_deref(), Some("gunzip -c"));
    let kb = words(Pipeline::Kb, "x_R1.fq.gz", "x_R2.fq.gz").unwrap();
    assert!(kb.contains(&"0,0,16:0,16,28:1,0,0".to_owned()), "{kb:?}");
    // Cell Ranger reads files by name: archive names are refused, Illumina names give a command.
    assert!(matches!(
        words(Pipeline::CellRanger, "x_R1.fq.gz", "x_R2.fq.gz"),
        Err(EmitRefusal::NotIlluminaNamed(_))
    ));
    let cr = words(
        Pipeline::CellRanger,
        "run/s_S1_L001_R1_001.fastq.gz",
        "run/s_S1_L001_R2_001.fastq.gz",
    )
    .unwrap();
    assert!(cr.windows(2).any(|w| w == ["--chemistry", "SC3Pv3"]), "{cr:?}");
    assert!(cr.windows(2).any(|w| w == ["--transcriptome", "${CELLRANGER_REF}"]), "{cr:?}");
    // Swapped: STARsolo still puts the cDNA first; Cell Ranger refuses.
    let swapped = detect(&r2, Some(&r1), &whitelists).unwrap();
    let emitted = swapped
        .emit(Pipeline::Starsolo, "x_R1.fq.gz".as_ref(), Some("x_R2.fq.gz".as_ref()), true)
        .unwrap();
    let first_read =
        emitted.argv.iter().skip_while(|a| **a != Arg::Text("--readFilesIn".into())).nth(1);
    assert_eq!(first_read, Some(&Arg::Text("x_R1.fq.gz".into())), "the labelled R1 holds cDNA");
    assert_eq!(
        swapped.emit(
            Pipeline::CellRanger,
            "s_S1_L001_R1_001.fastq.gz".as_ref(),
            Some("s_S1_L001_R2_001.fastq.gz".as_ref()),
            true
        ),
        Err(EmitRefusal::SwappedForCellRanger)
    );
    // One file: nothing to count against.
    assert_eq!(
        report.emit(Pipeline::Kb, "x_R1.fq.gz".as_ref(), None, true).map(|_| ()),
        Err(EmitRefusal::OneFile)
    );
    // Chaos probes: a newline in a path, or one file given twice, aborted every pipeline.
    for pipeline in Pipeline::ALL {
        // An empty path, which a library caller can pass, reached an assertion too.
        for (a, b) in [
            ("nl\ndir/x_R1.fq.gz", "x_R2.fq.gz"),
            ("x_R1.fq.gz", "x_R1.fq.gz"),
            ("", "x_R2.fq.gz"),
            ("x_R1.fq.gz", ""),
        ] {
            let refused = words(pipeline, a, b);
            assert!(matches!(refused, Err(EmitRefusal::InputFormat(_))), "{a:?} {refused:?}");
        }
    }
    // Chaos probes: a comma split one path in two where the option takes a list, and a path
    // starting with '-' became an option to kb.
    for pipeline in [Pipeline::Starsolo, Pipeline::CellRanger, Pipeline::Simpleaf] {
        let (a, b) = ("a,b/s_S1_L001_R1_001.fastq.gz", "a,b/s_S1_L001_R2_001.fastq.gz");
        let refused = words(pipeline, a, b);
        assert!(matches!(refused, Err(EmitRefusal::InputFormat(_))), "{pipeline:?} {refused:?}");
    }
    // A named pipe gives its reads once: --emit hung opening it again to read its format.
    let fifo = std::env::temp_dir().join(format!("tchem_fifo_{}", std::process::id()));
    let made = std::process::Command::new("mkfifo").arg(&fifo).status();
    if made.is_ok_and(|s| s.success()) {
        let refused = words(Pipeline::Starsolo, fifo.to_str().unwrap(), "x_R2.fq.gz");
        std::fs::remove_file(&fifo).unwrap();
        assert!(matches!(refused, Err(EmitRefusal::InputFormat(_))), "{refused:?}");
    }
    let kb = words(Pipeline::Kb, "-o_R1.fq.gz", "-o_R2.fq.gz").unwrap();
    assert!(kb.iter().all(|w| !w.starts_with("-o_")), "{kb:?}");
    assert!(kb.contains(&"./-o_R1.fq.gz".to_owned()), "{kb:?}");
    // The shell form fills placeholders from the environment, and the JSON names its contract.
    let shell =
        report.emit(Pipeline::Kb, "a b.fq.gz".as_ref(), Some("c.fq.gz".as_ref()), true).unwrap();
    assert!(shell.to_shell().contains("'a b.fq.gz'"), "{}", shell.to_shell());
    assert!(shell.to_json().contains("\"emit_schema_version\": \"1.0\""));
}

/// The no-match fingerprint goes into a public issue link, so it may hold no read. Chaos probe:
/// over 2 reads its sketch spelled the reads out base for base, the cDNA read included. (One
/// or two synthetic reads here match a layout at LOW, which offers no link.)
#[test]
fn a_fingerprint_over_few_reads_holds_no_read() {
    let mut bases = Bases(19);
    let v3: Vec<Vec<u8>> = (0..500).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("fingerprint", &v3, &mut bases);
    let whitelists = Whitelists::load(&dir).unwrap();
    let linker = b"GTCACTGCTAGCATCG";
    let sample = |n: usize, bases: &mut Bases| {
        let (mut r1, mut r2) = (Reads::with_capacity(n, n * 58), Reads::with_capacity(n, n * 90));
        let mut seen = Vec::new();
        for _ in 0..n {
            let read = [bases.next(12), linker.to_vec(), bases.next(30)].concat();
            let cdna = bases.next(CDNA_LENGTH);
            r1.push(&read);
            r2.push(&cdna);
            seen.extend([read, cdna]);
        }
        (detect(&r1, Some(&r2), &whitelists).unwrap(), seen)
    };
    for n in [3, 6, 999] {
        let (report, seen) = sample(n, &mut bases);
        assert!(
            matches!(report.verdict(), Verdict::NoMatch { .. }),
            "{n} reads: {:?}",
            report.verdict()
        );
        let json = report.to_json("R1.fq.gz".as_ref(), Some("R2.fq.gz".as_ref())).unwrap();
        let text = report.to_text(false);
        let link = text.lines().find(|l| l.contains("issues/new")).unwrap_or("");
        assert!(!link.is_empty(), "a no-match offers the kit-request link: {text}");
        for read in &seen {
            for piece in read.windows(12).filter(|w| !linker.windows(12).any(|l| l == *w)) {
                let piece = std::str::from_utf8(piece).unwrap();
                assert!(!json.contains(piece) && !link.contains(piece), "{n} reads: {piece}");
            }
        }
        assert!(json.contains("\"sketch\": null"), "{json}");
    }
    let (report, _) = sample(1_000, &mut bases);
    let json = report.to_json("R1.fq.gz".as_ref(), Some("R2.fq.gz".as_ref())).unwrap();
    assert!(json.contains("GTCACTGCTAGCATCG"), "the linker is shared from 1,000 reads: {json}");
}

/// Multiome is nuclei: `STARsolo` and kb counted exons only, while Cell Ranger counts introns too
/// (audit, 8 Oct 2026). Here the reads' barcodes are on the ARC list, not the v3 one.
#[test]
fn nuclei_commands_count_introns() {
    let mut bases = Bases(23);
    let arc: Vec<Vec<u8>> = (0..500).map(|_| bases.next(16)).collect();
    let dir = whitelist_dir("nuclei", &arc, &mut bases);
    let other: Vec<String> = (0..96).map(|_| String::from_utf8(bases.next(16)).unwrap()).collect();
    std::fs::rename(dir.join("10x_v3_permit.txt"), dir.join("10x_arc_v1_gex_permit.txt")).unwrap();
    std::fs::write(dir.join("10x_v3_permit.txt"), other.join("\n") + "\n").unwrap();
    let whitelists = Whitelists::load(&dir).unwrap();
    let (r1, r2) = tenx_v3_reads(&arc, &mut bases);
    let report = detect(&r1, Some(&r2), &whitelists).unwrap();
    assert_eq!(report.layout().map(|l| l.name), Some("10x Multiome GEX"));
    let words = |pipeline| {
        let emitted =
            report.emit(pipeline, "x_R1.fq.gz".as_ref(), Some("x_R2.fq.gz".as_ref()), true);
        let argv = emitted.unwrap().argv;
        argv.iter()
            .map(|a| match a {
                Arg::Text(t) => t.clone(),
                Arg::Var(p) => format!("${{{}}}", p.name),
            })
            .collect::<Vec<_>>()
    };
    let star = words(Pipeline::Starsolo);
    assert!(star.windows(2).any(|w| w == ["--soloFeatures", "GeneFull"]), "{star:?}");
    let kb = words(Pipeline::Kb);
    assert!(kb.windows(2).any(|w| w == ["--workflow", "nac"]), "{kb:?}");
    assert!(kb.contains(&"${KB_C2}".to_owned()) && kb.contains(&"--sum=total".to_owned()));
}

/// A stream reads like a file: the same reads from a pipe or download as from disk, plain or gzip.
#[test]
fn reads_a_stream_as_it_reads_a_file() {
    let mut bases = Bases(17);
    let mut text = String::new();
    for i in 0..300 {
        let seq = String::from_utf8(bases.next(28)).unwrap();
        std::fmt::Write::write_fmt(
            &mut text,
            format_args!("@r{i}\n{seq}\n+\n{}\n", "I".repeat(28)),
        )
        .unwrap();
    }
    let dir = Scratch::new("stream");
    let plain = dir.join("r.fq");
    std::fs::write(&plain, &text).unwrap();
    let gz = dir.join("r.fq.gz");
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    std::io::Write::write_all(&mut encoder, text.as_bytes()).unwrap();
    std::fs::write(&gz, encoder.finish().unwrap()).unwrap();
    for path in [&plain, &gz] {
        let from_file = Reads::from_file(path, 100).unwrap();
        let bytes = std::fs::read(path).unwrap();
        let from_stream = Reads::from_reader(bytes.as_slice(), path, 100).unwrap();
        assert_eq!(from_file.len(), 100, "stops at the limit");
        assert!(from_file.iter().eq(from_stream.iter()), "{}", path.display());
    }
    // A stream cut mid-record past the limit still gives the reads before the cut.
    let cut = &text.as_bytes()[..text.len() >> 1];
    assert_eq!(Reads::from_reader(cut, "cut".as_ref(), 100).unwrap().len(), 100);
}
