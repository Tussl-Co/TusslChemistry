// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `tchem`: the command line over the `tussl_chemistry` library. It parses arguments, samples the
//! files, calls the library and prints what it returns; every judgement lives in the library.

mod cli;
mod remote;

use cli::{Command, Detect, Format, Inputs};
use std::error::Error;
use std::fmt;
use std::path::Path;
use std::path::PathBuf;
use std::process::ExitCode;
use std::time::Instant;
use tussl_chemistry::{EmitRefusal, ReadKind, Reads, SAMPLE_READS_MAX, Whitelists, packed_path};

/// Exit status for a refused command line, as distinct from a refused input (1).
const EXIT_USAGE: u8 = 2;
/// Exit status when `--emit` finds no command that fits (below HIGH, not gene expression, a
/// layout the pipeline cannot take, files it would read wrongly).
const EXIT_EMIT_REFUSED: u8 = 3;
/// Deepest chain of causes printed under an error.
const CAUSES_MAX: usize = 8;

/// A refusal that is only a sentence.
#[derive(Debug)]
struct Message(String);

impl fmt::Display for Message {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

impl Error for Message {}

fn main() -> ExitCode {
    let env_whitelists = std::env::var_os(cli::WHITELISTS_ENV);
    let command = match cli::parse(std::env::args_os().skip(1), env_whitelists) {
        Ok(command) => command,
        Err(refusal) => {
            eprintln!("error: {refusal}\n\nFor more information, try '--help'.");
            return ExitCode::from(EXIT_USAGE);
        }
    };
    // The parser bounds the sample; the library relies on it.
    if let Command::Detect(detect) = &command {
        assert!(detect.reads >= 1, "a sample has at least one read");
        assert!(detect.reads <= SAMPLE_READS_MAX, "a sample is capped at SAMPLE_READS_MAX");
    }
    let outcome = match command {
        Command::Help => {
            print!("{}", cli::USAGE);
            Ok(())
        }
        Command::Version => {
            println!("tchem {}", env!("CARGO_PKG_VERSION"));
            Ok(())
        }
        Command::Pack { whitelists } => pack(&whitelists),
        Command::Detect(detect) => run(&detect),
    };
    match outcome {
        Ok(()) => ExitCode::SUCCESS,
        Err(error) => {
            print_error(error.as_ref());
            // No command fits this run: distinct from a refused input, so a pipeline can tell.
            if error.downcast_ref::<EmitRefusal>().is_some() {
                ExitCode::from(EXIT_EMIT_REFUSED)
            } else {
                ExitCode::FAILURE
            }
        }
    }
}

/// `Error: <what>`, then each cause indented beneath it, numbered when there are several.
fn print_error(error: &dyn Error) {
    let message = error.to_string();
    assert!(!message.is_empty(), "a refusal says why");
    eprintln!("Error: {message}");
    let mut causes: Vec<&dyn Error> = Vec::with_capacity(CAUSES_MAX);
    let mut next = error.source();
    while let Some(cause) = next {
        if causes.len() == CAUSES_MAX {
            break;
        }
        causes.push(cause);
        next = cause.source();
    }
    assert!(causes.len() <= CAUSES_MAX, "causes are capped at {CAUSES_MAX}");
    match causes.as_slice() {
        [] => {}
        [only] => eprintln!("\nCaused by:\n    {only}"),
        many => {
            eprintln!("\nCaused by:");
            for (index, cause) in many.iter().enumerate() {
                eprintln!("    {index}: {cause}");
            }
        }
    }
}

fn pack(dir: &Path) -> Result<(), Box<dyn Error>> {
    let whitelists = Whitelists::pack(dir)?;
    let lists = whitelists.lists();
    assert!(!lists.is_empty(), "packing writes every list");
    for (list, count) in lists {
        let packed = packed_path(&list);
        assert!(packed.exists(), "{} was just written", packed.display());
        println!("{}: {count} barcodes", packed.display());
    }
    Ok(())
}

/// The first `reads` reads of the file at `path`; a file with none is refused here, by name, before
/// the library would refuse it by position.
fn run_sample(path: &Path, reads: usize) -> Result<Reads, Box<dyn Error>> {
    let sample = Reads::from_file(path, reads)?;
    if sample.is_empty() {
        return Err(Message(format!("{} has no reads", path.display())).into());
    }
    assert!(!sample.is_empty(), "a refused empty file never reaches detection");
    assert!(sample.len() <= reads, "the sample holds at most --reads reads");
    Ok(sample)
}

/// The samples detection reads: R1, maybe R2, any index files set aside, and where they came
/// from when that is not the command line (a run streamed by accession).
struct Sampled {
    r1: (PathBuf, Reads),
    r2: Option<(PathBuf, Reads)>,
    index: Vec<PathBuf>,
    source: Option<Source>,
}

/// Where a run's reads came from when it was read by accession, and what the archive records.
struct Source {
    line: String,
    instrument_model: Option<String>,
}

/// Samples the input files. Labelled files are taken as given; files by position, or a run's
/// files streamed by accession, are each sampled, index reads set aside, and the one or two left
/// are R1 and R2 in the order given.
fn run_inputs(inputs: &Inputs, reads: usize) -> Result<Sampled, Box<dyn Error>> {
    let (samples, source) = match inputs {
        Inputs::Labelled { r1, r2 } => {
            // R1 first, so a refusal names the first broken file in the order given.
            let first = (r1.clone(), run_sample(r1, reads)?);
            let r2 =
                r2.as_ref().map(|p| run_sample(p, reads).map(|s| (p.clone(), s))).transpose()?;
            return Ok(Sampled { r1: first, r2, index: Vec::new(), source: None });
        }
        Inputs::Files(files) => {
            let samples: Result<Vec<_>, _> =
                files.iter().map(|p| run_sample(p, reads).map(|s| (p.clone(), s))).collect();
            (samples?, None)
        }
        Inputs::Accession(accession) => {
            let remote::Run { run, urls, instrument_model } = remote::fastq_urls(accession)?;
            let streamed =
                urls.iter().map(|u| remote::stream(u, reads)).collect::<Result<Vec<_>, _>>()?;
            let bytes: u64 = streamed.iter().map(|(_, _, b)| b).sum();
            let samples: Vec<(PathBuf, Reads)> =
                streamed.into_iter().map(|(path, sample, _)| (path, sample)).collect();
            let names: Vec<&str> = samples.iter().map(|(p, _)| remote::short_name(p)).collect();
            // Megabytes to one decimal: the bytes are measured, not estimated.
            let shift = remote::MEGABYTE_SHIFT;
            let tenths = ((bytes & ((1 << shift) - 1)) * 10) >> shift;
            let megabytes = format!("{}.{tenths}", bytes >> shift);
            let source = format!(
                "Run {run} ({accession}), first reads streamed from ENA ({megabytes} MB): {}",
                names.join(", ")
            );
            (samples, Some(Source { line: source, instrument_model }))
        }
    };
    assert!(!samples.is_empty(), "positional files or a run's files: at least one");
    assert!(samples.len() <= cli::FILES_MAX, "at most a lane's files");
    set_aside_index(samples, source)
}

/// Files by position: index reads are set aside, and the one or two left are R1 and R2.
fn set_aside_index(
    samples: Vec<(PathBuf, Reads)>,
    source: Option<Source>,
) -> Result<Sampled, Box<dyn Error>> {
    let total = samples.len();
    assert!((1..=cli::FILES_MAX).contains(&total), "one to a lane's files, not {total}");
    let (index, kept): (Vec<_>, Vec<_>) =
        samples.into_iter().partition(|(_, sample)| sample.looks_like_index());
    let index: Vec<PathBuf> = index.into_iter().map(|(path, _)| path).collect();
    let mut kept = kept.into_iter();
    let (Some(first), second, None) = (kept.next(), kept.next(), kept.next()) else {
        let message = format!(
            "after setting aside {} index file(s), {} read file(s) are left; expected one or two \
             (R1 and R2)",
            index.len(),
            total - index.len()
        );
        return Err(Message(message).into());
    };
    assert!(index.len() + 1 + usize::from(second.is_some()) == total, "files accounted for");
    Ok(Sampled { r1: first, r2: second, index, source })
}

/// The line naming index files set aside, or nothing.
fn run_index_note(index: &[PathBuf]) -> String {
    let names: Vec<String> = index.iter().map(|p| p.display().to_string()).collect();
    assert_eq!(names.len(), index.len(), "one name per index file");
    assert!(index.len() <= cli::FILES_MAX, "at most a lane's files are set aside");
    if names.is_empty() {
        String::new()
    } else {
        format!("Index reads, set aside: {}\n", names.join(", "))
    }
}

/// The line that names a long-read run's second file, which the long-read detector does not read.
fn long_pair_note(kind: ReadKind, r2: Option<&Path>) -> String {
    let note = match (kind, r2) {
        (ReadKind::Long, Some(path)) => format!(
            "Second file, long reads too, set aside: {} (the files of a long-read run are \
             samples of one library)\n",
            path.display()
        ),
        _ => String::new(),
    };
    assert_eq!(note.is_empty(), kind == ReadKind::Short || r2.is_none(), "a note for a long pair");
    assert!(note.is_empty() || note.ends_with('\n'), "the note is one whole line");
    note
}

fn run(detect: &Detect) -> Result<(), Box<dyn Error>> {
    let start = Instant::now();
    let whitelists = Whitelists::load(&detect.whitelists)?;
    let loaded = start.elapsed();

    let sampled = run_inputs(&detect.inputs, detect.reads)?;
    let (r1_path, r1) = &sampled.r1;
    let r2_path = sampled.r2.as_ref().map(|(path, _)| path.as_path());
    let r2 = sampled.r2.as_ref().map(|(_, reads)| reads);
    let mut report = tussl_chemistry::detect(r1, r2, &whitelists)?;
    if let Some(model) = sampled.source.as_ref().and_then(|s| s.instrument_model.as_deref()) {
        report.note_archive_model(model);
    }
    if report.swapped() {
        assert!(r2.is_some(), "a swap needs two files");
    }

    let output = match (detect.emit, detect.format) {
        (Some(pipeline), format) => {
            let emitted = report.emit(pipeline, r1_path, r2_path, detect.emit_anyway)?;
            match format {
                Format::Text => emitted.to_shell(),
                Format::Json => emitted.to_json() + "\n",
            }
        }
        (None, Format::Text) => {
            let source = sampled.source.as_ref().map_or(String::new(), |s| format!("{}\n", s.line));
            let second = long_pair_note(report.read_kind(), r2_path);
            source + &run_index_note(&sampled.index) + &second + &report.to_text(detect.all)
        }
        (None, Format::Json) => report.to_json(r1_path, r2_path)? + "\n",
    };
    assert!(output.ends_with('\n'), "every output ends its last line");
    print!("{output}");
    eprintln!(
        "\n[{} reads; {} barcodes across {} lists loaded in {:.1} s; total {:.1} s]",
        r1.len(),
        whitelists.barcode_count(),
        whitelists.lists().len(),
        loaded.as_secs_f64(),
        start.elapsed().as_secs_f64()
    );
    // Where the reads came from, on stderr too, so a JSON run still says what it streamed.
    if let Some(source) = &sampled.source {
        eprintln!("[{}]", source.line);
    }
    Ok(())
}
