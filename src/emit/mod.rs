// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `--emit`: a ready command for a downstream pipeline, built from a report. Each pipeline gets
//! only what the chemistry decides (barcode and UMI positions, the barcode list, the strand, the
//! file order); what only the user knows (a genome index, a reference) is a named placeholder,
//! which the shell form fills from the environment and stops on when it is unset.
//!
//! A wrong command fails the way a wrong chemistry does, silently, so emitting refuses rather
//! than guesses: below HIGH (unless asked), on anything but gene expression, on long reads, on a
//! layout the pipeline cannot express, and wherever a pipeline would read the files wrongly
//! (Cell Ranger, which finds R1 and R2 by name, refuses swapped or archive-named files).
//!
//! Every value printed was checked against the pipeline itself (bench/emit): Cell Ranger 10.1.0's
//! argument parser, STAR 2.7.11b's parameter parser, kb 0.30.2 and simpleaf 0.30.0 end to end.

mod cellranger;
mod kb;
mod shell;
mod simpleaf;
mod starsolo;

pub use shell::render_shell;

use crate::input::fastq::Format;
use crate::layouts::{LAYOUTS, Layout as Structure};
use crate::{Confidence, LibraryType, ReadKind, Report, Strand, Verdict};
use std::fmt;
use std::path::Path;

/// The pipelines `--emit` writes for.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Pipeline {
    /// `cellranger count` (10x kits).
    CellRanger,
    /// `STARsolo` (`STAR --soloType CB_UMI_Simple`).
    Starsolo,
    /// `simpleaf quant` (alevin-fry).
    Simpleaf,
    /// `kb count` (kallisto | bustools).
    Kb,
}

impl Pipeline {
    /// Every pipeline, in the order `--help` lists them.
    pub const ALL: [Pipeline; 4] =
        [Pipeline::CellRanger, Pipeline::Starsolo, Pipeline::Simpleaf, Pipeline::Kb];

    /// The name `--emit` takes.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Pipeline::CellRanger => "cellranger",
            Pipeline::Starsolo => "starsolo",
            Pipeline::Simpleaf => "simpleaf",
            Pipeline::Kb => "kb",
        }
    }

    /// The pipeline `--emit` names, or None for a name it does not know.
    #[must_use]
    pub fn from_name(name: &str) -> Option<Pipeline> {
        let found = Pipeline::ALL.into_iter().find(|p| p.name() == name);
        assert!(found.is_none_or(|p| p.name() == name), "a found pipeline has the name asked");
        found
    }
}

/// Something the user must supply, filled from the environment variable of the same name.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct Placeholder {
    /// The environment variable, e.g. `STAR_INDEX`.
    pub name: &'static str,
    /// What it must hold, in a few words; printed when it is unset.
    pub meaning: &'static str,
}

/// One argument of an emitted command.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Arg {
    /// Literal text, quoted for the shell when printed.
    Text(String),
    /// A value the user supplies.
    Var(Placeholder),
}

/// A command for one pipeline: its arguments, what the user must fill in, and notes.
#[derive(Clone, PartialEq, Debug)]
#[non_exhaustive]
pub struct Emitted {
    /// The pipeline.
    pub pipeline: Pipeline,
    /// The program and its arguments, in order.
    pub argv: Vec<Arg>,
    /// What the user must supply, once each, in first-use order.
    pub placeholders: Vec<Placeholder>,
    /// What else the command assumes, as sentences.
    pub notes: Vec<String>,
    /// The matched layout's name, for the header.
    pub layout: &'static str,
    /// The verdict's confidence: `HIGH`, `MEDIUM` or `LOW`.
    pub confidence: &'static str,
    /// The layout's held-out record, as the text report words it.
    pub validation: String,
}

/// Why no command was emitted.
#[derive(Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum EmitRefusal {
    /// The verdict is not a match: nothing to emit for.
    NotAMatch,
    /// The match is below HIGH and `--emit-anyway` was not given.
    BelowHigh(&'static str),
    /// Long reads: these pipelines take short reads.
    LongReads,
    /// Not a gene-expression library; these commands quantify gene expression.
    NotGeneExpression(&'static str),
    /// The pipeline cannot express this layout (or tchem does not write it yet).
    Unsupported {
        /// The pipeline asked for.
        pipeline: &'static str,
        /// The matched layout.
        layout: &'static str,
        /// Why the pipeline cannot take it, as a clause.
        why: &'static str,
    },
    /// The cDNA strand is not known, and the pipeline needs it.
    StrandUnknown,
    /// One file was given; the pipeline needs the barcode and the cDNA file.
    OneFile,
    /// The files were swapped, and Cell Ranger takes R1 and R2 by name.
    SwappedForCellRanger,
    /// Cell Ranger needs Illumina-named FASTQs (`<sample>_S1_L001_R1_001.fastq.gz`).
    NotIlluminaNamed(String),
    /// A BAM, or one gzipped file beside a plain one: the pipeline cannot read them as given.
    InputFormat(String),
    /// The layout reads probes, not cDNA (Flex).
    ProbeBased(&'static str),
}

impl Emitted {
    /// The command as a shell script: a header, then one flag per line.
    #[must_use]
    pub fn to_shell(&self) -> String {
        let script = render_shell(self);
        assert!(script.ends_with('\n'), "a script ends its last line");
        assert!(script.starts_with("# tchem "), "a script opens with its header");
        script
    }

    /// The command as one JSON object (schema/emit.v1.json).
    #[must_use]
    pub fn to_json(&self) -> String {
        let json = crate::output::json::pretty(&crate::output::json::emit_value(self));
        assert!(json.starts_with('{') && json.ends_with('}'), "one JSON object");
        assert!(json.contains("\"emit_schema_version\""), "the object names its contract");
        json
    }
}

impl fmt::Display for EmitRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NotAMatch => write!(f, "no command: the verdict is not a match"),
            Self::BelowHigh(c) => write!(
                f,
                "no command: confidence is {c}, and a wrong chemistry fails silently \
                 (--emit-anyway emits it all the same)"
            ),
            Self::LongReads => {
                write!(f, "no command: long reads (use the verdict's wf-single-cell kit)")
            }
            Self::NotGeneExpression(kind) => write!(
                f,
                "no command: the library is {kind}, and these commands count gene expression"
            ),
            Self::Unsupported { pipeline, layout, why } => {
                write!(f, "no {pipeline} command for {layout}: {why}")
            }
            Self::StrandUnknown => write!(f, "no command: the cDNA strand could not be measured"),
            Self::OneFile => write!(f, "no command: give both the barcode and the cDNA file"),
            Self::SwappedForCellRanger => write!(
                f,
                "no command: R1 and R2 are swapped, and Cell Ranger reads them by name; \
                 rename the files so R1 holds the barcode read"
            ),
            Self::NotIlluminaNamed(file) => write!(
                f,
                "no command: Cell Ranger needs Illumina-named FASTQs \
                 (<sample>_S1_L001_R1_001.fastq.gz), and {file} is not; rename R1 and R2 to match"
            ),
            Self::InputFormat(why) => write!(f, "no command: {why}"),
            Self::ProbeBased(layout) => write!(
                f,
                "no command: {layout} reads probes, not cDNA; count it with cellranger multi and \
                 its probe set"
            ),
        }
    }
}

impl std::error::Error for EmitRefusal {}

/// What every pipeline's writer is given: the matched layout and the files in their roles.
pub(crate) struct Matched<'a> {
    pub structure: &'static Structure,
    pub layout: crate::Layout,
    /// The file holding the layout's barcode read, and the other one, after swap detection.
    pub barcode_file: &'a Path,
    pub cdna_file: Option<&'a Path>,
    pub strand: Option<Strand>,
    pub swapped: bool,
    /// The verdict's confidence.
    pub confidence: Confidence,
    /// How the barcode and cDNA files are stored (a BAM never gets here).
    pub stored: (Stored, Option<Stored>),
}

/// The command for `pipeline` from `report`; `r1` and `r2` are the files as given.
///
/// # Errors
/// An [`EmitRefusal`] saying why no command fits.
pub fn emit(
    report: &Report,
    pipeline: Pipeline,
    r1: &Path,
    r2: Option<&Path>,
    below_high: bool,
) -> Result<Emitted, EmitRefusal> {
    let matched = matched(report, r1, r2, below_high)?;
    let mut emitted = match pipeline {
        Pipeline::CellRanger => cellranger::emit(&matched)?,
        Pipeline::Starsolo => starsolo::emit(&matched)?,
        Pipeline::Simpleaf => simpleaf::emit(&matched)?,
        Pipeline::Kb => kb::emit(&matched)?,
    };
    if matched.confidence != Confidence::High {
        emitted.notes.push(format!(
            "emitted at {} confidence because --emit-anyway was given: check the verdict first",
            matched.confidence.label()
        ));
    }
    assert_eq!(emitted.pipeline, pipeline, "the writer asked for wrote the command");
    assert!(
        emitted.argv.iter().all(|a| match a {
            Arg::Var(p) => emitted.placeholders.contains(p),
            Arg::Text(t) => !t.is_empty(),
        }),
        "every placeholder used is declared, and no argument is empty"
    );
    Ok(emitted)
}

/// The checks every pipeline shares, and the files in their roles.
fn matched<'a>(
    report: &Report,
    r1: &'a Path,
    r2: Option<&'a Path>,
    below_high: bool,
) -> Result<Matched<'a>, EmitRefusal> {
    let Verdict::Match { layout, confidence } = report.verdict() else {
        return Err(EmitRefusal::NotAMatch);
    };

    if confidence != Confidence::High && !below_high {
        return Err(EmitRefusal::BelowHigh(confidence.label()));
    }
    if layout.read_kind == ReadKind::Long {
        return Err(EmitRefusal::LongReads);
    }
    match report.library() {
        Some(LibraryType::GeneExpression) | None => {}
        Some(other) => return Err(EmitRefusal::NotGeneExpression(other.label())),
    }
    let Some(structure) = LAYOUTS.iter().find(|l| l.name == layout.name) else {
        unreachable!("a matched layout {} is in the table", layout.name)
    };
    if crate::layouts::probe_based(structure) {
        return Err(EmitRefusal::ProbeBased(structure.name));
    }
    let (barcode_file, cdna_file) = match layout.barcode_input {
        crate::Input::First => (r1, r2),
        crate::Input::Second => (r2.unwrap_or(r1), Some(r1)),
    };
    printable(barcode_file)?;
    if let Some(cdna) = cdna_file {
        printable(cdna)?;
        distinct(barcode_file, cdna)?;
    }
    assert!(cdna_file.is_none_or(|c| c != barcode_file), "the barcode and cDNA files differ");
    let stored = (stored(barcode_file)?, cdna_file.map(stored).transpose()?);
    let strand = report.cdna_strand();
    let swapped = report.swapped();
    let matched =
        Matched { structure, layout, barcode_file, cdna_file, strand, swapped, confidence, stored };
    assert_eq!(matched.layout.name, matched.structure.name, "one layout, two views");
    Ok(matched)
}

/// A path the command can print as one argument: not empty (a library caller may pass one), and
/// no control character (a newline above all would break the script's one-argument-per-line
/// shape, and is never meant).
fn printable(path: &Path) -> Result<(), EmitRefusal> {
    let text = path.to_string_lossy();
    if text.is_empty() {
        return Err(EmitRefusal::InputFormat("an empty path names no file".into()));
    }
    if text.chars().any(char::is_control) {
        return Err(EmitRefusal::InputFormat(format!(
            "the path {} holds a control character; rename the file",
            text.escape_debug()
        )));
    }
    Ok(())
}

/// The barcode and cDNA files are two files, by name and, where both resolve, on disk.
fn distinct(barcode: &Path, cdna: &Path) -> Result<(), EmitRefusal> {
    let same_on_disk = matches!(
        (std::fs::canonicalize(barcode), std::fs::canonicalize(cdna)),
        (Ok(a), Ok(b)) if a == b
    );
    if barcode == cdna || same_on_disk {
        return Err(EmitRefusal::InputFormat(format!(
            "R1 and R2 are the same file ({}); give the barcode and the cDNA file",
            barcode.display()
        )));
    }
    Ok(())
}

/// The header fields every writer fills the same way.
pub(crate) fn emitted(pipeline: Pipeline, m: &Matched<'_>, argv: Vec<Arg>) -> Emitted {
    let mut placeholders: Vec<Placeholder> = Vec::new();
    for arg in &argv {
        if let Arg::Var(p) = arg
            && !placeholders.contains(p)
        {
            placeholders.push(*p);
        }
    }
    assert!(matches!(argv.first(), Some(Arg::Text(_))), "a command starts with its program");
    // Two placeholders may share a name with different meanings (kb's nac index): never both.
    let distinct = placeholders
        .iter()
        .enumerate()
        .all(|(i, p)| placeholders.iter().take(i).all(|q| q.name != p.name));
    assert!(distinct, "placeholder names are distinct");
    Emitted {
        pipeline,
        argv,
        placeholders,
        notes: Vec::new(),
        layout: m.layout.name,
        confidence: m.confidence.label(),
        validation: crate::output::text::validation_sentence(m.layout.validation),
    }
}

/// The text of a path argument. Paths the run opened are Unicode on every platform the
/// pipelines run on; one that is not cannot be passed through a shell faithfully. A relative
/// path starting with `-` gets `./` in front, or a pipeline would read it as an option (kb took
/// `-o_S1_L001_R1_001.fastq.gz` as `-o`).
pub(crate) fn path_text(path: &Path) -> Result<String, EmitRefusal> {
    let given = path
        .to_str()
        .ok_or_else(|| EmitRefusal::InputFormat(format!("{} is not UTF-8", path.display())))?;
    if given.is_empty() {
        return Err(EmitRefusal::InputFormat("an empty path names no file".into()));
    }
    let text = if given.starts_with('-') { format!("./{given}") } else { given.to_owned() };
    assert!(!text.starts_with('-'), "a path argument is never an option");
    assert!(
        text.ends_with(given) && text.len() - given.len() <= 2,
        "the path, ./ at most in front"
    );
    Ok(text)
}

/// A path for an option that splits its value at commas (several lanes or samples): one holding
/// a comma would be read as two paths.
pub(crate) fn unsplit(text: String, option: &'static str) -> Result<String, EmitRefusal> {
    if text.contains(',') {
        return Err(EmitRefusal::InputFormat(format!(
            "{option} splits its value at commas, and {text} holds one; rename it"
        )));
    }
    Ok(text)
}

/// The cDNA file, which every pipeline here needs.
pub(crate) fn cdna<'a>(m: &Matched<'a>) -> Result<&'a Path, EmitRefusal> {
    let cdna = m.cdna_file.ok_or(EmitRefusal::OneFile)?;
    assert!(cdna != m.barcode_file, "the cDNA and barcode files differ");
    assert!(!cdna.as_os_str().is_empty(), "a given file has a name");
    Ok(cdna)
}

/// How a FASTQ file is stored, read from its first bytes (as the reader reads it), not its name.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub(crate) enum Stored {
    Plain,
    Gzip,
}

/// How `file` is stored, from the reader's own reading of its first bytes; a BAM or an SRA
/// archive is refused, since these pipelines take FASTQ. A library caller may emit for files it
/// has not written yet (reads held in memory): when the file cannot be read, the name decides.
pub(crate) fn stored(file: &Path) -> Result<Stored, EmitRefusal> {
    let refusal = |what: &str| {
        EmitRefusal::InputFormat(format!("{} is {what}: convert it to FASTQ first", file.display()))
    };
    // A pipe or device gives its reads once, and the pipeline reads the file again; opening a
    // drained pipe here would also wait for a writer that never comes (seen: a hang).
    if let Ok(meta) = std::fs::metadata(file)
        && !meta.is_file()
    {
        return Err(EmitRefusal::InputFormat(format!(
            "{} is not a regular file (a pipe or a device); the pipeline reads it again, so \
             give it a file",
            file.display()
        )));
    }
    let ext = |e: &str| file.extension().is_some_and(|x| x == e);
    let read = crate::input::fastq::format_of(file);
    let gzip = match read {
        Ok(Format::Fastq { gzip }) => gzip,
        Ok(Format::Bam) => return Err(refusal("BAM")),
        Ok(Format::Sra) => return Err(refusal("an SRA archive")),
        Err(_) if ext("bam") => return Err(refusal("BAM")),
        Err(_) => ext("gz"),
    };
    Ok(if gzip { Stored::Gzip } else { Stored::Plain })
}

/// The barcode and UMI of a layout whose barcode read is one barcode then a UMI, from its
/// first base: (barcode length, list or None for an unlisted barcode, UMI length). None for
/// any other layout (several barcodes, linkers, a variable start).
pub(crate) fn single_barcode(s: &Structure) -> Option<(usize, bool, usize)> {
    use crate::layouts::Segment;
    if s.offsets != [0] {
        return None;
    }
    let parsed = match s.segments {
        [Segment::Barcode { len, .. }, Segment::Umi(umi), rest @ ..] => {
            Some((*len, true, *umi, rest))
        }
        [Segment::Unlisted { len, .. }, Segment::Umi(umi), rest @ ..] => {
            Some((*len, false, *umi, rest))
        }
        _ => None,
    };
    let (barcode, listed, umi, rest) = parsed?;
    if !rest.iter().all(|seg| matches!(seg, Segment::Tail(_))) {
        return None;
    }
    assert!(barcode >= 1 && umi >= 1, "{}: a barcode and a UMI occupy bases", s.name);
    assert!(rest.len() <= 1, "{}: at most one tail", s.name);
    Some((barcode, listed, umi))
}

/// The refusal for a layout a pipeline cannot take yet.
pub(crate) fn unsupported(pipeline: Pipeline, m: &Matched<'_>, why: &'static str) -> EmitRefusal {
    EmitRefusal::Unsupported { pipeline: pipeline.name(), layout: m.layout.name, why }
}

/// The strand a pipeline needs, or the refusal when it was not measured.
pub(crate) fn strand(m: &Matched<'_>) -> Result<Strand, EmitRefusal> {
    let strand = m.strand.ok_or(EmitRefusal::StrandUnknown)?;
    assert!(m.layout.read_kind == ReadKind::Short, "emitting is for short reads");
    Ok(strand)
}

#[cfg(test)]
mod tests {
    use super::Pipeline;

    /// `--emit` takes these names, and the JSON reports print them: lower case ASCII, distinct,
    /// and each one names its own pipeline back.
    #[test]
    fn names_are_lower_case() {
        for pipeline in Pipeline::ALL {
            let name = pipeline.name();
            assert!(!name.is_empty() && name.bytes().all(|b| b.is_ascii_lowercase()), "{name}");
            assert_eq!(Pipeline::from_name(name), Some(pipeline), "{name} names its pipeline");
        }
    }
}
