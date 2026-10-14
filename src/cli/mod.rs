// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The command line: a handful of long flags, parsed by hand. Every flag takes its value as the
//! next argument or after `=`; anything unknown, repeated, missing or out of range is refused
//! with the reason, before any work starts.

mod refusal;
mod usage;

pub use refusal::ArgumentRefusal;
pub use usage::USAGE;

use std::ffi::OsString;
use std::path::PathBuf;
use tussl_chemistry::{Pipeline, SAMPLE_READS_DEFAULT, SAMPLE_READS_MAX};

/// Most arguments read: every flag once with its value, and room for `--flag=value` spellings.
pub(crate) const ARGUMENTS_MAX: usize = 32;
/// Most files by position: one lane's R1, R2, I1 and I2.
pub const FILES_MAX: usize = 4;
/// The environment variable naming the whitelist directory when `--whitelists` is absent.
pub const WHITELISTS_ENV: &str = "TUSSL_CHEMISTRY_WHITELISTS";

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    Text,
    Json,
}

/// What to do: detect a library's layout, pack the whitelists, or print help or the version.
#[derive(Debug)]
pub enum Command {
    Detect(Detect),
    Pack { whitelists: PathBuf },
    Help,
    Version,
}

#[derive(Debug)]
pub struct Detect {
    pub inputs: Inputs,
    pub whitelists: PathBuf,
    pub reads: usize,
    pub all: bool,
    pub format: Format,
    /// `--emit`: print a command for this pipeline instead of the report.
    pub emit: Option<Pipeline>,
    /// `--emit-anyway`: emit below HIGH too.
    pub emit_anyway: bool,
}

/// The read files, as given.
#[derive(Debug)]
pub enum Inputs {
    /// `--r1` and maybe `--r2`: the files as labelled; swapped files are still detected.
    Labelled { r1: PathBuf, r2: Option<PathBuf> },
    /// A sample's files by position (1 to `FILES_MAX`): index reads are set aside, and what is
    /// left (one or two files) is detected in the order given.
    Files(Vec<PathBuf>),
    /// A public run or experiment by accession (SRR, ERR, DRR, SRX, ERX, DRX): its FASTQ files are
    /// streamed from ENA and taken as files by position.
    Accession(String),
}

/// Flags as they arrive, each at most once. The switches are independent flags, as typed.
#[allow(clippy::struct_excessive_bools)]
#[derive(Default)]
struct Seen {
    r1: Option<PathBuf>,
    r2: Option<PathBuf>,
    /// Files given by position, in order.
    files: Vec<PathBuf>,
    whitelists: Option<PathBuf>,
    reads: Option<usize>,
    all: bool,
    format: Option<Format>,
    emit: Option<Pipeline>,
    emit_anyway: bool,
    pack: bool,
    help: bool,
    version: bool,
    /// Flags recorded so far: each `set`/`set_flag` adds one, and `finish` checks the fields
    /// agree.
    given: usize,
}

/// Parses the arguments after the program name. `env_whitelists` is the environment's
/// `TUSSL_CHEMISTRY_WHITELISTS`, used when `--whitelists` is absent.
pub fn parse(
    arguments: impl IntoIterator<Item = OsString>,
    env_whitelists: Option<OsString>,
) -> Result<Command, ArgumentRefusal> {
    let mut seen = Seen::default();
    let mut arguments = arguments.into_iter();
    let mut count = 0usize;
    while let Some(argument) = arguments.next() {
        count += 1;
        if count > ARGUMENTS_MAX {
            return Err(ArgumentRefusal::TooMany);
        }
        let text = argument.to_string_lossy().into_owned();
        if !text.starts_with('-') {
            if seen.files.len() == FILES_MAX {
                return Err(ArgumentRefusal::TooManyFiles);
            }
            seen.files.push(PathBuf::from(argument));
            continue;
        }
        let (flag, inline) = match text.split_once('=') {
            Some((flag, value)) if flag.starts_with("--") => (flag.to_owned(), Some(value)),
            _ => (text.clone(), None),
        };
        let mut value = |name: &'static str| -> Result<OsString, ArgumentRefusal> {
            match inline {
                Some(value) => Ok(OsString::from(value)),
                None => arguments.next().ok_or(ArgumentRefusal::MissingValue(name)),
            }
        };
        parse_flag(&mut seen, &flag, &mut value)?;
        assert!(count <= ARGUMENTS_MAX, "the argument budget is checked first");
    }
    assert!(count <= ARGUMENTS_MAX, "loop exits within the argument budget");
    // A flag consumes one argument, or two with its value, so never more flags than arguments.
    assert!(seen.given + seen.files.len() <= count, "{} flags from {count} arguments", seen.given);
    assert!(seen.files.len() <= FILES_MAX, "positional files are capped");
    finish(seen, env_whitelists)
}

/// Records one flag into `seen`, taking its value through `value` when it has one.
fn parse_flag(
    seen: &mut Seen,
    flag: &str,
    value: &mut dyn FnMut(&'static str) -> Result<OsString, ArgumentRefusal>,
) -> Result<(), ArgumentRefusal> {
    // `parse` split any `--flag=value` at its first `=`.
    if flag.starts_with("--") {
        assert!(!flag.contains('='), "a long flag reaches here without its inline value");
    }
    let before = seen.given;
    let given = &mut seen.given;
    let recorded = match flag {
        "--r1" => set(given, &mut seen.r1, "--r1", PathBuf::from(value("--r1")?)),
        "--r2" => set(given, &mut seen.r2, "--r2", PathBuf::from(value("--r2")?)),
        "--whitelists" => {
            let directory = PathBuf::from(value("--whitelists")?);
            set(given, &mut seen.whitelists, "--whitelists", directory)
        }
        "--reads" => set(given, &mut seen.reads, "--reads", parse_reads(&value("--reads")?)?),
        "--format" => {
            let format = value("--format")?;
            let parsed = match format.to_str() {
                Some("text") => Format::Text,
                Some("json") => Format::Json,
                _ => {
                    return Err(ArgumentRefusal::BadValue {
                        flag: "--format",
                        value: format.to_string_lossy().into_owned(),
                        expected: "text or json",
                    });
                }
            };
            set(given, &mut seen.format, "--format", parsed)
        }
        "--emit" => set(given, &mut seen.emit, "--emit", parse_pipeline(&value("--emit")?)?),
        "--emit-anyway" => set_flag(given, &mut seen.emit_anyway, "--emit-anyway"),
        "--all" => set_flag(given, &mut seen.all, "--all"),
        "--pack" => set_flag(given, &mut seen.pack, "--pack"),
        "--help" | "-h" => set_flag(given, &mut seen.help, "--help"),
        "--version" | "-V" => set_flag(given, &mut seen.version, "--version"),
        other => Err(ArgumentRefusal::Unknown(other.to_owned())),
    };
    // A recorded flag counts once; a refused one leaves the count alone.
    let expected = if recorded.is_ok() { before + 1 } else { before };
    assert_eq!(seen.given, expected, "{flag}: flag count went from {before} to {}", seen.given);
    recorded
}

/// Stores a flag's value, refusing a second one.
fn set<T>(
    given: &mut usize,
    slot: &mut Option<T>,
    flag: &'static str,
    value: T,
) -> Result<(), ArgumentRefusal> {
    assert!(flag.starts_with("--"), "{flag} is recorded under its long name");
    if slot.is_some() {
        return Err(ArgumentRefusal::Repeated(flag));
    }
    *slot = Some(value);
    *given += 1;
    assert!(*given <= ARGUMENTS_MAX, "{given} flags recorded past the argument budget");
    Ok(())
}

/// Turns a switch on, refusing a second one.
fn set_flag(given: &mut usize, slot: &mut bool, flag: &'static str) -> Result<(), ArgumentRefusal> {
    assert!(flag.starts_with("--"), "{flag} is recorded under its long name");
    if *slot {
        return Err(ArgumentRefusal::Repeated(flag));
    }
    *slot = true;
    *given += 1;
    assert!(*given <= ARGUMENTS_MAX, "{given} flags recorded past the argument budget");
    Ok(())
}

fn parse_pipeline(value: &OsString) -> Result<Pipeline, ArgumentRefusal> {
    let text = value.to_string_lossy();
    let pipeline = Pipeline::from_name(&text).ok_or_else(|| ArgumentRefusal::BadValue {
        flag: "--emit",
        value: text.clone().into_owned(),
        expected: "cellranger, starsolo, simpleaf or kb",
    })?;
    assert_eq!(pipeline.name(), text, "the pipeline parsed is the one named");
    assert!(Pipeline::ALL.contains(&pipeline), "a known pipeline");
    Ok(pipeline)
}

fn parse_reads(value: &OsString) -> Result<usize, ArgumentRefusal> {
    let text = value.to_string_lossy();
    let refusal = || ArgumentRefusal::BadValue {
        flag: "--reads",
        value: text.clone().into_owned(),
        expected: "a whole number from 1 to 1000000",
    };
    let reads: usize = text.parse().map_err(|_| refusal())?;
    if !(1..=SAMPLE_READS_MAX).contains(&reads) {
        return Err(refusal());
    }
    assert!(reads >= 1, "a sample has at least one read");
    assert!(reads <= SAMPLE_READS_MAX, "a sample is capped at SAMPLE_READS_MAX reads");
    Ok(reads)
}

/// Checks the flags fit together and fills in defaults.
fn finish(seen: Seen, env_whitelists: Option<OsString>) -> Result<Command, ArgumentRefusal> {
    let recorded = [
        seen.r1.is_some(),
        seen.r2.is_some(),
        seen.whitelists.is_some(),
        seen.reads.is_some(),
        seen.all,
        seen.format.is_some(),
        seen.emit.is_some(),
        seen.emit_anyway,
        seen.pack,
        seen.help,
        seen.version,
    ];
    let recorded = recorded.iter().filter(|&&flag| flag).count();
    assert_eq!(recorded, seen.given, "every recorded flag was counted once");
    if seen.help {
        return Ok(Command::Help);
    }
    if seen.version {
        return Ok(Command::Version);
    }
    let env = env_whitelists.filter(|value| !value.is_empty()).map(PathBuf::from);
    let whitelists = seen.whitelists.or(env).ok_or(ArgumentRefusal::Missing("--whitelists"))?;
    if seen.pack {
        if seen.r1.is_some() || seen.r2.is_some() || !seen.files.is_empty() || seen.emit.is_some() {
            return Err(ArgumentRefusal::PackWithReads);
        }
        return Ok(Command::Pack { whitelists });
    }
    let accession = match seen.files.as_slice() {
        // A file on disk is read as a file, whatever its name.
        [only] if !only.exists() => {
            only.to_str().filter(|t| crate::remote::is_accession(t)).map(str::to_owned)
        }
        _ => None,
    };
    let inputs = match (seen.r1, seen.r2, seen.files.is_empty()) {
        (Some(r1), r2, true) => Inputs::Labelled { r1, r2 },
        (None, None, false) if accession.is_some() => {
            if seen.emit.is_some() {
                return Err(ArgumentRefusal::EmitWithAccession);
            }
            Inputs::Accession(accession.unwrap_or_default())
        }
        (None, None, false) => Inputs::Files(seen.files),
        (None, _, true) => return Err(ArgumentRefusal::Missing("--r1 (or files by position)")),
        (_, _, false) => return Err(ArgumentRefusal::MixedInputs),
    };
    let reads = seen.reads.unwrap_or(SAMPLE_READS_DEFAULT);
    assert!((1..=SAMPLE_READS_MAX).contains(&reads), "reads {reads} escaped its range check");
    let format = seen.format.unwrap_or(Format::Text);
    if seen.emit_anyway && seen.emit.is_none() {
        return Err(ArgumentRefusal::Missing("--emit (--emit-anyway qualifies it)"));
    }
    Ok(Command::Detect(Detect {
        inputs,
        whitelists,
        reads,
        all: seen.all,
        format,
        emit: seen.emit,
        emit_anyway: seen.emit_anyway,
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(arguments: &[&str], env: Option<&str>) -> Result<Command, ArgumentRefusal> {
        parse(arguments.iter().map(OsString::from), env.map(OsString::from))
    }

    #[test]
    fn detect_with_defaults_and_both_spellings() {
        let Ok(Command::Detect(d)) = run(&["--r1", "a.fq", "--r2=b.fq", "--whitelists", "w"], None)
        else {
            panic!("expected a detect command")
        };
        let Inputs::Labelled { r1, r2 } = d.inputs else { panic!("expected labelled files") };
        assert_eq!((r1, r2), (PathBuf::from("a.fq"), Some(PathBuf::from("b.fq"))));
        assert_eq!((d.reads, d.all, d.format), (SAMPLE_READS_DEFAULT, false, Format::Text));
        let Ok(Command::Detect(d)) =
            run(&["--r1", "a", "--reads", "5", "--all", "--format", "json"], Some("env"))
        else {
            panic!("expected a detect command")
        };
        assert_eq!(
            (d.whitelists, d.reads, d.all, d.format),
            (PathBuf::from("env"), 5, true, Format::Json)
        );
    }

    #[test]
    fn files_by_position() {
        let Ok(Command::Detect(d)) = run(&["--whitelists", "w", "I1.fq", "R1.fq", "R2.fq"], None)
        else {
            panic!("expected a detect command")
        };
        let Inputs::Files(files) = d.inputs else { panic!("expected positional files") };
        assert_eq!(files, ["I1.fq", "R1.fq", "R2.fq"].map(PathBuf::from));
        assert_eq!(run(&["--r1", "a", "b"], Some("w")).unwrap_err(), ArgumentRefusal::MixedInputs);
        assert_eq!(
            run(&["a", "b", "c", "d", "e"], Some("w")).unwrap_err(),
            ArgumentRefusal::TooManyFiles
        );
        assert_eq!(run(&["--pack", "a"], Some("w")).unwrap_err(), ArgumentRefusal::PackWithReads);
    }

    #[test]
    fn refusals() {
        assert_eq!(run(&["--r1"], Some("w")).unwrap_err(), ArgumentRefusal::MissingValue("--r1"));
        assert_eq!(
            run(&["--r1", "a"], None).unwrap_err(),
            ArgumentRefusal::Missing("--whitelists")
        );
        assert_eq!(
            run(&["--whitelists", "w"], None).unwrap_err(),
            ArgumentRefusal::Missing("--r1 (or files by position)")
        );
        assert_eq!(
            run(&["--r1", "a", "--r1", "b"], Some("w")).unwrap_err(),
            ArgumentRefusal::Repeated("--r1")
        );
        assert_eq!(
            run(&["--pack", "--r1", "a"], Some("w")).unwrap_err(),
            ArgumentRefusal::PackWithReads
        );
        assert_eq!(
            run(&["--bogus"], Some("w")).unwrap_err(),
            ArgumentRefusal::Unknown("--bogus".into())
        );
        for bad in ["0", "1000001", "-3", "x"] {
            assert!(matches!(
                run(&["--r1", "a", "--reads", bad], Some("w")),
                Err(ArgumentRefusal::BadValue { .. })
            ));
        }
        assert!(matches!(
            run(&["--r1", "a", "--format", "xml"], Some("w")),
            Err(ArgumentRefusal::BadValue { .. })
        ));
        let many: Vec<&str> = std::iter::repeat_n("--all", ARGUMENTS_MAX + 1).collect();
        assert!(matches!(
            run(&many, Some("w")),
            Err(ArgumentRefusal::Repeated(_) | ArgumentRefusal::TooMany)
        ));
    }

    #[test]
    fn emit_flags() {
        let Ok(Command::Detect(d)) = run(&["--r1", "a", "--emit", "starsolo"], Some("w")) else {
            panic!("expected a detect command")
        };
        assert_eq!((d.emit, d.emit_anyway), (Some(Pipeline::Starsolo), false));
        let Ok(Command::Detect(d)) = run(&["--r1", "a", "--emit=kb", "--emit-anyway"], Some("w"))
        else {
            panic!("expected a detect command")
        };
        assert_eq!((d.emit, d.emit_anyway), (Some(Pipeline::Kb), true));
        assert!(matches!(
            run(&["--r1", "a", "--emit", "salmon"], Some("w")),
            Err(ArgumentRefusal::BadValue { flag: "--emit", .. })
        ));
        assert!(matches!(
            run(&["--r1", "a", "--emit-anyway"], Some("w")),
            Err(ArgumentRefusal::Missing(_))
        ));
    }

    #[test]
    fn an_accession_alone_is_a_remote_run() {
        let Ok(Command::Detect(d)) = run(&["SRR19374511"], Some("w")) else {
            panic!("expected a detect command")
        };
        assert!(matches!(d.inputs, Inputs::Accession(ref a) if a == "SRR19374511"));
        // Beside another file it is a file name; with --emit it is refused (pipelines need files).
        assert!(matches!(run(&["SRR19374511", "b"], Some("w")), Ok(Command::Detect(_))));
        assert_eq!(
            run(&["--emit", "kb", "SRR19374511"], Some("w")).unwrap_err(),
            ArgumentRefusal::EmitWithAccession
        );
    }

    #[test]
    fn pack_help_and_version() {
        assert!(matches!(run(&["--pack", "--whitelists", "w"], None), Ok(Command::Pack { .. })));
        assert!(matches!(run(&["--help"], None), Ok(Command::Help)));
        assert!(matches!(run(&["--version", "--r1", "a"], None), Ok(Command::Version)));
    }
}
