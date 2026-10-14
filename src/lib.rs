// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Tussl Chemistry: names a single-cell RNA library's layout (kit, chemistry, barcode list,
//! geometry, cDNA strand and library type) from its first reads.
//!
//! The `tchem` command is a thin wrapper over this library. Call it directly to configure a mapper
//! in the same process, from reads already in memory, with no subprocess and no temporary files:
//!
//! ```no_run
//! use tussl_chemistry::{Reads, Verdict, Whitelists, detect};
//!
//! # fn main() -> Result<(), Box<dyn std::error::Error>> {
//! // Load the barcode lists once (text, or packed by `tchem --pack`).
//! let whitelists = Whitelists::load("whitelists".as_ref())?;
//!
//! // A sample of reads: the head of a FASTQ file, or reads pushed in from memory.
//! let r1 = Reads::from_file("R1.fastq.gz".as_ref(), 100_000)?;
//! let mut r2 = Reads::with_capacity(100_000, 100_000 * 151);
//! r2.push(b"ACGTTGCAAGGCTTACCGATCGATTGCA");
//!
//! let report = detect(&r1, Some(&r2), &whitelists)?;
//! match report.verdict() {
//!     Verdict::Match { layout, confidence } => {
//!         println!("{} ({confidence:?}): geometry {}", layout.name, layout.geometry);
//!     }
//!     Verdict::Ambiguous { candidates, reason } => println!("{candidates:?}: {reason}"),
//!     Verdict::NoMatch { hint } => println!("no known layout: {hint:?}"),
//!     _ => {}
//! }
//! println!("{:?} {:?}", report.cdna_strand(), report.library());
//! # Ok(())
//! # }
//! ```
//!
//! Everything else (scoring, decision rules, k-mer tables) is private and may change between
//! releases without notice; the items exported here are the stable surface.
//!
//! # Panics
//!
//! Every function asserts its own invariants (TigerStyle): a failed assertion is a bug in this
//! library, never a reaction to input, and stops the program rather than return a wrong answer.
//! Bad input is a typed refusal (`ReadsRefusal`, `WhitelistRefusal`, `DetectRefusal`). The only
//! panics a caller can cause are breaking a documented limit, listed under each function's own
//! `# Panics` (`Reads::with_capacity`, `Reads::push`, `Reads::from_file`).

// Per-function `# Panics` sections would repeat the paragraph above on every function; the
// caller-caused panics are documented where they arise.
#![allow(clippy::missing_panics_doc)]

mod constants;
mod detect;
mod emit;
mod hash;
mod input;
mod layouts;
mod lists;
mod long_read;
mod output;
mod signals;

use constants::LONG_READS_EXAMINED;
use detect::score::{Lists, Score};
use input::whitelist::Whitelist;
use layouts::Mate;
use lists::ListId;
use long_read::detect::LongVerdict;
use std::fmt;
use std::path::{Path, PathBuf};

pub use constants::{
    FASTQ_LINE_BYTES_MAX, SAMPLE_BASES_MAX, SAMPLE_READS_DEFAULT, SAMPLE_READS_MAX,
};
pub use detect::decide::Confidence;
pub use emit::{Arg, EmitRefusal, Emitted, Pipeline, Placeholder, render_shell};
pub use input::fastq::{Reads, ReadsRefusal};
pub use input::whitelist::{WhitelistRefusal, packed_path};
pub use layouts::validation::{Level, Validation};
pub use output::json::{JsonRefusal, SCHEMA_VERSION};
pub use signals::instrument::{Evidence, Instrument, Platform};
pub use signals::library::LibraryType;

/// Every vendor barcode list the layouts refer to, loaded once and shared by every detection.
pub struct Whitelists {
    lists: Lists,
    dir: PathBuf,
}

impl Whitelists {
    /// Loads every list from `dir` (files as `scripts/fetch_whitelists.sh` names them), from the
    /// packed form where `--pack` wrote one, one thread per list. An optional list (one only a
    /// Cell Ranger install holds) that is absent is skipped: the layouts that need it are then
    /// reported as not checked.
    ///
    /// # Errors
    /// A missing required list, or a malformed, oversized, stale or corrupt list file, with its
    /// path and reason.
    pub fn load(dir: &Path) -> Result<Self, WhitelistRefusal> {
        let whitelists = Self::read(dir, Whitelist::load)?;
        assert!(
            whitelists.lists.loaded().all(|(id, _)| dir.join(id.file_name()).exists()),
            "every loaded list is on disk"
        );
        assert_eq!(whitelists.dir, dir, "the lists remember where they came from");
        Ok(whitelists)
    }

    /// As [`Whitelists::load`], from the text lists, and writes each list's packed form beside it
    /// (`<list>.packed`, atomically) so later loads take milliseconds.
    ///
    /// # Errors
    /// As [`Whitelists::load`], or a packed file that could not be written.
    pub fn pack(dir: &Path) -> Result<Self, WhitelistRefusal> {
        let whitelists = Self::read(dir, Whitelist::pack)?;
        assert!(whitelists.lists.loaded().count() >= 1, "at least one list is packed");
        let packed = |id: ListId| packed_path(&dir.join(id.file_name())).exists();
        assert!(
            whitelists.lists.loaded().all(|(id, _)| packed(id)),
            "every loaded list's packed form is on disk"
        );
        Ok(whitelists)
    }

    fn read(
        dir: &Path,
        read: fn(&Path) -> Result<Whitelist, WhitelistRefusal>,
    ) -> Result<Self, WhitelistRefusal> {
        let loaded: Vec<Option<(ListId, Whitelist)>> = std::thread::scope(|scope| {
            let handles: Vec<_> = ListId::ALL
                .iter()
                .map(|&id| {
                    let path = dir.join(id.file_name());
                    scope.spawn(move || {
                        if id.optional() && !path.exists() {
                            return Ok(None);
                        }
                        let list = read(&path)?;
                        let expected = id.barcode_length();
                        if list.barcode_length() != expected {
                            let len = list.barcode_length();
                            return Err(WhitelistRefusal::KitMismatch { path, len, expected });
                        }
                        Ok(Some((id, list)))
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|handle| {
                    // A loader that panicked broke an invariant: re-raise it here.
                    handle.join().unwrap_or_else(|payload| std::panic::resume_unwind(payload))
                })
                .collect::<Result<_, _>>()
        })?;
        assert_eq!(loaded.len(), ListId::ALL.len(), "one result per list id");
        assert!(
            loaded.iter().zip(ListId::ALL).all(|(slot, want)| slot.is_some() || want.optional()),
            "only an optional list may be absent"
        );
        // Lists::new asserts each list lands in its own id's slot.
        let lists = Lists::new(loaded.into_iter().flatten());
        Ok(Self { lists, dir: dir.to_path_buf() })
    }

    /// The directory the lists came from.
    #[must_use]
    pub fn dir(&self) -> &Path {
        assert!(self.lists.loaded().count() >= 1, "loaded lists come from a directory");
        assert!(!self.dir.as_os_str().is_empty(), "the directory has a name");
        &self.dir
    }

    /// Each loaded list's text file and how many distinct barcodes it holds, in a fixed order.
    #[must_use]
    pub fn lists(&self) -> Vec<(PathBuf, usize)> {
        let out: Vec<(PathBuf, usize)> = self
            .lists
            .loaded()
            .map(|(id, list)| (self.dir.join(id.file_name()), list.count()))
            .collect();
        assert!(out.len() <= ListId::ALL.len(), "at most one entry per list");
        assert!(out.iter().all(|(_, count)| *count >= 1), "a loaded list holds a barcode");
        out
    }

    /// Distinct barcodes across every loaded list.
    #[must_use]
    pub fn barcode_count(&self) -> usize {
        let loaded = self.lists.loaded().count();
        let total: usize = self.lists.loaded().map(|(_, list)| list.count()).sum();
        assert!(total >= loaded, "every list holds at least one barcode");
        assert!(
            total <= loaded * constants::BARCODES_MAX,
            "no list holds more than BARCODES_MAX barcodes"
        );
        total
    }
}

/// One of the two samples given to [`detect`], in the order given.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Input {
    /// The first sample (`r1`).
    First,
    /// The second sample (`r2`).
    Second,
}

/// Which kind of sequencing the sample came from, and so which detector read it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum ReadKind {
    /// Short reads (Illumina and similar): barcode and cDNA in separate reads, at fixed positions.
    Short,
    /// Long reads (Oxford Nanopore, PacBio): whole molecules, found by searching for the adapter.
    Long,
}

/// The strand a mapper should count the cDNA read on.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Strand {
    /// The cDNA read is the transcript's sense strand (3' kits).
    Forward,
    /// The cDNA read is antisense (5' kits).
    Reverse,
}

/// A matched layout, as a mapper needs it.
#[derive(Clone, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub struct Layout {
    /// The kit and chemistry, e.g. `10x 3' v3/v3.1`.
    pub name: &'static str,
    /// The read geometry in simpleaf / piscem / alevin-fry notation.
    pub geometry: &'static str,
    /// The barcode list files to correct against, under [`Whitelists::dir`]; empty for kits
    /// without a public list.
    pub barcode_lists: Vec<PathBuf>,
    /// Which given sample holds the barcode read, after swap detection.
    pub barcode_input: Input,
    /// Processing advice for this kit (merging barcode pairs, sublibraries), or empty.
    pub note: &'static str,
    /// The kind of reads the layout was found in.
    pub read_kind: ReadKind,
    /// The kit as Oxford Nanopore's `wf-single-cell` names it (`--kit`), for 10x kits; None for
    /// kits it does not support.
    pub long_read_kit: Option<&'static str>,
    /// How far this layout's HIGH calls have been checked on held-out runs, for this kind of
    /// reads. It describes the layout, not this run: it vouches for HIGH calls only.
    pub validation: Validation,
}

/// What the reads show.
#[derive(Clone, PartialEq, Debug)]
#[non_exhaustive]
pub enum Verdict {
    /// One layout fits.
    Match {
        /// The layout.
        layout: Layout,
        /// How firmly the evidence supports it.
        confidence: Confidence,
    },
    /// Layouts the reads cannot separate that would be processed differently.
    Ambiguous {
        /// The candidates' names.
        candidates: Vec<&'static str>,
        /// Why the reads cannot separate them.
        reason: String,
    },
    /// No known layout fits.
    NoMatch {
        /// What the files do show (a dropped barcode read, long reads, an unknown kit), if
        /// anything.
        hint: Option<String>,
    },
}

/// Why detection was refused.
#[derive(Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum DetectRefusal {
    /// A sample holds no reads.
    NoReads(Input),
}

impl fmt::Display for DetectRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NoReads(Input::First) => write!(f, "the first sample has no reads"),
            Self::NoReads(Input::Second) => write!(f, "the second sample has no reads"),
        }
    }
}

impl std::error::Error for DetectRefusal {}

/// Detects the layout of a library from a sample of its reads: `r1` and, when the library has
/// one, `r2`. Labels need not be right: swapped files are detected. A sample of long reads
/// (most longer than any short read) goes to the long-read detector, which finds the 10x kit by
/// searching each read for the adapter. When both files are long reads, it reads `r1` only.
///
/// # Errors
/// A sample with no reads.
pub fn detect(
    r1: &Reads,
    r2: Option<&Reads>,
    whitelists: &Whitelists,
) -> Result<Report, DetectRefusal> {
    if r1.is_empty() {
        return Err(DetectRefusal::NoReads(Input::First));
    }
    if r2.is_some_and(Reads::is_empty) {
        return Err(DetectRefusal::NoReads(Input::Second));
    }
    // ENA serves some long-read runs as two files. Each is a sample of the same library, so the
    // long-read detector reads the first and the caller sets the second aside.
    let long_pair = r2.is_some_and(long_read::is_long);
    if (r2.is_none() || long_pair) && long_read::is_long(r1) {
        let profile = detect::profile::Profile::of(r1);
        let long = long_read::detect(r1, profile, &whitelists.lists);
        assert_eq!(long.evidence.reads, r1.len().min(LONG_READS_EXAMINED), "reads examined");
        let instrument = signals::instrument::of(r1.origin());
        return Ok(Report {
            inner: Inner::Long(long),
            whitelists: whitelists.dir.clone(),
            instrument,
            archive_model: None,
        });
    }
    assert!(
        !(long_read::is_long(r1) && r2.is_some_and(long_read::is_long)),
        "two files of long reads go to the long-read detector"
    );
    let short = detect::detect(r1, r2, &whitelists.lists);
    assert_eq!(short.r1_profile.reads, r1.len(), "the report covers every sampled read");
    assert_eq!(short.r2_profile.is_some(), r2.is_some(), "the second sample is profiled if given");
    // The barcode read's name usually carries it; the cDNA read's when the first was stripped.
    let mut instrument = signals::instrument::of(r1.origin());
    if instrument.platform.is_none()
        && let Some(r2) = r2
    {
        instrument = signals::instrument::of(r2.origin());
    }
    Ok(Report {
        inner: Inner::Short(short),
        whitelists: whitelists.dir.clone(),
        instrument,
        archive_model: None,
    })
}

/// The outcome of [`detect`]: the verdict and the evidence behind it.
pub struct Report {
    inner: Inner,
    whitelists: PathBuf,
    /// The instrument the reads' source records (read names, BAM read group).
    instrument: Instrument,
    /// The instrument model an archive records for the run, when the reads came from one.
    archive_model: Option<String>,
}

enum Inner {
    Short(detect::Report),
    Long(long_read::LongReport),
}

impl Report {
    /// Which detector read the sample.
    #[must_use]
    pub fn read_kind(&self) -> ReadKind {
        let kind = match &self.inner {
            Inner::Short(_) => ReadKind::Short,
            Inner::Long(_) => ReadKind::Long,
        };
        if let Inner::Long(long) = &self.inner {
            assert!(long.profile.reads > 0, "a long report covers reads");
        }
        assert!(!self.whitelists.as_os_str().is_empty(), "a report knows its lists' directory");
        kind
    }

    /// The verdict, with the matched layout, the candidates, or a hint.
    #[must_use]
    pub fn verdict(&self) -> Verdict {
        let verdict = match &self.inner {
            Inner::Short(short) => self.verdict_short(short),
            Inner::Long(long) => self.verdict_long(long),
        };
        if let Verdict::Ambiguous { candidates, .. } = &verdict {
            assert!(candidates.len() >= 2, "an ambiguity names two layouts at least");
        }
        if let Verdict::Match { layout, .. } = &verdict {
            assert!(!layout.name.is_empty(), "a matched layout has a name");
        }
        verdict
    }

    fn verdict_short(&self, short: &detect::Report) -> Verdict {
        let verdict = match &short.verdict {
            detect::decide::Verdict::Match { best, confidence } => Verdict::Match {
                layout: self.layout_of(short, score_of(short, *best)),
                confidence: *confidence,
            },
            detect::decide::Verdict::Ambiguous { candidates, reason } => Verdict::Ambiguous {
                candidates: candidates.iter().map(|&i| score_of(short, i).layout.name).collect(),
                reason: reason.clone(),
            },
            detect::decide::Verdict::NoMatch => {
                Verdict::NoMatch { hint: output::json::hint_sentence(short) }
            }
        };
        assert!(!short.scores.is_empty(), "every layout is scored");
        assert!(short.r1_profile.reads > 0, "a report covers at least one read");
        verdict
    }

    fn verdict_long(&self, long: &long_read::LongReport) -> Verdict {
        let verdict = match &long.verdict {
            LongVerdict::Match { kit, confidence, .. } => {
                Verdict::Match { layout: self.layout_long(kit), confidence: *confidence }
            }
            LongVerdict::Ambiguous { kits, reason } => Verdict::Ambiguous {
                candidates: kits.iter().map(|k| k.layout.name).collect(),
                reason: reason.clone(),
            },
            LongVerdict::NoMatch { reason } => Verdict::NoMatch { hint: Some(reason.clone()) },
        };
        assert!(long.evidence.reads > 0, "a long report examined reads");
        assert!(long.evidence.adapter <= long.evidence.reads, "adapter reads are a subset");
        verdict
    }

    /// The matched layout, or None when the verdict is not a match.
    #[must_use]
    pub fn layout(&self) -> Option<Layout> {
        let layout = match self.verdict() {
            Verdict::Match { layout, .. } => Some(layout),
            _ => None,
        };
        if let Some(layout) = &layout {
            assert!(!layout.geometry.is_empty(), "a matched layout gives the mapper a geometry");
            assert_eq!(layout.read_kind, self.read_kind(), "the layout came from this detector");
        }
        layout
    }

    /// The strand to count the cDNA read on: the measured strand when clear, else the matched
    /// layout's (polyT means forward, TSO reverse). None unless the verdict is a match, and for
    /// long reads, which come in both orientations.
    #[must_use]
    pub fn cdna_strand(&self) -> Option<Strand> {
        let Inner::Short(short) = &self.inner else { return None };
        let matched = match &short.verdict {
            detect::decide::Verdict::Match { best, .. } => Some(score_of(short, *best)),
            _ => None,
        };
        let word = matched.and_then(|score| output::json::mapper_strand(short, Some(score)));
        let strand = word.map(|w| if w == "forward" { Strand::Forward } else { Strand::Reverse });
        if strand.is_some() {
            assert!(matched.is_some(), "only a match gives the mapper a strand");
        }
        assert_eq!(strand.is_some(), word.is_some(), "every strand word maps to a strand");
        strand
    }

    /// What kind of library the cDNA read shows, when a layout matched or tied and there is a cDNA
    /// read to look at. Short reads only, for now.
    #[must_use]
    pub fn library(&self) -> Option<LibraryType> {
        let Inner::Short(short) = &self.inner else { return None };
        let library = short.library.as_ref().map(signals::library::LibraryEvidence::call);
        if library.is_some() {
            assert!(short.r2_profile.is_some(), "library typing needs a cDNA read");
            assert!(
                !matches!(short.verdict, detect::decide::Verdict::NoMatch),
                "a library is typed only behind a layout"
            );
        }
        library
    }

    /// Layouts not checked because a list they need was not fetched (lists only a Cell Ranger
    /// install holds; `scripts/fetch_whitelists.sh --cellranger <dir>` fetches them).
    pub fn not_checked(&self) -> impl Iterator<Item = &'static str> {
        let unchecked = match &self.inner {
            Inner::Short(short) => short.unchecked,
            // Long reads are tested against lists every install has.
            Inner::Long(_) => layouts::LayoutSet::EMPTY,
        };
        if let Inner::Short(short) = &self.inner {
            assert!(
                short.scores.iter().all(|s| unchecked.iter().all(|l| l.name != s.layout.name)),
                "a layout not checked has no score"
            );
        }
        assert!(
            unchecked.iter().all(|l| l.lists().iter().any(ListId::optional)),
            "only a layout reading an optional list goes unchecked"
        );
        unchecked.iter().map(|l| l.name)
    }

    /// Whether the samples were given the wrong way round (barcode read as `r2`, or the other way).
    /// Never for long reads, which come as one sample.
    #[must_use]
    pub fn swapped(&self) -> bool {
        let Inner::Short(short) = &self.inner else { return false };
        if short.swapped {
            assert!(short.r2_profile.is_some(), "a swap needs two samples");
        }
        assert!(short.r1_profile.reads > 0, "a report covers at least one read");
        short.swapped
    }

    /// The report for people, exactly as `tchem` prints it; `all` lists every candidate.
    #[must_use]
    pub fn to_text(&self, all: bool) -> String {
        let body = match &self.inner {
            Inner::Short(short) => output::text::render(short, all),
            Inner::Long(long) => output::text::render_long(long),
        };
        let text =
            output::text::instrument_line(&self.instrument, self.archive_model.as_deref()) + &body;
        assert!(text.contains("VERDICT  "), "every report ends in a verdict");
        assert!(text.ends_with('\n'), "every report line is terminated");
        text
    }

    /// The instrument the reads' source records: read names, or a BAM read group's tags. Detection
    /// never reads it.
    #[must_use]
    pub fn instrument(&self) -> &Instrument {
        let instrument = &self.instrument;
        assert!(
            instrument.platform.is_some() == instrument.evidence.is_some(),
            "evidence iff known"
        );
        assert!(
            instrument.platform.is_some() || instrument.model.is_none(),
            "no model without one"
        );
        instrument
    }

    /// Records the instrument model an archive states for the run (ENA's `instrument_model`),
    /// reported beside what the reads say: the two disagree often enough to show both. The text
    /// is the caller's: it is kept as a read name is (bounded, control characters made visible),
    /// and a second call replaces the first.
    pub fn note_archive_model(&mut self, model: &str) {
        let kept = input::fastq::kept_name(model.trim().as_bytes());
        self.archive_model = Some(kept).filter(|m| !m.is_empty());
        // Library callers pass any text.
        assert!(self.archive_model.as_ref().is_none_or(|m| !m.contains('\n')), "one line");
        assert!(
            self.archive_model
                .as_ref()
                .is_none_or(|m| m.len() <= input::fastq::KEPT_NAME_BYTES_MAX),
            "a noted model is bounded"
        );
    }

    /// The JSON object `tchem --format json` prints. `r1` and `r2` are the file names to report
    /// for the two samples.
    ///
    /// # Errors
    /// A path that is not valid UTF-8, which JSON cannot carry faithfully.
    pub fn to_json(&self, r1: &Path, r2: Option<&Path>) -> Result<String, JsonRefusal> {
        let value = match &self.inner {
            Inner::Short(short) => output::json::build(short, r1, r2, &self.whitelists)?,
            Inner::Long(long) => output::json::build_long(long, r1, &self.whitelists)?,
        };
        let value =
            output::json::with_instrument(value, &self.instrument, self.archive_model.as_deref());
        let json = output::json::pretty(&value);
        assert!(json.starts_with('{'), "the output is one JSON object");
        assert!(json.ends_with('}'), "the output is one JSON object, closed");
        Ok(json)
    }

    /// A ready command for `pipeline` (`tchem --emit`): what the chemistry decides, with what only
    /// the user knows left as named placeholders. `r1` and `r2` are the files as given, before
    /// swap detection. Below HIGH it refuses unless `below_high`.
    ///
    /// # Errors
    /// An [`EmitRefusal`] saying why no command fits: not a HIGH match, not gene expression, a
    /// layout the pipeline cannot express, or files the pipeline would read wrongly.
    pub fn emit(
        &self,
        pipeline: Pipeline,
        r1: &Path,
        r2: Option<&Path>,
        below_high: bool,
    ) -> Result<Emitted, EmitRefusal> {
        let emitted = emit::emit(self, pipeline, r1, r2, below_high)?;
        assert_eq!(emitted.pipeline, pipeline, "the command is for the pipeline asked");
        assert!(!emitted.argv.is_empty(), "a command has a program");
        Ok(emitted)
    }

    fn layout_of(&self, short: &detect::Report, score: &Score) -> Layout {
        let barcode_lists: Vec<PathBuf> =
            score.layout.lists().iter().map(|id| self.whitelists.join(id.file_name())).collect();
        // The barcode read is the layout's mate, in the orientation detection chose.
        let first = matches!(score.layout.read, Mate::R1) != short.swapped;
        if !first {
            assert!(short.r2_profile.is_some(), "a barcode read in the second sample");
        }
        assert!(barcode_lists.len() <= score.layout.segments.len(), "one list per segment at most");
        Layout {
            name: score.layout.name,
            geometry: score.layout.geometry,
            barcode_lists,
            barcode_input: if first { Input::First } else { Input::Second },
            note: score.layout.note,
            read_kind: ReadKind::Short,
            long_read_kit: score.layout.long_read_kit,
            validation: layouts::validation::of(score.layout.name, ReadKind::Short),
        }
    }

    fn layout_long(&self, kit: &long_read::detect::Kit) -> Layout {
        let barcode_lists = vec![self.whitelists.join(kit.list.file_name())];
        assert!(kit.layout.long_read_kit.is_some(), "a long-read kit names its wf-single-cell kit");
        assert_eq!(
            kit.layout.lists(),
            lists::ListSet::EMPTY.with(kit.list),
            "a 10x kit has one barcode list"
        );
        Layout {
            name: kit.layout.name,
            geometry: kit.layout.geometry,
            barcode_lists,
            barcode_input: Input::First,
            note: kit.layout.note,
            read_kind: ReadKind::Long,
            long_read_kit: kit.layout.long_read_kit,
            validation: layouts::validation::of(kit.layout.name, ReadKind::Long),
        }
    }
}

/// The short-read score a verdict index names.
fn score_of(short: &detect::Report, index: usize) -> &Score {
    let Some(score) = short.scores.get(index) else {
        unreachable!("verdict index {index} past {} scores", short.scores.len())
    };
    assert!(index < short.scores.len(), "the index names a score");
    score
}
