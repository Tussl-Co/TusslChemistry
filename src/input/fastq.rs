// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Reads the first sequences of a FASTQ. Layout detection needs bases only, so names and
//! qualities are checked for shape and dropped. Reading stops at `limit`, so a 30 GB file costs
//! the same as a 10 MB one.
//!
//! The file is untrusted input: every line is read through a cap, total bases are capped, and a
//! malformed record is a typed refusal, never an assertion.

use crate::constants::{
    FASTQ_LINE_BYTES_MAX, INDEX_COVER_MIN, INDEX_READ_BASES_MAX, INDEX_READS_EXAMINED,
    INDEX_TOP_SEQUENCES, SAMPLE_BASES_MAX, SAMPLE_READS_MAX,
};
use crate::detect::profile::ratio;
use crate::hash::{FxHashMap, map_with_capacity};
use flate2::read::MultiGzDecoder;
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

/// Read buffer over the (decompressed) file: a mebibyte keeps reads from the decoder large.
const READ_BUFFER_BYTES: usize = 1 << 20;
/// A gzip stream's first two bytes (RFC 1952).
const GZIP_MAGIC: [u8; 2] = [0x1f, 0x8b];
/// Starting room for one line: a short-read line fits without growing; longer lines grow once.
const LINE_CAPACITY_BYTES: usize = 512;

/// Sampled reads stored back to back: one buffer for every base and one for where each read ends,
/// instead of a heap allocation per read. Invariant: `ends` ascends and its last entry is
/// `bases.len()`.
#[derive(Debug, PartialEq, Eq, Default)]
pub struct Reads {
    bases: Vec<u8>,
    ends: Vec<usize>,
    /// What the source says about the instrument: the first record's name, and for BAM the read
    /// group's platform and model tags. Never read by detection.
    origin: Origin,
}

/// What a sample's source says about the instrument that read it.
#[derive(Debug, PartialEq, Eq, Default, Clone)]
pub struct Origin {
    /// The first record's name (FASTQ header without `@`, or BAM read name), as `kept_name`
    /// keeps it: at most `KEPT_NAME_BYTES_MAX` bytes.
    pub first_name: Option<String>,
    /// A BAM `@RG` line's `PL` (platform) and `PM` (model) tags.
    pub platform_tag: Option<String>,
    pub model_tag: Option<String>,
}

/// Bytes of a read name kept: instrument, run and flowcell come first in every format seen.
pub const NAME_BYTES_MAX: usize = 256;
/// Bytes a kept name may hold: each kept byte may become U+FFFD, three bytes in UTF-8.
pub const KEPT_NAME_BYTES_MAX: usize = NAME_BYTES_MAX * REPLACEMENT_BYTES;
const REPLACEMENT_BYTES: usize = 3;
const _: () = assert!(char::REPLACEMENT_CHARACTER.len_utf8() == REPLACEMENT_BYTES);

/// A name as kept: its first `NAME_BYTES_MAX` raw bytes as lossy UTF-8 (a character cut at the end
/// becomes U+FFFD), and every control character but tab replaced by U+FFFD, as invalid UTF-8 is.
/// Names are file bytes: a BAM name can hold a newline, a NUL or a terminal escape, and none of
/// them may reach a report.
pub fn kept_name(raw: &[u8]) -> String {
    let raw = raw.get(..raw.len().min(NAME_BYTES_MAX)).unwrap_or(raw);
    let visible =
        |c: char| if c.is_control() && c != '\t' { char::REPLACEMENT_CHARACTER } else { c };
    let name: String = String::from_utf8_lossy(raw).chars().map(visible).collect();
    assert!(name.len() <= KEPT_NAME_BYTES_MAX, "a kept name is bounded");
    assert!(name.chars().all(|c| !c.is_control() || c == '\t'), "a kept name is visible: {name:?}");
    name
}

impl Reads {
    /// The first `limit` reads of a FASTQ or BAM file, plain or gzip (any number of gzip members,
    /// so BGZF too), sequences only. The format is read from the file's first bytes, not its name.
    /// In BAM, secondary and supplementary alignments are skipped (they repeat a read). Fewer
    /// when the file is shorter, or when the sample reaches [`SAMPLE_BASES_MAX`] bases.
    ///
    /// # Errors
    /// An unreadable file, or a malformed record, naming the file and record.
    ///
    /// # Panics
    /// When `limit` is 0 or past [`SAMPLE_READS_MAX`].
    pub fn from_file(path: &Path, limit: usize) -> Result<Self, ReadsRefusal> {
        let reads = head_sequences(path, limit)?;
        assert!(reads.len() <= limit, "sampled {} reads past the limit {limit}", reads.len());
        assert!(reads.len() <= SAMPLE_READS_MAX, "a sample is capped at SAMPLE_READS_MAX");
        Ok(reads)
    }

    /// The first `limit` reads of a stream (a pipe, a download) holding FASTQ or BAM, plain or
    /// gzip, as [`Reads::from_file`] reads a file. Reading stops after `limit` records, so the
    /// stream need not end. `name` names the source in refusals.
    ///
    /// # Errors
    /// An unreadable stream, or a malformed record, naming the source and record.
    ///
    /// # Panics
    /// When `limit` is 0 or past [`SAMPLE_READS_MAX`].
    pub fn from_reader(source: impl Read, name: &Path, limit: usize) -> Result<Self, ReadsRefusal> {
        let reads = head_sequences_from(source, name, limit)?;
        assert!(reads.len() <= limit, "sampled {} reads past the limit {limit}", reads.len());
        assert!(reads.len() <= SAMPLE_READS_MAX, "a sample is capped at SAMPLE_READS_MAX");
        Ok(reads)
    }

    /// What the source says about the instrument (see [`Origin`]).
    #[must_use]
    pub fn origin(&self) -> &Origin {
        let origin = &self.origin;
        assert!(
            origin.first_name.as_ref().is_none_or(|n| n.len() <= KEPT_NAME_BYTES_MAX),
            "bounded"
        );
        let tags = [&origin.platform_tag, &origin.model_tag];
        assert!(
            tags.iter().all(|t| t.as_ref().is_none_or(|t| t.len() <= KEPT_NAME_BYTES_MAX)),
            "bounded"
        );
        origin
    }

    /// Records where the reads came from; the readers call it once.
    pub(crate) fn set_origin(&mut self, origin: Origin) {
        assert!(
            origin.first_name.as_ref().is_none_or(|n| !n.is_empty()),
            "a kept name says something"
        );
        self.origin = origin;
        assert!(
            self.origin.first_name.as_ref().is_none_or(|n| n.len() <= KEPT_NAME_BYTES_MAX),
            "bounded"
        );
    }

    /// An empty sample with room for `reads` reads holding `bases` bases in all, allocated once
    /// up front so filling it does not reallocate.
    ///
    /// # Panics
    /// When `reads` exceeds [`SAMPLE_READS_MAX`] or `bases` exceeds [`SAMPLE_BASES_MAX`].
    #[must_use]
    pub fn with_capacity(reads: usize, bases: usize) -> Self {
        assert!(reads <= SAMPLE_READS_MAX, "{reads} reads past the sample cap {SAMPLE_READS_MAX}");
        assert!(bases <= SAMPLE_BASES_MAX, "{bases} bases past the sample cap {SAMPLE_BASES_MAX}");
        Self {
            bases: Vec::with_capacity(bases),
            ends: Vec::with_capacity(reads),
            origin: Origin::default(),
        }
    }

    /// Appends one read's bases, upper-cased on the way in. Any byte is accepted; bases other
    /// than A, C, G and T count as unknown.
    ///
    /// # Panics
    /// When the sample already holds [`SAMPLE_READS_MAX`] reads, the read is longer than
    /// [`FASTQ_LINE_BYTES_MAX`], or the bases pass [`SAMPLE_BASES_MAX`] by more than one read:
    /// limits a caller checks with [`Reads::len`] before pushing.
    pub fn push(&mut self, sequence: &[u8]) {
        // A read is one FASTQ line, and the sample holds at most the cap of reads; the last read
        // may overshoot the base budget by one line, which `head_sequences` then stops on.
        assert!(
            sequence.len() <= FASTQ_LINE_BYTES_MAX,
            "a {}-base read is longer than a line",
            sequence.len()
        );
        assert!(
            self.ends.len() < SAMPLE_READS_MAX,
            "a read past the sample cap {SAMPLE_READS_MAX}"
        );
        let start = self.bases.len();
        self.bases.extend(sequence.iter().map(u8::to_ascii_uppercase));
        self.ends.push(self.bases.len());
        assert_eq!(self.bases.len(), start + sequence.len(), "a read is stored base for base");
        assert!(
            self.bases.len() < SAMPLE_BASES_MAX + FASTQ_LINE_BYTES_MAX,
            "bases past the sample budget"
        );
    }

    /// Bases held across every read.
    pub(crate) fn bases_len(&self) -> usize {
        assert!(
            self.ends.last().is_none_or(|&end| end == self.bases.len()),
            "the last read ends at the last base"
        );
        assert!(
            self.bases.len() < SAMPLE_BASES_MAX + FASTQ_LINE_BYTES_MAX,
            "bases within the budget"
        );
        self.bases.len()
    }

    /// How many reads the sample holds.
    #[must_use]
    pub fn len(&self) -> usize {
        let count = self.ends.len();
        assert!(count <= SAMPLE_READS_MAX, "{count} reads past the sample cap {SAMPLE_READS_MAX}");
        assert_eq!(
            self.ends.last().copied().unwrap_or(0),
            self.bases.len(),
            "the last read ends at the last base"
        );
        count
    }

    /// Whether the sample holds no reads.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        let empty = self.ends.is_empty();
        // Reads may be empty, so bases without reads is the broken case, not reads without bases.
        if empty {
            assert!(self.bases.is_empty(), "no reads, yet {} bases stored", self.bases.len());
        }
        assert!(
            self.ends.len() <= SAMPLE_READS_MAX,
            "reads past the sample cap {SAMPLE_READS_MAX}"
        );
        empty
    }

    /// Whether these look like sample-index reads (I1/I2): reads no longer than
    /// `INDEX_READ_BASES_MAX`, most of them one of a handful of sequences. Barcode and cDNA reads
    /// are longer and nearly all distinct.
    #[must_use]
    pub fn looks_like_index(&self) -> bool {
        let examined: Vec<&[u8]> = self.iter().take(INDEX_READS_EXAMINED).collect();
        if examined.is_empty() || examined.iter().any(|read| read.len() > INDEX_READ_BASES_MAX) {
            return false;
        }
        let mut counts: FxHashMap<&[u8], usize> = map_with_capacity(examined.len());
        for read in &examined {
            *counts.entry(read).or_insert(0) += 1;
        }
        let mut top: Vec<usize> = counts.into_values().collect();
        top.sort_unstable_by(|a, b| b.cmp(a));
        let covered: usize = top.iter().take(INDEX_TOP_SEQUENCES).sum();
        assert!(covered <= examined.len(), "the top sequences cover at most every read");
        assert!(examined.len() <= INDEX_READS_EXAMINED, "the test examines a bounded sample");
        ratio(covered, examined.len()) >= INDEX_COVER_MIN
    }

    /// Every read, in file order.
    pub fn iter(&self) -> impl Iterator<Item = &[u8]> + '_ {
        // The slicing below relies on `push`'s invariant: ends ascend and stop at bases.len().
        assert_eq!(
            self.ends.last().copied().unwrap_or(0),
            self.bases.len(),
            "the last read ends at the last base"
        );
        let starts = std::iter::once(0).chain(self.ends.iter().copied());
        starts.zip(&self.ends).map(|(start, &end)| {
            let Some(sequence) = self.bases.get(start..end) else {
                unreachable!("ends ascend and stop at bases.len(), so {start}..{end} is in range")
            };
            sequence
        })
    }

    /// The same reads in reverse order, for the debug build's order-independence check.
    #[cfg(debug_assertions)]
    pub(crate) fn reversed(&self) -> Self {
        let mut out = Self::with_capacity(self.len(), self.bases.len());
        let all: Vec<&[u8]> = self.iter().collect();
        for sequence in all.into_iter().rev() {
            out.push(sequence);
        }
        assert_eq!(out.ends.len(), self.ends.len(), "reversing keeps every read");
        assert_eq!(out.bases.len(), self.bases.len(), "reversing keeps every base");
        out
    }

    fn check(&self) {
        assert!(self.ends.windows(2).all(|w| w.first() <= w.last()), "read ends must ascend");
        assert_eq!(
            self.ends.last().copied().unwrap_or(0),
            self.bases.len(),
            "last end must be the base count"
        );
    }
}

/// Why a reads file (FASTQ or BAM) was refused. Messages name the file and record so the user
/// can find the problem.
#[derive(Debug)]
#[non_exhaustive]
pub enum ReadsRefusal {
    /// The file could not be opened.
    Open {
        /// The file.
        path: PathBuf,
        /// The operating system's error.
        source: io::Error,
    },
    /// The file could not be read or decompressed.
    Read {
        /// The operating system's error.
        source: io::Error,
    },
    /// A record that does not start with `@`.
    NoHeader {
        /// The file.
        path: PathBuf,
        /// The record (read), from 1.
        record: usize,
    },
    /// A record whose third line is not `+`.
    NoSeparator {
        /// The file.
        path: PathBuf,
        /// The record (read), from 1.
        record: usize,
    },
    /// A record whose quality line differs in length from its sequence.
    QualityLength {
        /// The file.
        path: PathBuf,
        /// The record (read), from 1.
        record: usize,
    },
    /// A line past `FASTQ_LINE_BYTES_MAX`.
    LineTooLong {
        /// The file.
        path: PathBuf,
        /// The record (read), from 1.
        record: usize,
    },
    /// An NCBI SRA archive, which needs converting to FASTQ first.
    SraArchive {
        /// The file.
        path: PathBuf,
    },
    /// A BAM file whose header breaks the format.
    BamHeader {
        /// The file.
        path: PathBuf,
        /// What is wrong.
        reason: &'static str,
    },
    /// A BAM record that breaks the format.
    BamRecord {
        /// The file.
        path: PathBuf,
        /// The record, from 1.
        record: usize,
        /// What is wrong.
        reason: &'static str,
    },
}

impl fmt::Display for ReadsRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Open { path, .. } => write!(f, "open {}", path.display()),
            Self::Read { source } => write!(f, "{source}"),
            Self::NoHeader { path, record } => {
                write!(f, "{}: record {record} does not start with '@'", path.display())
            }
            Self::NoSeparator { path, record } => {
                write!(f, "{}: record {record} has no '+' separator", path.display())
            }
            Self::QualityLength { path, record } => {
                write!(
                    f,
                    "{}: record {record} quality length differs from sequence",
                    path.display()
                )
            }
            Self::LineTooLong { path, record } => {
                write!(
                    f,
                    "{}: record {record} has a line over {FASTQ_LINE_BYTES_MAX} bytes",
                    path.display()
                )
            }
            Self::SraArchive { path } => write!(
                f,
                "{}: an NCBI SRA archive, not FASTQ or BAM; convert it with fasterq-dump (SRA \
                 Toolkit)",
                path.display()
            ),
            Self::BamHeader { path, reason } => write!(f, "{}: {reason}", path.display()),
            Self::BamRecord { path, record, reason } => {
                write!(f, "{}: record {record}: {reason}", path.display())
            }
        }
    }
}

impl std::error::Error for ReadsRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Open { source, .. } | Self::Read { source } => Some(source),
            _ => None,
        }
    }
}

/// Reads one line into `line` (cleared first), at most `FASTQ_LINE_BYTES_MAX` plus the newline.
/// Returns the bytes read; a line that hits the cap without ending is refused.
fn read_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    path: &Path,
    record: usize,
) -> Result<usize, ReadsRefusal> {
    assert!(record >= 1, "records are numbered from 1");
    line.clear();
    let cap = u64::try_from(FASTQ_LINE_BYTES_MAX + 1).unwrap_or(u64::MAX);
    let n =
        reader.take(cap).read_until(b'\n', line).map_err(|source| ReadsRefusal::Read { source })?;
    // The reader stops at the cap, and appends exactly what it reports into the cleared line.
    assert_eq!(n, line.len(), "read_until reports the bytes it appended");
    assert!(
        n <= FASTQ_LINE_BYTES_MAX + 1,
        "a {n}-byte line escaped the {FASTQ_LINE_BYTES_MAX}-byte cap"
    );
    if n > FASTQ_LINE_BYTES_MAX && line.last() != Some(&b'\n') {
        return Err(ReadsRefusal::LineTooLong { path: path.to_path_buf(), record });
    }
    Ok(n)
}

pub fn head_sequences(path: &Path, limit: usize) -> Result<Reads, ReadsRefusal> {
    assert!(limit >= 1, "a sample holds at least one read");
    assert!(limit <= SAMPLE_READS_MAX, "limit {limit} past the sample cap {SAMPLE_READS_MAX}");
    let file = File::open(path)
        .map_err(|source| ReadsRefusal::Open { path: path.to_path_buf(), source })?;
    let reads = head_sequences_from(file, path, limit)?;
    assert!(reads.len() <= limit, "sampled {} reads past the limit {limit}", reads.len());
    Ok(reads)
}

/// `head_sequences` over any stream (a file, a pipe, a download): the format is read from its
/// first bytes; `path` only names the source in refusals. Stops after `limit` records without
/// reading further, so a stream need not end.
pub fn head_sequences_from(
    source: impl Read,
    path: &Path,
    limit: usize,
) -> Result<Reads, ReadsRefusal> {
    assert!(limit >= 1, "a sample holds at least one read");
    assert!(limit <= SAMPLE_READS_MAX, "limit {limit} past the sample cap {SAMPLE_READS_MAX}");
    let (mut reader, format) = decoded(source)?;
    match format {
        Format::Sra => return Err(ReadsRefusal::SraArchive { path: path.to_path_buf() }),
        Format::Bam => return crate::input::bam::head_sequences(path, &mut reader, limit),
        Format::Fastq { .. } => {}
    }
    let mut reads = Reads::default();
    let mut line = Vec::with_capacity(LINE_CAPACITY_BYTES);
    let mut sequence = Vec::with_capacity(LINE_CAPACITY_BYTES);
    let mut first_name: Option<String> = None;

    // Bounded by `limit` records and by `SAMPLE_BASES_MAX`; each line by `FASTQ_LINE_BYTES_MAX`.
    while reads.len() < limit && reads.bases.len() < SAMPLE_BASES_MAX {
        let record = reads.len() + 1;
        if read_line(&mut reader, &mut line, path, record)? == 0 {
            break;
        }
        if line.first() != Some(&b'@') {
            return Err(ReadsRefusal::NoHeader { path: path.to_path_buf(), record });
        }
        if record == 1 {
            let name = trim_newline(line.get(1..).unwrap_or_default());
            first_name = Some(kept_name(name)).filter(|n| !n.is_empty());
        }
        read_line(&mut reader, &mut line, path, record)?;
        sequence.clear();
        sequence.extend_from_slice(trim_newline(&line));
        read_line(&mut reader, &mut line, path, record)?;
        if line.first() != Some(&b'+') {
            return Err(ReadsRefusal::NoSeparator { path: path.to_path_buf(), record });
        }
        read_line(&mut reader, &mut line, path, record)?;
        if trim_newline(&line).len() != sequence.len() {
            return Err(ReadsRefusal::QualityLength { path: path.to_path_buf(), record });
        }
        if reads.is_empty() {
            // Size both buffers once, from the first read, within the base budget.
            let expected = sequence.len().max(1).saturating_mul(limit).min(SAMPLE_BASES_MAX);
            reads = Reads::with_capacity(limit, expected);
        }
        reads.push(&sequence);
    }
    reads.set_origin(Origin { first_name, ..Origin::default() });
    assert!(reads.len() <= limit, "sampled {} reads past the limit {limit}", reads.len());
    assert!(
        reads.bases.len() < SAMPLE_BASES_MAX + FASTQ_LINE_BYTES_MAX,
        "sampled {} bases, past the {SAMPLE_BASES_MAX} budget by more than one line",
        reads.bases.len()
    );
    reads.check();
    Ok(reads)
}

/// What a source's first bytes show it to be. The reader and `--emit` both ask this, so a file
/// is never read as one format and emitted for as another.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Format {
    /// FASTQ, plain or gzip (any number of members, so BGZF too).
    Fastq { gzip: bool },
    /// BAM, plain or BGZF.
    Bam,
    /// An NCBI SRA archive, which holds reads but is no read format.
    Sra,
}

/// `source` decompressed when it is gzip, and its format, read from its first bytes.
fn decoded<'a>(
    source: impl Read + 'a,
) -> Result<(BufReader<Box<dyn Read + 'a>>, Format), ReadsRefusal> {
    let mut file = BufReader::with_capacity(READ_BUFFER_BYTES, source);
    let sra = opens_with(peek(&mut file)?, SRA_MAGIC);
    let gzip = opens_with(peek(&mut file)?, &GZIP_MAGIC);
    let inner: Box<dyn Read + 'a> = if gzip {
        // Multi-member: 10x and most sequencers write BGZF, many small gzip members; BAM is BGZF.
        Box::new(MultiGzDecoder::new(file))
    } else {
        Box::new(file)
    };
    let mut reader = BufReader::with_capacity(READ_BUFFER_BYTES, inner);
    let format = if sra {
        Format::Sra
    } else if opens_with(peek(&mut reader)?, &crate::input::bam::MAGIC) {
        Format::Bam
    } else {
        Format::Fastq { gzip }
    };
    Ok((reader, format))
}

/// The format of the file at `path`, read from its first bytes as the reader reads them.
///
/// # Errors
/// The file cannot be opened or its first bytes read.
pub fn format_of(path: &Path) -> Result<Format, ReadsRefusal> {
    let file = File::open(path)
        .map_err(|source| ReadsRefusal::Open { path: path.to_path_buf(), source })?;
    let (_, format) = decoded(file)?;
    Ok(format)
}

/// Bytes a format is recognised by: gzip's 2-byte magic, BAM's 4, an SRA archive's 8.
const SNIFF_BYTES: usize = 8;
/// An NCBI SRA archive (the `.sra` files SRA serves) opens with these bytes.
const SRA_MAGIC: &[u8] = b"NCBI.sra";
// The magics differ in their first byte, so at most one of them opens a file.
const _: () = assert!(
    matches!((SRA_MAGIC.first(), GZIP_MAGIC.first()), (Some(a), Some(b)) if *a != *b),
    "an SRA archive and gzip open with different bytes"
);

/// The next `SNIFF_BYTES` bytes `reader` holds (fewer at the stream's end), without consuming them.
fn peek(reader: &mut impl BufRead) -> Result<([u8; SNIFF_BYTES], usize), ReadsRefusal> {
    let bytes = reader.fill_buf().map_err(|source| ReadsRefusal::Read { source })?;
    let length = bytes.len().min(SNIFF_BYTES);
    let mut prefix = [0u8; SNIFF_BYTES];
    let (Some(target), Some(source)) = (prefix.get_mut(..length), bytes.get(..length)) else {
        unreachable!("length is at most both buffers' lengths")
    };
    target.copy_from_slice(source);
    assert!(length <= SNIFF_BYTES, "a peek copies at most SNIFF_BYTES");
    assert_eq!(prefix.get(..length), bytes.get(..length), "the prefix is the buffered start");
    Ok((prefix, length))
}

/// Whether a peeked prefix opens with `magic`.
fn opens_with((prefix, length): ([u8; SNIFF_BYTES], usize), magic: &[u8]) -> bool {
    assert!(magic.len() <= SNIFF_BYTES, "a magic fits the sniffed bytes");
    assert!(length <= SNIFF_BYTES, "a prefix holds at most SNIFF_BYTES");
    prefix.get(..length).is_some_and(|p| p.starts_with(magic))
}

fn trim_newline(line: &[u8]) -> &[u8] {
    let mut rest = line;
    while let [head @ .., b'\n' | b'\r'] = rest {
        rest = head;
    }
    assert!(line.starts_with(rest), "trimming keeps the line's start");
    assert!(!matches!(rest.last(), Some(b'\n' | b'\r')), "no line ending survives trimming");
    rest
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn file(name: &str, bytes: &[u8]) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tchem-fastq-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::File::create(&path).unwrap().write_all(bytes).unwrap();
        path
    }

    #[test]
    fn reads_records_and_strips_crlf() {
        let p = file("ok.fq", b"@a\r\nacgt\r\n+\r\nFFFF\r\n@b\nGG\n+\nFF\n");
        let reads = head_sequences(&p, 10).unwrap();
        assert_eq!(reads.iter().collect::<Vec<_>>(), vec![&b"ACGT"[..], &b"GG"[..]]);
        assert_eq!(head_sequences(&p, 1).unwrap().len(), 1);
    }

    #[test]
    fn refuses_a_line_past_the_cap_without_reading_it_all() {
        let p = file("long.fq", &vec![b'A'; FASTQ_LINE_BYTES_MAX + 10]);
        assert!(matches!(head_sequences(&p, 10), Err(ReadsRefusal::LineTooLong { record: 1, .. })));
    }

    #[test]
    fn refuses_malformed_records() {
        assert!(matches!(
            head_sequences(&file("h.fq", b">a\nAC\n"), 5),
            Err(ReadsRefusal::NoHeader { .. })
        ));
        assert!(matches!(
            head_sequences(&file("s.fq", b"@a\nAC\n-\nFF\n"), 5),
            Err(ReadsRefusal::NoSeparator { .. })
        ));
        assert!(matches!(
            head_sequences(&file("q.fq", b"@a\nACG\n+\nFF\n"), 5),
            Err(ReadsRefusal::QualityLength { .. })
        ));
    }

    #[test]
    fn index_reads_are_short_and_few() {
        let mut state = 3u64;
        let mut next = |n: u64| {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            usize::try_from((state >> 33) % n).unwrap()
        };
        let sample = |reads: &[Vec<u8>]| {
            let mut out = Reads::with_capacity(reads.len(), reads.len() * 30);
            for read in reads {
                out.push(read);
            }
            out
        };
        // A 10x single-index set: four 8-base oligos, with a sequencing error now and then.
        let oligos: [&[u8]; 4] = [b"AAACGGCG", b"CCTACCAT", b"GGCGTTTC", b"TTGTAAGA"];
        let index: Vec<Vec<u8>> = (0..2_000)
            .map(|i| {
                let mut read = oligos[i % 4].to_vec();
                if i % 50 == 0 {
                    read[next(8)] = b'N';
                }
                read
            })
            .collect();
        assert!(sample(&index).looks_like_index());
        let random = |length: usize, next: &mut dyn FnMut(u64) -> usize| -> Vec<u8> {
            (0..length).map(|_| b"ACGT"[next(4)]).collect()
        };
        let diverse: Vec<Vec<u8>> = (0..2_000).map(|_| random(8, &mut next)).collect();
        assert!(!sample(&diverse).looks_like_index(), "short but diverse reads are not an index");
        let barcodes: Vec<Vec<u8>> = (0..2_000).map(|_| random(28, &mut next)).collect();
        assert!(!sample(&barcodes).looks_like_index(), "barcode reads are not an index");
    }

    #[test]
    fn an_sra_archive_is_refused_by_name() {
        let path = std::env::temp_dir().join(format!("tchem-sra-{}.sra", std::process::id()));
        std::fs::write(&path, b"NCBI.sra\x88\x19\x03\x05 rest of the archive").unwrap();
        let refusal = head_sequences(&path, 10).unwrap_err();
        assert!(matches!(refusal, ReadsRefusal::SraArchive { .. }), "{refusal}");
        assert!(refusal.to_string().contains("fasterq-dump"), "{refusal}");
    }
}
