// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! A list of valid barcodes of one length (up to 16 bases), packed 2 bits per base into a u32.
//! A barcode counts as a hit when it is listed exactly or one substitution away, the same
//! tolerance Cell Ranger and alevin-fry use to correct sequencing errors.
//!
//! The list is held as a strictly increasing array of packed barcodes plus a bucket index on the
//! top bits, so a lookup reads one bucket of a few entries. Against a hash set this needs no
//! hashing at load (the 10x lists ship sorted) and a third of the memory, fixed at load.
//!
//! `--pack` writes that array beside each text list (`<list>.packed`), so later runs read 4 bytes
//! a barcode instead of parsing 17. A packed file is untrusted like the text: its header, size,
//! order and range are checked, it must match the text list's size when that is present, and any
//! fault is a refusal naming the fix (rerun `--pack`), never a silent fall back to the text.
use crate::constants::{BARCODES_MAX, WHITELIST_LINE_BYTES_MAX};
use std::fmt;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

pub const LEN_MAX: usize = 16;

/// The command that fetches every list: the repository's script, `tchem-fetch-whitelists` in the
/// conda package.
const FETCH_LISTS: &str = "scripts/fetch_whitelists.sh";
/// Most top bits the bucket index keys on: 2^20 buckets cost 4 MB and leave about 7 barcodes
/// per bucket on the largest 10x list.
const BUCKET_BITS_MAX: u32 = 20;
/// Packed file: magic, format version, barcode length, count, source text size, then the barcodes;
/// every field little-endian.
const PACKED_MAGIC: [u8; 4] = *b"TPWL";
const PACKED_VERSION: u32 = 1;
const PACKED_HEADER: usize = 24;
/// Where each header field sits: magic at 0, then version, barcode length, count (u32 each) and
/// the source text's size (u64).
const PACKED_OFFSET_VERSION: usize = 4;
const PACKED_OFFSET_LENGTH: usize = 8;
const PACKED_OFFSET_COUNT: usize = 12;
const PACKED_OFFSET_SOURCE: usize = 16;
const _: () = assert!(PACKED_OFFSET_VERSION == PACKED_MAGIC.len(), "the version follows the magic");
const _: () = assert!(PACKED_OFFSET_SOURCE + 8 == PACKED_HEADER, "the source size ends the header");
/// Read buffer for a text list: a mebibyte keeps syscalls rare on a 120 MB list.
const READ_BUFFER_BYTES: usize = 1 << 20;
/// Bytes per text line assumed when sizing the first allocation (a 16-base barcode and newline
/// is 17); only a guess, so it can be rough.
const LINE_BYTES_GUESS: u64 = 16;
/// Bytes per packed barcode: one little-endian u32.
const BARCODE_BYTES: usize = std::mem::size_of::<u32>();
/// Bytes read per step when streaming a packed list in.
const PACKED_CHUNK_BYTES: usize = 1 << 14;

/// Why a whitelist file was refused.
#[derive(Debug)]
#[non_exhaustive]
pub enum WhitelistRefusal {
    /// The file could not be opened.
    Open {
        /// The file.
        path: PathBuf,
        /// The operating system's error.
        source: io::Error,
    },
    /// The file could not be read.
    Read {
        /// The file.
        path: PathBuf,
        /// The operating system's error.
        source: io::Error,
    },
    /// Barcodes longer than the 16 bases a list can hold.
    TooLong {
        /// The file.
        path: PathBuf,
        /// The barcodes' length.
        len: usize,
    },
    /// A barcode of a different length from the ones before it.
    LengthMismatch {
        /// The file.
        path: PathBuf,
        /// The line, from 1.
        line: usize,
        /// This barcode's length.
        len: usize,
        /// The length of the barcodes before it.
        expected: usize,
    },
    /// A barcode with a base other than A, C, G or T.
    NonAcgt {
        /// The file.
        path: PathBuf,
        /// The line, from 1.
        line: usize,
    },
    /// A line longer than any barcode line.
    LineTooLong {
        /// The file.
        path: PathBuf,
        /// The line, from 1.
        line: usize,
    },
    /// More distinct barcodes than a list may hold.
    TooMany {
        /// The file.
        path: PathBuf,
    },
    /// No barcodes at all.
    Empty {
        /// The file.
        path: PathBuf,
    },
    /// A packed list that is corrupt, truncated, unsorted or stale.
    Packed {
        /// The file.
        path: PathBuf,
        /// What is wrong with it.
        reason: &'static str,
    },
    /// A list whose barcodes are not the length its kit's barcodes are.
    KitMismatch {
        /// The file.
        path: PathBuf,
        /// The list's barcode length.
        len: usize,
        /// The length the kit's barcodes are.
        expected: usize,
    },
    /// A packed list that could not be written.
    Write {
        /// The file.
        path: PathBuf,
        /// The operating system's error.
        source: io::Error,
    },
}

impl fmt::Display for WhitelistRefusal {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Open { path, .. } => {
                let message = format!(
                    "cannot open barcode list {}. Fetch the lists with {FETCH_LISTS} <dir>, or \
                     give --whitelists the directory that holds them",
                    path.display()
                );
                // A missing list is a missing fetch, almost always: the message names the fix.
                assert!(message.contains(FETCH_LISTS), "an open failure names {FETCH_LISTS}");
                f.write_str(&message)
            }
            Self::Read { path, source } => {
                write!(f, "read barcode list {}: {source}", path.display())
            }
            Self::TooLong { path, len } => write!(
                f,
                "{}: barcodes are {len} bases, at most {LEN_MAX} supported",
                path.display()
            ),
            Self::LengthMismatch { path, line, len, expected } => {
                write!(
                    f,
                    "{} line {line}: barcode is {len} bases, earlier ones {expected}",
                    path.display()
                )
            }
            Self::NonAcgt { path, line } => {
                write!(f, "{} line {line}: barcode has a non-ACGT base", path.display())
            }
            Self::LineTooLong { path, line } => {
                write!(f, "{} line {line}: over {WHITELIST_LINE_BYTES_MAX} bytes", path.display())
            }
            Self::TooMany { path } => {
                write!(f, "{}: more than {BARCODES_MAX} barcodes", path.display())
            }
            Self::Empty { path } => write!(f, "{}: no barcodes", path.display()),
            Self::Packed { path, reason } => {
                write!(f, "{}: {reason}; rerun with --pack to rebuild it", path.display())
            }
            Self::KitMismatch { path, len, expected } => write!(
                f,
                "{}: barcodes are {len} bases, but this kit's are {expected}",
                path.display()
            ),
            Self::Write { path, source } => write!(f, "write {}: {source}", path.display()),
        }
    }
}

impl std::error::Error for WhitelistRefusal {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Open { source, .. } | Self::Read { source, .. } | Self::Write { source, .. } => {
                Some(source)
            }
            _ => None,
        }
    }
}

pub struct Whitelist {
    len: usize,
    /// `len` as a u32 for bit arithmetic; at most `LEN_MAX` (16), so 2 bits per base fit a u32.
    len_bits: u32,
    /// Listed barcodes, packed, strictly increasing.
    sorted: Vec<u32>,
    /// `sorted[starts[b]..starts[b + 1]]` holds the barcodes whose top bits are `b`.
    starts: Vec<u32>,
    /// A barcode's bucket is `value >> bucket_shift`.
    bucket_shift: u32,
    /// For small lists only (BD, SPLiT-seq, Singleron), every barcode one substitution from a
    /// listed one and not itself listed, sorted: a one-off lookup is then one search, not up to
    /// 48. Empty for lists past `NEIGHBOURS_LIST_MAX`, which search variant by variant.
    neighbours: Vec<u32>,
}

/// Lists up to this many barcodes precompute their one-substitution neighbours: at most
/// 4096 x 16 x 3 = 196,608 values, built once at load.
const NEIGHBOURS_LIST_MAX: usize = 4096;

const _: () = assert!(2 * LEN_MAX <= u32::BITS as usize, "a barcode must pack into a u32");
const _: () = assert!(BARCODES_MAX < u32::MAX as usize, "bucket starts index the array as u32");

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Hit {
    Exact,
    OneOff,
    Miss,
}

fn code(base: u8) -> Option<u32> {
    match base {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' => Some(3),
        _ => None,
    }
}

/// Packs `sequence` (exactly the barcode). Returns the packed value and the position of a single
/// non-ACGT base (stored as A), or None when there are two or more.
fn pack(sequence: &[u8]) -> Option<(u32, Option<u32>)> {
    assert!(sequence.len() <= LEN_MAX, "barcode of {} bases past LEN_MAX", sequence.len());
    let mut value = 0u32;
    let mut unknown = None;
    for (i, &base) in (0u32..).zip(sequence) {
        let c = match code(base) {
            Some(c) => c,
            None if unknown.is_none() => {
                unknown = Some(i);
                0
            }
            None => return None,
        };
        value = (value << 2) | c;
    }
    // Hot: checked build only. The value fits the barcode's bits, and an N lies inside it.
    debug_assert!(
        sequence.len() == LEN_MAX || value >> (2 * sequence.len()) == 0,
        "packed past the bases"
    );
    if let Some(i) = unknown {
        debug_assert!((i as usize) < sequence.len(), "the N at {i} lies inside the barcode");
    }
    Some((value, unknown))
}

/// Where `--pack` (or [`crate::Whitelists::pack`]) puts the packed form of the text list at
/// `path`: beside it, as `<list>.packed`.
#[must_use]
pub fn packed_path(path: &Path) -> PathBuf {
    let mut name = path.as_os_str().to_owned();
    name.push(".packed");
    let packed = PathBuf::from(name);
    // `pack` writes here: it must never be the text list itself, nor leave its directory.
    assert_ne!(packed.as_path(), path, "the packed list is a separate file");
    assert_eq!(packed.parent(), path.parent(), "the packed list sits beside its text list");
    packed
}

impl Whitelist {
    /// Loads the list at `path`, from its packed form when `--pack` has written one.
    pub fn load(path: &Path) -> Result<Self, WhitelistRefusal> {
        // Callers name the text list; its packed form is found from it, never the other way.
        assert!(
            path.extension().is_none_or(|e| e != "packed"),
            "load takes the text list's path, not {}",
            path.display()
        );
        let packed = packed_path(path);
        let list =
            if packed.exists() { Self::load_packed(&packed, path) } else { Self::load_text(path) };
        if let Ok(list) = &list {
            assert!((1..=LEN_MAX).contains(&list.len), "a loaded list has 1 to 16 base barcodes");
        }
        list
    }

    /// Parses the text list at `path` and writes its packed form beside it, atomically: the bytes
    /// go to a temporary file that is renamed into place, so a reader never sees half a file.
    pub fn pack(path: &Path) -> Result<Self, WhitelistRefusal> {
        let list = Self::load_text(path)?;
        let source = std::fs::metadata(path)
            .map_err(|source| WhitelistRefusal::Read { path: path.to_path_buf(), source })?
            .len();
        let Ok(len) = u32::try_from(list.len) else {
            unreachable!("len {} is at most LEN_MAX", list.len)
        };
        let Ok(count) = u32::try_from(list.sorted.len()) else {
            unreachable!("count past BARCODES_MAX")
        };
        let mut bytes = Vec::with_capacity(PACKED_HEADER + BARCODE_BYTES * list.sorted.len());
        bytes.extend_from_slice(&PACKED_MAGIC);
        for field in [PACKED_VERSION, len, count] {
            bytes.extend_from_slice(&field.to_le_bytes());
        }
        bytes.extend_from_slice(&source.to_le_bytes());
        assert_eq!(bytes.len(), PACKED_HEADER, "header layout");
        for value in &list.sorted {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        let target = packed_path(path);
        let mut temporary = target.as_os_str().to_owned();
        temporary.push(".tmp");
        let temporary = PathBuf::from(temporary);
        let written = |source| WhitelistRefusal::Write { path: target.clone(), source };
        std::fs::write(&temporary, &bytes).map_err(written)?;
        std::fs::rename(&temporary, &target).map_err(written)?;
        Ok(list)
    }

    /// Reads a packed list, refusing anything `pack` would not have written for the text list at
    /// `text` (when that is present).
    fn load_packed(path: &Path, text: &Path) -> Result<Self, WhitelistRefusal> {
        let refuse = |reason| WhitelistRefusal::Packed { path: path.to_path_buf(), reason };
        let size = std::fs::metadata(path)
            .map_err(|source| WhitelistRefusal::Open { path: path.to_path_buf(), source })?
            .len();
        let most = u64::try_from(PACKED_HEADER + BARCODE_BYTES * BARCODES_MAX).unwrap_or(u64::MAX);
        if size > most {
            return Err(refuse("larger than any list this build accepts"));
        }
        let read_error = |source| WhitelistRefusal::Read { path: path.to_path_buf(), source };
        let mut file = File::open(path)
            .map_err(|source| WhitelistRefusal::Open { path: path.to_path_buf(), source })?;
        let mut header = [0u8; PACKED_HEADER];
        if size < u64::try_from(PACKED_HEADER).unwrap_or(u64::MAX) {
            return Err(refuse("shorter than its header"));
        }
        file.read_exact(&mut header).map_err(read_error)?;
        let word = |at: usize| {
            header
                .get(at..at + 4)
                .and_then(|b| <[u8; 4]>::try_from(b).ok())
                .map_or(0, u32::from_le_bytes)
        };
        if header.get(..4) != Some(&PACKED_MAGIC[..]) {
            return Err(refuse("not a packed barcode list"));
        }
        if word(PACKED_OFFSET_VERSION) != PACKED_VERSION {
            return Err(refuse("written by a different format version"));
        }
        let (len, count) =
            (word(PACKED_OFFSET_LENGTH) as usize, word(PACKED_OFFSET_COUNT) as usize);
        let source = header
            .get(PACKED_OFFSET_SOURCE..PACKED_HEADER)
            .and_then(|b| <[u8; 8]>::try_from(b).ok())
            .map_or(0, u64::from_le_bytes);
        if !(1..=LEN_MAX).contains(&len) || !(1..=BARCODES_MAX).contains(&count) {
            return Err(refuse("barcode length or count out of range"));
        }
        if usize::try_from(size).ok() != Some(PACKED_HEADER + BARCODE_BYTES * count) {
            return Err(refuse("size disagrees with its barcode count"));
        }
        if std::fs::metadata(text).is_ok_and(|m| m.len() != source) {
            return Err(refuse("stale: the text list changed after it was packed"));
        }
        // Streamed through one fixed buffer, so the file's bytes and the array never coexist.
        let mut sorted = Vec::with_capacity(count);
        let mut chunk = [0u8; PACKED_CHUNK_BYTES];
        let mut left = BARCODE_BYTES * count;
        while left > 0 {
            let Some(part) = chunk.get_mut(..left.min(PACKED_CHUNK_BYTES)) else {
                unreachable!("min keeps the part in the buffer")
            };
            file.read_exact(part).map_err(read_error)?;
            sorted.extend(
                part.chunks_exact(BARCODE_BYTES)
                    .filter_map(|b| <[u8; 4]>::try_from(b).ok())
                    .map(u32::from_le_bytes),
            );
            left -= part.len();
        }
        assert_eq!(sorted.len(), count, "4-byte reads over a 4-byte-aligned length");
        if !sorted.iter().zip(sorted.iter().skip(1)).all(|(a, b)| a < b) {
            return Err(refuse("barcodes not strictly increasing"));
        }
        if len < LEN_MAX && sorted.last().is_some_and(|&v| v >> (2 * len) != 0) {
            return Err(refuse("a barcode packs past the list's length"));
        }
        Ok(Self::from_sorted(len, sorted))
    }

    /// Loads a text list from an untrusted file: every line and the barcode count are capped, and
    /// any malformed line is a typed refusal.
    fn load_text(path: &Path) -> Result<Self, WhitelistRefusal> {
        let owned = || path.to_path_buf();
        let file =
            File::open(path).map_err(|source| WhitelistRefusal::Open { path: owned(), source })?;
        let size = file.metadata().map_or(0, |m| m.len());
        let mut reader = BufReader::with_capacity(READ_BUFFER_BYTES, file);
        // A 16-base list is 17 bytes a line; the guess (16 bytes a line) only sizes the first
        // allocation.
        let guess = usize::try_from(size.div_euclid(LINE_BYTES_GUESS))
            .unwrap_or(BARCODES_MAX)
            .min(BARCODES_MAX);
        let mut values: Vec<u32> = Vec::with_capacity(guess);
        let mut len = 0;
        let mut raw = Vec::with_capacity(WHITELIST_LINE_BYTES_MAX + 1);
        let cap = u64::try_from(WHITELIST_LINE_BYTES_MAX + 1).unwrap_or(u64::MAX);
        // Bounded: each iteration consumes at least one byte of a finite file, and the count of
        // kept barcodes is capped by BARCODES_MAX.
        for line in 1.. {
            raw.clear();
            let n = (&mut reader)
                .take(cap)
                .read_until(b'\n', &mut raw)
                .map_err(|source| WhitelistRefusal::Read { path: owned(), source })?;
            if n == 0 {
                break;
            }
            if n > WHITELIST_LINE_BYTES_MAX && raw.last() != Some(&b'\n') {
                return Err(WhitelistRefusal::LineTooLong { path: owned(), line });
            }
            if let Some(value) = Self::load_text_line(&raw, &mut len, path, line)? {
                values.push(value);
            }
            // Duplicates count once: past twice the cap, compact, and refuse only when the
            // distinct barcodes themselves pass it. Memory stays under 2 * BARCODES_MAX entries.
            if values.len() > 2 * BARCODES_MAX {
                sort_dedup(&mut values);
                if values.len() > BARCODES_MAX {
                    return Err(WhitelistRefusal::TooMany { path: owned() });
                }
            }
            assert!(values.len() <= 2 * BARCODES_MAX, "parsed barcodes past twice the cap");
        }
        if values.is_empty() {
            return Err(WhitelistRefusal::Empty { path: owned() });
        }
        sort_dedup(&mut values);
        if values.len() > BARCODES_MAX {
            return Err(WhitelistRefusal::TooMany { path: owned() });
        }
        assert!((1..=LEN_MAX).contains(&len), "accepted barcodes of length {len}");
        Ok(Self::from_sorted(len, values))
    }

    /// One line of a text list: the packed barcode, None for a blank line, or the refusal. The
    /// first barcode sets `len`, and every later one must match it.
    fn load_text_line(
        raw: &[u8],
        len: &mut usize,
        path: &Path,
        line: usize,
    ) -> Result<Option<u32>, WhitelistRefusal> {
        assert!(
            raw.len() <= WHITELIST_LINE_BYTES_MAX + 1,
            "the reader caps a line at WHITELIST_LINE_BYTES_MAX + 1"
        );
        assert!(*len <= LEN_MAX, "the barcode length is set from an accepted barcode");
        let barcode = raw.trim_ascii();
        if barcode.is_empty() {
            return Ok(None);
        }
        if *len == 0 {
            if barcode.len() > LEN_MAX {
                return Err(WhitelistRefusal::TooLong {
                    path: path.to_path_buf(),
                    len: barcode.len(),
                });
            }
            *len = barcode.len();
        }
        if barcode.len() != *len {
            return Err(WhitelistRefusal::LengthMismatch {
                path: path.to_path_buf(),
                line,
                len: barcode.len(),
                expected: *len,
            });
        }
        let Some((value, None)) = pack(barcode) else {
            return Err(WhitelistRefusal::NonAcgt { path: path.to_path_buf(), line });
        };
        assert!((1..=LEN_MAX).contains(len), "an accepted barcode has 1 to 16 bases");
        Ok(Some(value))
    }

    /// Builds the bucket index over `sorted`: distinct barcodes of `len` bases, increasing.
    fn from_sorted(len: usize, sorted: Vec<u32>) -> Self {
        assert!((1..=LEN_MAX).contains(&len), "barcode length {len} outside 1..={LEN_MAX}");
        assert!(
            (1..=BARCODES_MAX).contains(&sorted.len()),
            "{} barcodes outside 1..={BARCODES_MAX}",
            sorted.len()
        );
        assert!(
            sorted.iter().zip(sorted.iter().skip(1)).all(|(a, b)| a < b),
            "barcodes must be strictly increasing"
        );
        let Ok(len_bits) = u32::try_from(len) else { unreachable!("len {len} is at most LEN_MAX") };
        let value_bits = 2 * len_bits;
        assert!(
            sorted.last().is_some_and(|&v| value_bits == u32::BITS || v >> value_bits == 0),
            "a barcode packs past {len} bases"
        );
        // About one bucket per barcode, capped; a 97-barcode list needs 128 buckets, not 2^20.
        let bucket_bits =
            sorted.len().next_power_of_two().ilog2().min(BUCKET_BITS_MAX).min(value_bits);
        let bucket_shift = value_bits - bucket_bits;
        let buckets = 1usize << bucket_bits;
        let mut starts = Vec::with_capacity(buckets + 1);
        let mut at = 0usize;
        for bucket in 0..buckets {
            // Advance past every barcode of an earlier bucket; `at` only moves forward.
            while sorted.get(at).is_some_and(|&v| bucket_of(v, bucket_shift) < bucket) {
                at += 1;
            }
            let Ok(start) = u32::try_from(at) else { unreachable!("{at} past BARCODES_MAX") };
            starts.push(start);
        }
        let Ok(end) = u32::try_from(sorted.len()) else { unreachable!("count past BARCODES_MAX") };
        starts.push(end);
        assert_eq!(starts.len(), buckets + 1, "one start per bucket plus the end");
        assert!(
            starts.iter().zip(starts.iter().skip(1)).all(|(a, b)| a <= b),
            "bucket starts must not decrease"
        );
        let mut list = Self { len, len_bits, sorted, starts, bucket_shift, neighbours: Vec::new() };
        if list.sorted.len() <= NEIGHBOURS_LIST_MAX {
            list.neighbours = list.neighbours_of_listed();
        }
        list
    }

    fn contains(&self, value: u32) -> bool {
        // Hot: checked build only. Probes are packed barcodes of this list's length, and the
        // index relies on `from_sorted`'s order: buckets ascend, so each run is a valid range.
        debug_assert!(
            self.len == LEN_MAX || value >> (2 * self.len_bits) == 0,
            "probe packs past the list's length"
        );
        let bucket = bucket_of(value, self.bucket_shift);
        let (Some(&start), Some(&end)) = (self.starts.get(bucket), self.starts.get(bucket + 1))
        else {
            unreachable!("bucket {bucket} of a {}-base barcode past the index", self.len)
        };
        let Some(run) = self.sorted.get(start as usize..end as usize) else {
            unreachable!("bucket run {start}..{end} past {} barcodes", self.sorted.len())
        };
        debug_assert!(start <= end, "bucket starts ascend"); // Hot: checked build only.
        run.binary_search(&value).is_ok()
    }

    /// A list built in memory, for tests.
    #[cfg(test)]
    pub(crate) fn from_barcodes(barcodes: &[&[u8]]) -> Self {
        let mut values: Vec<u32> =
            barcodes.iter().filter_map(|b| pack(b)).map(|(v, _)| v).collect();
        sort_dedup(&mut values);
        Self::from_sorted(barcodes.first().map_or(1, |b| b.len()), values)
    }

    /// The length of the list's barcodes.
    pub fn barcode_length(&self) -> usize {
        assert!((1..=LEN_MAX).contains(&self.len), "a list's barcodes are 1 to 16 bases");
        assert_eq!(self.len_bits as usize, self.len, "len_bits mirrors len");
        self.len
    }

    pub fn count(&self) -> usize {
        let count = self.sorted.len();
        // Pairs with `from_sorted`: a list holds 1 to BARCODES_MAX distinct barcodes.
        assert!(count >= 1, "a list holds at least one barcode");
        assert!(count <= BARCODES_MAX, "a list holds at most BARCODES_MAX barcodes");
        count
    }

    /// `sequence` is exactly one barcode's bases.
    pub fn hit(&self, sequence: &[u8]) -> Hit {
        assert_eq!(sequence.len(), self.len, "barcode window must match the list's length");
        // Hot: checked build only. `len_bits` mirrors `len` for the shift arithmetic below.
        debug_assert_eq!(self.len_bits as usize, self.len, "len_bits mirrors len");
        let shift = |i: u32| 2 * (self.len_bits - 1 - i);
        let Some((value, unknown)) = pack(sequence) else {
            return Hit::Miss;
        };
        if let Some(i) = unknown {
            // The N could be any base; one listed candidate makes it a one-off hit.
            let cleared = value & !(3 << shift(i));
            let listed = (0..4).any(|c| self.contains(cleared | (c << shift(i))));
            return if listed { Hit::OneOff } else { Hit::Miss };
        }
        if self.contains(value) {
            return Hit::Exact;
        }
        if !self.neighbours.is_empty() {
            let one_off = self.neighbours.binary_search(&value).is_ok();
            return if one_off { Hit::OneOff } else { Hit::Miss };
        }
        self.hit_by_search(value)
    }

    /// Every value one substitution from a listed barcode and not itself listed, sorted.
    fn neighbours_of_listed(&self) -> Vec<u32> {
        let shift = |i: u32| 2 * (self.len_bits - 1 - i);
        let mut out = Vec::with_capacity(self.sorted.len() * self.len * 3);
        for &value in &self.sorted {
            for i in 0..self.len_bits {
                let s = shift(i);
                let (original, cleared) = ((value >> s) & 3, value & !(3 << s));
                out.extend((0..4).filter(|&c| c != original).map(|c| cleared | (c << s)));
            }
        }
        out.sort_unstable();
        out.dedup();
        out.retain(|&v| !self.contains(v));
        assert!(out.len() <= self.sorted.len() * self.len * 3, "at most 3 variants per base");
        assert!(out.iter().zip(out.iter().skip(1)).all(|(a, b)| a < b), "neighbours sorted");
        out
    }

    /// Whether some single substitution of `value` is listed, trying each in turn.
    fn hit_by_search(&self, value: u32) -> Hit {
        let shift = |i: u32| 2 * (self.len_bits - 1 - i);
        // Hot: checked build only.
        debug_assert!(!self.contains(value), "an exact hit is answered before the search");
        debug_assert!(self.neighbours.is_empty(), "a list with a neighbourhood answers from it");
        for i in 0..self.len_bits {
            let s = shift(i);
            let original = (value >> s) & 3;
            let cleared = value & !(3 << s);
            for c in 0..4 {
                if c != original && self.contains(cleared | (c << s)) {
                    return Hit::OneOff;
                }
            }
        }
        Hit::Miss
    }
}

/// A barcode's bucket. A one-bucket index shifts by the whole width, which `>>` would refuse.
fn bucket_of(value: u32, shift: u32) -> usize {
    // Hot: checked build only. A shift covers at most the whole 32-bit barcode.
    debug_assert!(shift <= u32::BITS, "bucket shift {shift} past a u32");
    let bucket = value.checked_shr(shift).unwrap_or(0) as usize;
    // A probe that fits its list lands in one of at most 2^BUCKET_BITS_MAX buckets.
    debug_assert!(bucket >> BUCKET_BITS_MAX == 0, "bucket {bucket} past the index");
    bucket
}

/// Sorts and drops repeats; the 10x lists ship sorted, so the sort is usually a linear check.
fn sort_dedup(values: &mut Vec<u32>) {
    let before = values.len();
    if !values.is_sorted() {
        values.sort_unstable();
    }
    values.dedup();
    // Pairs with `from_sorted`'s precondition. Load-time only, so the linear check is cheap.
    assert!(values.len() <= before, "deduplication never adds barcodes");
    let increasing = values.iter().zip(values.iter().skip(1)).all(|(a, b)| a < b);
    assert!(increasing, "deduplicated barcodes are strictly increasing");
}

#[cfg(test)]
mod tests {
    use super::*;

    fn list(barcodes: &[&str]) -> Whitelist {
        let bytes: Vec<&[u8]> = barcodes.iter().map(|b| b.as_bytes()).collect();
        Whitelist::from_barcodes(&bytes)
    }

    #[test]
    fn exact_one_off_n_and_miss() {
        let whitelist = list(&["AAACCTGAGAAACCAT"]);
        assert_eq!(whitelist.hit(b"AAACCTGAGAAACCAT"), Hit::Exact);
        assert_eq!(whitelist.hit(b"AAACCTGAGAAACCAG"), Hit::OneOff);
        assert_eq!(whitelist.hit(b"NAACCTGAGAAACCAT"), Hit::OneOff);
        assert_eq!(whitelist.hit(b"AAACCTGAGAAACCGG"), Hit::Miss);
        assert_eq!(whitelist.hit(b"NNACCTGAGAAACCAT"), Hit::Miss);
    }

    fn list_file(name: &str, body: &str) -> PathBuf {
        let dir = std::env::temp_dir().join(format!("tchem-whitelist-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join(name);
        std::fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn load_refuses_bad_lists() {
        assert!(matches!(
            Whitelist::load(&list_file("e", "\n\n")),
            Err(WhitelistRefusal::Empty { .. })
        ));
        assert!(matches!(
            Whitelist::load(&list_file("n", "ACGT\nACGN\n")),
            Err(WhitelistRefusal::NonAcgt { line: 2, .. })
        ));
        assert!(matches!(
            Whitelist::load(&list_file("m", "ACGT\nACG\n")),
            Err(WhitelistRefusal::LengthMismatch { line: 2, len: 3, expected: 4, .. })
        ));
        assert!(matches!(
            Whitelist::load(&list_file("l", &"A".repeat(17))),
            Err(WhitelistRefusal::LineTooLong { .. } | WhitelistRefusal::TooLong { .. })
        ));
        assert_eq!(Whitelist::load(&list_file("ok", " ACGT \r\nTTTT\n")).unwrap().count(), 2);
    }

    /// The bucketed array answers membership exactly as a hash set would, for every bucket width
    /// (tiny lists, short barcodes, full 16-base barcodes).
    #[test]
    fn membership_matches_a_hash_set() {
        let mut state = 7u64;
        let mut next = || {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            u32::try_from(state >> 32).unwrap()
        };
        for (len, n) in [(16usize, 50_000usize), (16, 1), (9, 97), (8, 3_000), (1, 3)] {
            let mask = if len == 16 { u32::MAX } else { (1u32 << (2 * len)) - 1 };
            let mut values: Vec<u32> = (0..n).map(|_| next() & mask).collect();
            let set: crate::hash::FxHashSet<u32> = values.iter().copied().collect();
            sort_dedup(&mut values);
            let whitelist = Whitelist::from_sorted(len, values);
            assert_eq!(whitelist.count(), set.len());
            assert!(set.iter().all(|&v| whitelist.contains(v)), "a listed barcode is missing");
            for _ in 0..100_000 {
                let probe = next() & mask;
                assert_eq!(
                    whitelist.contains(probe),
                    set.contains(&probe),
                    "len {len} probe {probe:#x}"
                );
            }
        }
    }

    /// A packed list answers like its text, and every way a packed file can be wrong is refused.
    #[test]
    fn packed_round_trip_and_refusals() {
        let text = list_file("p.txt", "TTTT\nACGT\nACGT\nGGCA\n");
        let from_text = Whitelist::load(&text).unwrap();
        Whitelist::pack(&text).unwrap();
        let packed = packed_path(&text);
        let from_packed = Whitelist::load(&text).unwrap();
        assert_eq!((from_packed.len, &from_packed.sorted), (from_text.len, &from_text.sorted));
        let good = std::fs::read(&packed).unwrap();
        let refused = |bytes: &[u8]| {
            std::fs::write(&packed, bytes).unwrap();
            match Whitelist::load(&text) {
                Err(WhitelistRefusal::Packed { reason, .. }) => reason,
                other => panic!("expected a packed refusal, got {:?}", other.map(|w| w.count())),
            }
        };
        assert_eq!(refused(&good[..10]), "shorter than its header");
        assert_eq!(refused(&good[..good.len() - 1]), "size disagrees with its barcode count");
        let mut bad = good.clone();
        bad[0] = b'X';
        assert_eq!(refused(&bad), "not a packed barcode list");
        let mut bad = good.clone();
        bad[4] = 9;
        assert_eq!(refused(&bad), "written by a different format version");
        let mut bad = good.clone();
        bad.swap(24, 28);
        bad.swap(25, 29);
        bad.swap(26, 30);
        bad.swap(27, 31);
        assert_eq!(refused(&bad), "barcodes not strictly increasing");
        let mut bad = good.clone();
        let last = bad.len() - 1;
        bad[last] = 0xff;
        assert_eq!(refused(&bad), "a barcode packs past the list's length");
        std::fs::write(&packed, &good).unwrap();
        std::fs::write(&text, "TTTT\nACGT\nGGCA\n").unwrap();
        assert_eq!(refused(&good), "stale: the text list changed after it was packed");
    }

    #[test]
    fn short_barcodes() {
        let whitelist = list(&["GTCGCTATA"]);
        assert_eq!(whitelist.hit(b"GTCGCTATA"), Hit::Exact);
        assert_eq!(whitelist.hit(b"GTCGCTATT"), Hit::OneOff);
        assert_eq!(whitelist.hit(b"GTCGCTAAT"), Hit::Miss);
    }

    /// The precomputed neighbourhood answers every query exactly as the variant-by-variant
    /// search does, including Ns and sequences far from the list.
    #[test]
    fn neighbours_answer_like_the_search() {
        let mut state = 41u64;
        let mut next = || {
            state = state.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            usize::try_from(state >> 33).unwrap()
        };
        for length in [8usize, 9, 16] {
            let barcodes: Vec<Vec<u8>> =
                (0..96).map(|_| (0..length).map(|_| b"ACGT"[next() % 4]).collect()).collect();
            let slices: Vec<&[u8]> = barcodes.iter().map(Vec::as_slice).collect();
            let list = Whitelist::from_barcodes(&slices);
            assert!(!list.neighbours.is_empty(), "a small list precomputes its neighbours");
            let mut by_search = Whitelist::from_barcodes(&slices);
            by_search.neighbours.clear();
            for round in 0..20_000 {
                let query: Vec<u8> = if round % 3 == 0 {
                    let mut q = barcodes[next() % barcodes.len()].clone();
                    let at = next() % length;
                    q[at] = b"ACGTN"[next() % 5];
                    q
                } else {
                    (0..length).map(|_| b"ACGTN"[next() % 5]).collect()
                };
                assert_eq!(
                    list.hit(&query),
                    by_search.hit(&query),
                    "{}",
                    String::from_utf8_lossy(&query)
                );
            }
        }
    }
}
