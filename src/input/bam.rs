// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The sequences of a BAM file's first records. Long-read platforms deliver BAM, not FASTQ: PacBio
//! always (HiFi reads, segmented reads), Oxford Nanopore's basecaller by default. Only what the
//! detector needs is read: each record's bases, 4 bits each. Secondary and supplementary
//! alignments repeat a read already seen, so they are skipped.
//!
//! The format (SAM/BAM specification, section 4.2), little-endian: magic `BAM\1`; the header text
//! and reference names; then records, each a 4-byte length and the record. Every length is checked
//! against a cap before anything is read, and a file that breaks the format is a typed refusal.

use crate::constants::{FASTQ_LINE_BYTES_MAX, SAMPLE_BASES_MAX, SAMPLE_READS_MAX};
use crate::input::fastq::{Origin, Reads, ReadsRefusal, kept_name};
use std::io::{self, Read};
use std::path::Path;

/// The four bytes every BAM stream opens with.
pub const MAGIC: [u8; 4] = *b"BAM\x01";
/// Longest header text read: SAM headers hold program lines and read groups, rarely a megabyte.
const HEADER_TEXT_BYTES_MAX: usize = 1 << 26;
/// Most reference sequences named: a draft assembly can have millions of contigs.
const REFERENCES_MAX: usize = 1 << 24;
/// Longest reference name: SAM allows 255, with room for its terminator.
const REFERENCE_NAME_BYTES_MAX: usize = 1 << 10;
/// Longest record: a read's bases (half a byte each), its qualities (a byte each) and its tags,
/// which for PacBio subreads include per-base kinetics; 32 MiB covers a read of several megabases.
const RECORD_BYTES_MAX: usize = 1 << 25;
/// Fixed fields before the read name: reference, position, name length, quality, bin, operation
/// count, flags, sequence length, mate reference, mate position, template length.
const RECORD_FIXED_BYTES: usize = 32;
/// Flags of records that repeat a read: secondary (0x100) and supplementary (0x800) alignments.
const FLAG_REPEATS_READ: u16 = 0x100 | 0x800;
/// Bases by their 4-bit code (section 4.2.3).
const BASES: &[u8; 16] = b"=ACMGRSVTWYHKDBN";

const _: () = assert!(RECORD_FIXED_BYTES < RECORD_BYTES_MAX);
const _: () = assert!(FASTQ_LINE_BYTES_MAX + RECORD_FIXED_BYTES < RECORD_BYTES_MAX);

/// The bases of the first `limit` reads in a decompressed BAM stream that opens with `MAGIC`.
pub fn head_sequences(
    path: &Path,
    stream: &mut impl Read,
    limit: usize,
) -> Result<Reads, ReadsRefusal> {
    assert!(limit >= 1, "a sample holds at least one read");
    assert!(limit <= SAMPLE_READS_MAX, "limit {limit} past the sample cap {SAMPLE_READS_MAX}");
    let (platform_tag, model_tag) = head_sequences_header(path, stream)?;
    let mut first_name = None;
    let mut reads = Reads::default();
    let mut record = Vec::with_capacity(1 << 12);
    let mut sequence = Vec::with_capacity(1 << 12);
    let mut index = 0usize;
    // Bounded by `limit` reads, `SAMPLE_BASES_MAX` bases, and the file's end.
    while reads.len() < limit && reads.bases_len() < SAMPLE_BASES_MAX {
        index += 1;
        if !head_sequences_record(path, stream, index, &mut record)? {
            break;
        }
        // The first record's name only: a name built per record would allocate on every read.
        if index == 1 {
            first_name = record_name(&record);
        }
        if head_sequences_bases(path, index, &record, &mut sequence)? {
            if reads.is_empty() {
                let expected = sequence.len().max(1).saturating_mul(limit).min(SAMPLE_BASES_MAX);
                reads = Reads::with_capacity(limit, expected);
            }
            reads.push(&sequence);
        }
    }
    reads.set_origin(Origin { first_name, platform_tag, model_tag });
    assert!(reads.len() <= limit, "sampled {} reads past the limit {limit}", reads.len());
    assert!(reads.len() <= index, "one read per record at most");
    Ok(reads)
}

/// A record's read name (after the fixed fields, `l_read_name` bytes ending in a NUL), as kept.
fn record_name(record: &[u8]) -> Option<String> {
    let length = usize::from(*record.get(READ_NAME_LENGTH_AT)?);
    let raw = record.get(RECORD_FIXED_BYTES..RECORD_FIXED_BYTES + length)?;
    let raw = raw.strip_suffix(&[0]).unwrap_or(raw);
    Some(kept_name(raw)).filter(|n| !n.is_empty() && n != "*")
}

/// Offset of `l_read_name`, the read name's length, in a BAM record after its size (SAM spec).
const READ_NAME_LENGTH_AT: usize = 8;
const _: () = assert!(READ_NAME_LENGTH_AT < RECORD_FIXED_BYTES, "the length is a fixed field");

/// Header text scanned for `@RG` tags: read groups sit near the top of unaligned BAMs, and the
/// rest (reference lines, often megabytes) is skipped unread.
const READ_GROUP_SCAN_BYTES: usize = 1 << 20;
const _: () = assert!(READ_GROUP_SCAN_BYTES <= RECORD_BYTES_MAX, "read_exact takes the scan");

/// The first `@RG` line's `PL` (platform) and `PM` (model) tags, from the start of the text.
fn read_group_tags(text: &[u8]) -> (Option<String>, Option<String>) {
    assert!(text.len() <= READ_GROUP_SCAN_BYTES, "only the scanned start is parsed");
    let text = String::from_utf8_lossy(text);
    let line = text.lines().find(|l| l.starts_with("@RG\t"));
    let tag = |key: &str| {
        line.and_then(|l| l.split('\t').find_map(|f| f.strip_prefix(key)))
            .map(|v| kept_name(v.as_bytes()))
            .filter(|v| !v.is_empty())
    };
    let tags = (tag("PL:"), tag("PM:"));
    assert!(line.is_some() || tags == (None, None), "tags come from an @RG line");
    tags
}

fn refuse(path: &Path, reason: &'static str) -> ReadsRefusal {
    assert!(!reason.is_empty(), "a refusal says why");
    assert!(!path.as_os_str().is_empty(), "a refusal names its file");
    ReadsRefusal::BamHeader { path: path.to_path_buf(), reason }
}

fn read_exact(stream: &mut impl Read, buffer: &mut [u8]) -> io::Result<()> {
    let length = buffer.len();
    stream.read_exact(buffer)?;
    assert_eq!(buffer.len(), length, "read_exact fills the buffer");
    assert!(length <= RECORD_BYTES_MAX, "no read is larger than a record");
    Ok(())
}

/// A little-endian i32 length field, refused when negative or past `max`.
fn length(stream: &mut impl Read, max: usize) -> io::Result<Option<usize>> {
    assert!(max <= HEADER_TEXT_BYTES_MAX, "a header length cap is at most the text's");
    let mut word = [0u8; 4];
    read_exact(stream, &mut word)?;
    let value = usize::try_from(i32::from_le_bytes(word)).ok().filter(|&v| v <= max);
    assert!(value.is_none_or(|v| v <= max), "a length is within its cap");
    Ok(value)
}

/// Checks the magic and skips the header text and reference names.
fn head_sequences_header(
    path: &Path,
    stream: &mut impl Read,
) -> Result<(Option<String>, Option<String>), ReadsRefusal> {
    let truncated = |_| refuse(path, "the header is truncated");
    let mut magic = [0u8; 4];
    read_exact(stream, &mut magic).map_err(truncated)?;
    if magic != MAGIC {
        return Err(refuse(path, "not a BAM file: no BAM magic"));
    }
    let text = length(stream, HEADER_TEXT_BYTES_MAX).map_err(truncated)?;
    let text = text.ok_or_else(|| refuse(path, "a header text length out of range"))?;
    let scanned = text.min(READ_GROUP_SCAN_BYTES);
    let mut start = vec![0u8; scanned];
    read_exact(stream, &mut start).map_err(truncated)?;
    let tags = read_group_tags(&start);
    skip(stream, text - scanned).map_err(truncated)?;
    let references = length(stream, REFERENCES_MAX).map_err(truncated)?;
    let references = references.ok_or_else(|| refuse(path, "a reference count out of range"))?;
    for _ in 0..references {
        let name = length(stream, REFERENCE_NAME_BYTES_MAX).map_err(truncated)?;
        let name = name.ok_or_else(|| refuse(path, "a reference name length out of range"))?;
        skip(stream, name + 4).map_err(truncated)?; // The name, then the reference's length.
    }
    assert!(text <= HEADER_TEXT_BYTES_MAX, "the header text is within its cap");
    assert!(references <= REFERENCES_MAX, "the references are within their cap");
    Ok(tags)
}

/// Reads and discards `count` bytes.
fn skip(stream: &mut impl Read, count: usize) -> io::Result<()> {
    let wanted = u64::try_from(count).unwrap_or(u64::MAX);
    let skipped = io::copy(&mut stream.by_ref().take(wanted), &mut io::sink())?;
    if skipped != wanted {
        return Err(io::Error::from(io::ErrorKind::UnexpectedEof));
    }
    assert_eq!(skipped, wanted, "every byte asked for was skipped");
    assert!(count <= HEADER_TEXT_BYTES_MAX + 4, "a skip is within the header's caps");
    Ok(())
}

/// Reads record `index` into `record`; false at a clean end of file (between records).
fn head_sequences_record(
    path: &Path,
    stream: &mut impl Read,
    index: usize,
    record: &mut Vec<u8>,
) -> Result<bool, ReadsRefusal> {
    let mut word = [0u8; 4];
    let mut filled = 0;
    while filled < word.len() {
        let Some(rest) = word.get_mut(filled..) else { unreachable!("filled < 4") };
        match stream.read(rest) {
            Ok(0) if filled == 0 => return Ok(false),
            Ok(0) => return Err(record_refusal(path, index, "the file ends inside a record")),
            Ok(n) => filled += n,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(source) => return Err(ReadsRefusal::Read { source }),
        }
    }
    let size = usize::try_from(i32::from_le_bytes(word)).ok();
    let size = size.filter(|&s| (RECORD_FIXED_BYTES..=RECORD_BYTES_MAX).contains(&s));
    let size = size.ok_or_else(|| record_refusal(path, index, "a record length out of range"))?;
    record.clear();
    record.resize(size, 0);
    read_exact(stream, record)
        .map_err(|_| record_refusal(path, index, "the file ends inside a record"))?;
    assert_eq!(record.len(), size, "the record holds its stated length");
    assert!(size >= RECORD_FIXED_BYTES, "a record holds its fixed fields");
    Ok(true)
}

fn record_refusal(path: &Path, record: usize, reason: &'static str) -> ReadsRefusal {
    assert!(record >= 1, "records are numbered from 1");
    assert!(!reason.is_empty(), "a refusal says why");
    ReadsRefusal::BamRecord { path: path.to_path_buf(), record, reason }
}

/// The `N` bytes at `at` among a record's fixed fields, which every record holds.
fn field<const N: usize>(record: &[u8], at: usize) -> [u8; N] {
    assert!(at + N <= RECORD_FIXED_BYTES, "a fixed field lies in the fixed part");
    let Some(bytes) = record.get(at..at + N).and_then(|b| <[u8; N]>::try_from(b).ok()) else {
        unreachable!("a record holds its {RECORD_FIXED_BYTES} fixed bytes, checked on reading")
    };
    assert_eq!(bytes.len(), N, "a field is its width");
    bytes
}

/// Decodes a record's bases into `sequence`; false for a record that repeats a read.
fn head_sequences_bases(
    path: &Path,
    index: usize,
    record: &[u8],
    sequence: &mut Vec<u8>,
) -> Result<bool, ReadsRefusal> {
    let name = usize::from(field::<1>(record, 8)[0]);
    let operations = usize::from(u16::from_le_bytes(field::<2>(record, 12)));
    let flags = u16::from_le_bytes(field::<2>(record, 14));
    if flags & FLAG_REPEATS_READ != 0 {
        return Ok(false);
    }
    let bases = usize::try_from(i32::from_le_bytes(field::<4>(record, 16)));
    let bases = bases.ok().filter(|&b| b <= FASTQ_LINE_BYTES_MAX);
    let bases = bases.ok_or_else(|| record_refusal(path, index, "a read length out of range"))?;
    let start = RECORD_FIXED_BYTES + name + 4 * operations;
    let packed = record
        .get(start..start + bases.div_ceil(2))
        .ok_or_else(|| record_refusal(path, index, "the bases run past the record"))?;
    sequence.clear();
    sequence.extend(
        packed
            .iter()
            .flat_map(|&byte| [byte >> 4, byte & 0xf])
            .take(bases)
            .map(|code| BASES.get(usize::from(code)).copied().unwrap_or(b'N')),
    );
    assert_eq!(sequence.len(), bases, "every base is decoded");
    assert!(start + packed.len() <= record.len(), "the bases lie inside the record");
    Ok(true)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use std::io::Write;

    /// One record: name, flags, bases (4-bit), qualities; no CIGAR, no tags.
    fn record(name: &str, flags: u16, bases: &[u8]) -> Vec<u8> {
        let code = |b: u8| u8::try_from(BASES.iter().position(|&c| c == b).unwrap()).unwrap();
        let mut body = Vec::new();
        body.extend((-1i32).to_le_bytes()); // reference
        body.extend((-1i32).to_le_bytes()); // position
        body.push(u8::try_from(name.len() + 1).unwrap());
        body.push(255); // mapping quality
        body.extend(4680u16.to_le_bytes()); // bin
        body.extend(0u16.to_le_bytes()); // CIGAR operations
        body.extend(flags.to_le_bytes());
        body.extend(i32::try_from(bases.len()).unwrap().to_le_bytes());
        body.extend((-1i32).to_le_bytes());
        body.extend((-1i32).to_le_bytes());
        body.extend(0i32.to_le_bytes());
        body.extend(name.as_bytes());
        body.push(0);
        for pair in bases.chunks(2) {
            body.push(code(pair[0]) << 4 | pair.get(1).map_or(0, |&b| code(b)));
        }
        body.extend(std::iter::repeat_n(30u8, bases.len()));
        let mut out = i32::try_from(body.len()).unwrap().to_le_bytes().to_vec();
        out.extend(body);
        out
    }

    fn bam(records: &[Vec<u8>]) -> Vec<u8> {
        let mut raw = MAGIC.to_vec();
        let text = b"@HD\tVN:1.6\n";
        raw.extend(i32::try_from(text.len()).unwrap().to_le_bytes());
        raw.extend(text);
        raw.extend(1i32.to_le_bytes()); // one reference
        raw.extend(5i32.to_le_bytes());
        raw.extend(b"chr1\0");
        raw.extend(1000i32.to_le_bytes());
        for r in records {
            raw.extend(r);
        }
        raw
    }

    fn sample(raw: &[u8], limit: usize) -> Result<Reads, ReadsRefusal> {
        // Through gzip, as BAM is stored.
        let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
        encoder.write_all(raw).unwrap();
        let compressed = encoder.finish().unwrap();
        let mut stream = flate2::read::MultiGzDecoder::new(&compressed[..]);
        head_sequences(Path::new("t.bam"), &mut stream, limit)
    }

    #[test]
    fn reads_bases_and_skips_repeated_reads() {
        let raw = bam(&[
            record("a", 0, b"ACGTN"),
            record("b", 0x100, b"GGGG"),
            record("c", 0x800, b"TTTT"),
            record("d", 0x4, b"CA"),
        ]);
        let reads = sample(&raw, 10).unwrap();
        assert_eq!(reads.iter().collect::<Vec<_>>(), vec![&b"ACGTN"[..], &b"CA"[..]]);
        assert_eq!(sample(&raw, 1).unwrap().len(), 1);
    }

    #[test]
    fn refuses_broken_files() {
        let good = bam(&[record("a", 0, b"ACGT")]);
        let reason = |raw: &[u8]| match sample(raw, 10) {
            Err(
                ReadsRefusal::BamHeader { reason, .. } | ReadsRefusal::BamRecord { reason, .. },
            ) => reason,
            other => panic!("expected a BAM refusal, got {:?}", other.map(|r| r.len())),
        };
        assert_eq!(reason(&good[..good.len() - 3]), "the file ends inside a record");
        assert_eq!(reason(&good[..10]), "the header is truncated");
        let mut bad = good.clone();
        bad[3] = 2;
        assert_eq!(reason(&bad), "not a BAM file: no BAM magic");
        let mut bad = good.clone();
        bad[4..8].copy_from_slice(&(-5i32).to_le_bytes());
        assert_eq!(reason(&bad), "a header text length out of range");
        let mut bad = bam(&[]);
        bad.extend(8i32.to_le_bytes());
        bad.extend([0u8; 8]);
        assert_eq!(reason(&bad), "a record length out of range");
    }

    /// Reading allocates the same for 100 and 1,000 records: buffers are sized once, and only the
    /// first record's name is kept.
    #[test]
    fn allocates_independently_of_read_count() {
        let file =
            |n: usize| bam(&(0..n).map(|_| record("*", 0, b"ACGTACGTAC")).collect::<Vec<_>>());
        let (small, large) = (file(100), file(1_000));
        let count = |raw: &[u8], n: usize| {
            allocation_counter::measure(|| {
                std::hint::black_box(sample(raw, n).unwrap());
            })
            .count_total
        };
        assert_eq!(count(&small, 100), count(&large, 1_000), "the BAM reader allocates per read");
    }

    /// Chaos probe: a newline or an interior NUL in a read name aborted the process. A name is
    /// file bytes, so it is kept with its control characters made visible.
    #[test]
    fn keeps_hostile_names_visible() {
        for (raw, kept) in [
            ("a\nb", "a\u{fffd}b"),
            ("a\0b", "a\u{fffd}b"),
            ("a\x1b[31mb", "a\u{fffd}[31mb"),
            ("a\tb", "a\tb"),
        ] {
            let reads = sample(&bam(&[record(raw, 0, b"ACGT")]), 10).unwrap();
            assert_eq!(reads.origin().first_name.as_deref(), Some(kept), "{raw:?}");
        }
    }
}
