// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Shared 2-bit k-mer packing, and the one table of reference 25-mers the cDNA read is checked
//! against: chrM forward and reverse, housekeeping genes sense and antisense, TCR/BCR constant
//! regions. One probe per window answers every set, and one pass over the reads feeds both the
//! strand check and library typing.

use crate::constants::{KMER_LENGTH, KMER_STEP, KMER_WINDOWS_MIN, LOW_COMPLEXITY_RUN};
use crate::detect::profile::base_index;
use crate::hash::{FxHashMap, FxHashSet};
use crate::input::fastq::Reads;

/// Bases a u64 holds at 2 bits each.
pub const PACKED_BASES_MAX: usize = 32;
const _: () = assert!(2 * PACKED_BASES_MAX == u64::BITS as usize, "two bits per base fill a u64");

pub fn complement(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        _ => b'N',
    }
}

/// 2-bit packs up to 32 bases; None when any base is not A/C/G/T.
pub fn pack_kmer(k: &[u8]) -> Option<u64> {
    assert!(k.len() <= PACKED_BASES_MAX, "k-mer of {} bases does not fit a u64", k.len());
    let packed = k.iter().try_fold(0u64, |acc, &b| match base_index(b) {
        i @ 0..=3 => Some((acc << 2) | i as u64),
        _ => None,
    });
    if let Some(value) = packed
        && k.len() < 32
    {
        assert!(value >> (2 * k.len()) == 0, "{} bases packed past their 2 bits each", k.len());
    }
    packed
}

/// True when a k-mer holds a run of `run` identical bases: such k-mers (A-rich transcript ends)
/// match polyT stretches on the other strand, so they never vote.
pub fn low_complexity(k: &[u8], run: usize) -> bool {
    // A run of one base is every base, so the test only means something from two.
    assert!(run >= 2, "a homopolymer run is at least 2 bases, not {run}");
    let mut n = 1;
    let found = k.windows(2).any(|w| {
        n = if w.first() == w.last() { n + 1 } else { 1 };
        n >= run
    });
    if found {
        assert!(k.len() >= run, "a run of {run} found in only {} bases", k.len());
    }
    found
}

const _: () = assert!(
    KMER_LENGTH <= PACKED_BASES_MAX && KMER_STEP >= 1 && KMER_STEP <= KMER_LENGTH,
    "a 25-mer packs into a u64 and windows overlap"
);

/// Bits of a table entry: which reference sets hold the 25-mer.
const MITO_FORWARD: u8 = 1 << 0;
const MITO_REVERSE: u8 = 1 << 1;
const HOUSEKEEPING_SENSE: u8 = 1 << 2;
const HOUSEKEEPING_ANTISENSE: u8 = 1 << 3;
const CONSTANT: u8 = 1 << 4;
const SETS: usize = 5;

/// Every reference 25-mer, mapped to the sets that hold it:
/// - chrM (human and mouse) forward and reverse complement; most mitochondrial transcripts come
///   from one strand,
/// - ~40 highly expressed human nuclear genes' spliced transcripts, sense and antisense, with
///   low-complexity and both-strand k-mers removed,
/// - TCR/BCR constant regions, either strand.
fn table() -> &'static FxHashMap<u64, u8> {
    static TABLE: std::sync::OnceLock<FxHashMap<u64, u8>> = std::sync::OnceLock::new();
    TABLE.get_or_init(|| {
        let mut table: FxHashMap<u64, u8> = FxHashMap::default();
        let mut mark = |kmers: &mut dyn Iterator<Item = u64>, bit: u8| {
            for kmer in kmers {
                *table.entry(kmer).or_insert(0) |= bit;
            }
        };
        for genome in [
            &include_bytes!("../../data/reference/human_chrM.seq")[..],
            &include_bytes!("../../data/reference/mouse_chrM.seq")[..],
        ] {
            let reverse_complement: Vec<u8> = genome.iter().rev().map(|&b| complement(b)).collect();
            mark(&mut genome.windows(KMER_LENGTH).filter_map(pack_kmer), MITO_FORWARD);
            mark(&mut reverse_complement.windows(KMER_LENGTH).filter_map(pack_kmer), MITO_REVERSE);
        }
        let (sense, anti) =
            fasta_kmers(include_bytes!("../../data/reference/housekeeping_human.fa"), KMER_LENGTH);
        mark(&mut sense.into_iter(), HOUSEKEEPING_SENSE);
        mark(&mut anti.into_iter(), HOUSEKEEPING_ANTISENSE);
        let (sense, anti) = fasta_kmers(
            include_bytes!("../../data/reference/tcr_bcr_constant_human.fa"),
            KMER_LENGTH,
        );
        mark(&mut sense.into_iter().chain(anti), CONSTANT);
        assert!(table.values().all(|&bits| bits != 0), "every entry names at least one set");
        assert!(table.values().all(|&bits| bits >> SETS == 0), "entries name only the five sets");
        // A missing or empty data file would silently drop a whole set and every call it feeds.
        let named = table.values().fold(0u8, |all, &bits| all | bits);
        assert_eq!(named, (1 << SETS) - 1, "every reference set holds at least one 25-mer");
        table
    })
}

/// Per-read reference hits over the cDNA read: how many reads have at least two windows in each
/// set. Library typing looks at the first `crate::constants::LIBRARY_READS` reads only.
pub struct CdnaHits {
    pub reads: usize,
    pub mito_forward: usize,
    pub mito_reverse: usize,
    pub nuclear_sense: usize,
    pub nuclear_antisense: usize,
    /// Reads library typing looked at, and of those, reads with constant-region or (either strand)
    /// housekeeping hits.
    pub lib_reads: usize,
    pub constant: usize,
    pub housekeeping: usize,
}

impl CdnaHits {
    pub fn of(cdna: &Reads) -> Self {
        let table = table();
        let mut out = Self {
            reads: 0,
            mito_forward: 0,
            mito_reverse: 0,
            nuclear_sense: 0,
            nuclear_antisense: 0,
            lib_reads: 0,
            constant: 0,
            housekeeping: 0,
        };
        for sequence in cdna.iter() {
            let mut hits = [0usize; SETS];
            for_each_window(sequence, |kmer| {
                let bits = table.get(&kmer).copied().unwrap_or(0);
                // Hot: checked build only.
                debug_assert!(bits >> SETS == 0, "a table entry names only the five sets");
                for (set, count) in hits.iter_mut().enumerate() {
                    *count += usize::from(bits >> set & 1 == 1);
                }
            });
            let [mito_forward, mito_reverse, housekeeping_sense, housekeeping_antisense, constant] =
                hits;
            out.reads += 1;
            out.mito_forward += usize::from(mito_forward >= KMER_WINDOWS_MIN);
            out.mito_reverse += usize::from(mito_reverse >= KMER_WINDOWS_MIN);
            out.nuclear_sense += usize::from(housekeeping_sense >= KMER_WINDOWS_MIN);
            out.nuclear_antisense += usize::from(housekeeping_antisense >= KMER_WINDOWS_MIN);
            if out.lib_reads < crate::constants::LIBRARY_READS {
                out.lib_reads += 1;
                out.constant += usize::from(constant >= KMER_WINDOWS_MIN);
                // Either strand: library typing asks how much housekeeping signal there is.
                out.housekeeping +=
                    usize::from(housekeeping_sense + housekeeping_antisense >= KMER_WINDOWS_MIN);
            }
        }
        assert_eq!(out.reads, cdna.len(), "one tally per read");
        assert_eq!(
            out.lib_reads,
            out.reads.min(crate::constants::LIBRARY_READS),
            "library typing sees the first reads"
        );
        for (count, of) in [
            (out.mito_forward, out.reads),
            (out.mito_reverse, out.reads),
            (out.nuclear_sense, out.reads),
            (out.nuclear_antisense, out.reads),
            (out.constant, out.lib_reads),
            (out.housekeeping, out.lib_reads),
        ] {
            assert!(count <= of, "{count} hit reads out of {of}");
        }
        out
    }
}

/// Calls `f` with every 25-mer of `sequence` starting at a multiple of `KMER_STEP` and holding only
/// A/C/G/T: the windows of `sequence.windows(KMER_LENGTH).step_by(KMER_STEP)` that `pack_kmer`
/// packs, packed by rolling each base in once instead of re-reading 25.
fn for_each_window(sequence: &[u8], mut f: impl FnMut(u64)) {
    const MASK: u64 = if KMER_LENGTH == 32 { u64::MAX } else { (1 << (2 * KMER_LENGTH)) - 1 };
    let (mut value, mut valid, mut next_start) = (0u64, 0usize, 0usize);
    for (end, &base) in (1usize..).zip(sequence) {
        match base_index(base) {
            i @ 0..=3 => {
                value = ((value << 2) | i as u64) & MASK;
                valid += 1;
            }
            _ => valid = 0,
        }
        // Hot: checked build only. The run of clean bases ends here, so it is at most `end`.
        debug_assert!(valid <= end, "{valid} clean bases in the first {end}");
        // The window ending here starts at end - KMER_LENGTH; windows start only at multiples of
        // KMER_STEP.
        if end >= KMER_LENGTH && end - KMER_LENGTH == next_start {
            if valid >= KMER_LENGTH {
                f(value);
            }
            next_start += KMER_STEP;
        }
        // Hot: checked build only. The next window to emit never starts past this one's start.
        debug_assert!(
            next_start <= end.saturating_sub(KMER_LENGTH) + KMER_STEP,
            "skipped a window start"
        );
    }
    // Every window start up to the last whole window was visited once.
    assert!(
        next_start <= sequence.len().saturating_sub(KMER_LENGTH) + KMER_STEP,
        "window starts past the read"
    );
    if sequence.len() >= KMER_LENGTH {
        assert!(
            next_start > sequence.len() - KMER_LENGTH,
            "the last window start {} was not reached",
            sequence.len() - KMER_LENGTH
        );
    }
}

/// Sense and antisense k-mer sets of every record in an embedded FASTA, low-complexity k-mers and
/// k-mers present on both strands removed.
pub fn fasta_kmers(fasta: &[u8], k: usize) -> (FxHashSet<u64>, FxHashSet<u64>) {
    assert!(k >= 1, "a k-mer has at least one base");
    assert!(k <= PACKED_BASES_MAX, "a {k}-mer does not pack into a u64");
    let (mut sense, mut anti) = (FxHashSet::default(), FxHashSet::default());
    let mut sequence: Vec<u8> = Vec::new();
    let flush = |sequence: &mut Vec<u8>, sense: &mut FxHashSet<u64>, anti: &mut FxHashSet<u64>| {
        let reverse_complement: Vec<u8> = sequence.iter().rev().map(|&b| complement(b)).collect();
        for (s, set) in [(&sequence[..], &mut *sense), (&reverse_complement[..], &mut *anti)] {
            set.extend(
                s.windows(k)
                    .filter(|w| !low_complexity(w, LOW_COMPLEXITY_RUN))
                    .filter_map(pack_kmer),
            );
        }
        sequence.clear();
    };
    for line in fasta.split(|&b| b == b'\n') {
        if line.first() == Some(&b'>') {
            flush(&mut sequence, &mut sense, &mut anti);
        } else {
            sequence.extend(
                line.iter().map(u8::to_ascii_uppercase).filter(|b| !b.is_ascii_whitespace()),
            );
        }
    }
    flush(&mut sequence, &mut sense, &mut anti);
    let both: FxHashSet<u64> = sense.intersection(&anti).copied().collect();
    sense.retain(|x| !both.contains(x));
    anti.retain(|x| !both.contains(x));
    // Once at startup over the embedded references: a k-mer on both strands cannot vote.
    assert!(sense.is_disjoint(&anti), "a k-mer kept on both strands");
    if fasta.contains(&b'>') {
        assert!(!sense.is_empty(), "an embedded FASTA with records yields k-mers");
    }
    (sense, anti)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Rolling windows are exactly the windows the straightforward packing yields, Ns included.
    #[test]
    fn rolling_windows_match_packing_each_window() {
        let mut state = 5u64;
        for len in [0, 1, 24, 25, 26, 29, 30, 31, 90, 151] {
            for _ in 0..200 {
                let sequence: Vec<u8> = (0..len)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        b"ACGTACGTACGTACGTN"[usize::try_from(state >> 60).unwrap()]
                    })
                    .collect();
                let expected: Vec<u64> = sequence
                    .windows(KMER_LENGTH)
                    .step_by(KMER_STEP)
                    .filter_map(pack_kmer)
                    .collect();
                let mut rolled = Vec::new();
                for_each_window(&sequence, |v| rolled.push(v));
                assert_eq!(rolled, expected, "{}", String::from_utf8_lossy(&sequence));
            }
        }
    }
}
