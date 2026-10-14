// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Reference-free test for whether a read file has positional structure (linkers, barcodes, a
//! polyT tail at a consistent place), as any barcode read does even for a kit no layout describes.
//! A pure cDNA read has none: transcripts start everywhere. Used only to word the hint when no
//! layout matches, so it runs on that path alone.
//!
//! Thresholds measured on 7 Oct 2026 (scripts/discover.py's anchor idea): the commonest
//! non-adapter, non-low-complexity 12-mer, counted once per read, covers 0.3-6% of reads in cDNA,
//! bulk RNA, WGS and bisulfite files, and 27-94% of Singleron, inDrops and Parse barcode reads.
//! The `ScaleBio` 7-base linker anchors no 12-mer, but its barcode read has 28 dominant columns
//! against 0-3 for the cDNA-like files.

use crate::constants::{
    ANCHOR_KMER_LENGTH, ANCHOR_READS, ANCHOR_RUN_MAX, ANCHOR_SHARE_MIN, ANCHOR_SPAN,
    FIXED_COLUMNS_MIN,
};
use crate::detect::profile::Profile;
use crate::hash::{FxHashMap, map_with_capacity};
use crate::input::fastq::Reads;

/// Adapter and template-switch sequences: read-through and SMART cDNA ends, never an anchor.
const ADAPTERS: [&[u8]; 6] = [
    b"CTGTCTCTTATACACATCT",
    b"AGATGTGTATAAGAGACAG",
    b"AGATCGGAAGAGC",
    b"GCTCTTCCGATCT",
    // SMART template-switch oligo, both strands: it ends many cDNA reads.
    b"AAGCAGTGGTATCAACGCAGAGTACATGGG",
    b"CCCATGTACTCTGCGTTGATACCACTGCTT",
];

pub struct StructureSignal {
    pub fixed_columns: usize,
    /// The commonest qualifying 12-mer and the share of scanned reads carrying it.
    pub anchor: Option<(String, f64)>,
}

impl StructureSignal {
    /// Strength for choosing which file to name: the anchor share, or the fixed-column share.
    pub fn strength(&self) -> f64 {
        let anchor = self.anchor.as_ref().map_or(0.0, |(_, s)| *s);
        assert!((0.0..=1.0).contains(&anchor), "anchor share {anchor} outside 0..=1");
        let columns =
            crate::detect::profile::ratio(self.fixed_columns, FIXED_COLUMNS_MIN) * ANCHOR_SHARE_MIN;
        let strength = anchor.max(columns);
        assert!(strength.is_finite(), "structure strength {strength} is not a number");
        strength
    }

    pub fn structured(&self) -> bool {
        let structured = self.fixed_columns >= FIXED_COLUMNS_MIN
            || self.anchor.as_ref().is_some_and(|(_, s)| *s >= ANCHOR_SHARE_MIN);
        assert!(
            self.fixed_columns <= crate::constants::PROFILE_POSITIONS,
            "{} fixed columns out of {} profiled",
            self.fixed_columns,
            crate::constants::PROFILE_POSITIONS
        );
        // Paired with strength(): the bar for structure is exactly a strength of
        // ANCHOR_SHARE_MIN, so the file a hint names is always one that passed it.
        assert_eq!(
            structured,
            self.strength() >= ANCHOR_SHARE_MIN,
            "structured and strength disagree"
        );
        structured
    }
}

fn low_complexity(kmer: &[u8]) -> bool {
    assert_eq!(kmer.len(), ANCHOR_KMER_LENGTH, "anchors are {ANCHOR_KMER_LENGTH}-mers");
    let mut run = 1;
    let found = kmer.windows(2).any(|w| {
        run = if w.first() == w.last() { run + 1 } else { 1 };
        run > ANCHOR_RUN_MAX
    });
    if found {
        assert!(run > ANCHOR_RUN_MAX, "a low-complexity k-mer holds a run past {ANCHOR_RUN_MAX}");
    }
    found
}

/// Shortest adapter: a k-mer longer than it could never be found, silently.
const ADAPTER_LEN_MIN: usize = 13;
const _: () = assert!(ANCHOR_KMER_LENGTH <= ADAPTER_LEN_MIN, "an anchor fits inside every adapter");
#[allow(clippy::indexing_slicing)] // Const evaluation: an index past the array fails the build.
const _: () = {
    let mut index = 0;
    while index < ADAPTERS.len() {
        assert!(ADAPTERS[index].len() >= ADAPTER_LEN_MIN, "an adapter is shorter than the minimum");
        index += 1;
    }
};

/// Distinct anchors a read adds to the count map, on average, for sizing it up front: most of a
/// cDNA read's anchors are shared with other reads, so a few per read.
const ANCHORS_PER_READ_GUESS: usize = 8;

/// Whether a k-mer is bases only. Reads keep what the file holds, so an IUPAC code (R, Y, K...)
/// can sit where N usually does; an anchor holding one would be no sequence a kit writes.
fn acgt(kmer: &[u8]) -> bool {
    kmer.iter().all(|b| b"ACGT".contains(b))
}

fn in_adapter(kmer: &[u8]) -> bool {
    assert_eq!(kmer.len(), ANCHOR_KMER_LENGTH, "anchors are {ANCHOR_KMER_LENGTH}-mers");
    // The caller keeps only ACGT k-mers; adapters hold nothing else, so other bases here are a
    // filter bug.
    assert!(acgt(kmer), "a k-mer with a base outside ACGT reached the adapter scan");
    ADAPTERS.iter().any(|a| a.windows(kmer.len()).any(|w| w == kmer))
}

pub fn signal(reads: &Reads, profile: &Profile) -> StructureSignal {
    assert_eq!(profile.reads, reads.len(), "the profile must come from these reads");
    let mut counts: FxHashMap<&[u8], usize> =
        map_with_capacity(ANCHOR_READS * ANCHORS_PER_READ_GUESS);
    // At most one anchor per window of the scanned span.
    let mut seen: Vec<&[u8]> = Vec::with_capacity(ANCHOR_SPAN - ANCHOR_KMER_LENGTH + 1);
    let mut scanned = 0usize;
    for sequence in reads.iter().take(ANCHOR_READS) {
        scanned += 1;
        seen.clear();
        let head = sequence.get(..ANCHOR_SPAN).unwrap_or(sequence);
        seen.extend(
            head.windows(ANCHOR_KMER_LENGTH)
                .filter(|kmer| acgt(kmer) && !low_complexity(kmer) && !in_adapter(kmer)),
        );
        // Counted once per read: a linker repeated in one read still marks one read.
        seen.sort_unstable();
        seen.dedup();
        for kmer in &seen {
            *counts.entry(kmer).or_insert(0) += 1;
        }
    }
    assert!(scanned <= ANCHOR_READS, "scanned {scanned} reads past the anchor budget");
    let anchor =
        counts.iter().max_by(|a, b| a.1.cmp(b.1).then_with(|| b.0.cmp(a.0))).map(|(k, &n)| {
            (String::from_utf8_lossy(k).into_owned(), crate::detect::profile::ratio(n, scanned))
        });
    if let Some((kmer, share)) = &anchor {
        assert_eq!(kmer.len(), ANCHOR_KMER_LENGTH, "the anchor is a {ANCHOR_KMER_LENGTH}-mer");
        assert!((0.0..=1.0).contains(share), "anchor share {share} outside 0..=1");
        assert!(acgt(kmer.as_bytes()), "the anchor {kmer} is ACGT only");
    }
    StructureSignal { fixed_columns: profile.fixed_columns(), anchor }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reads_of(sequences: &[&[u8]]) -> Reads {
        let mut r = Reads::with_capacity(sequences.len(), sequences.iter().map(|s| s.len()).sum());
        for s in sequences {
            r.push(s);
        }
        r
    }

    #[test]
    fn filters_low_complexity_and_adapters() {
        assert!(low_complexity(b"ACGTTTTTTTTA"));
        assert!(!low_complexity(b"ACGTTTTTTTAC"));
        assert!(in_adapter(b"CAACGCAGAGTA"));
        assert!(in_adapter(b"CTGTCTCTTATA"));
        assert!(!in_adapter(b"GTGGCCGATGTT"));
    }

    #[test]
    fn a_shared_linker_is_structure_and_shifting_sequence_is_not() {
        // A linker at varying positions (as BC1 length varies) still anchors every read.
        let linked: Vec<Vec<u8>> = (0..400u32)
            .map(|i| {
                let pad = b"ACGTACGTAC".iter().take(usize::try_from(8 + i % 4).unwrap()).copied();
                let tail = format!("{i:016b}")
                    .into_bytes()
                    .iter()
                    .map(|&b| if b == b'1' { b'G' } else { b'C' })
                    .collect::<Vec<_>>();
                pad.chain(b"GAGTGATTGCTTGTGACGCCTT".iter().copied()).chain(tail).collect()
            })
            .collect();
        let refs: Vec<&[u8]> = linked.iter().map(Vec::as_slice).collect();
        let r = reads_of(&refs);
        let s = signal(&r, &Profile::of(&r));
        assert!(s.structured(), "a linker in every read must count as structure");
        assert!(s.anchor.as_ref().is_some_and(|(_, share)| *share > 0.99));

        // Pseudo-random reads: no k-mer repeats across reads, no fixed column.
        let mut state = 9u64;
        let random: Vec<Vec<u8>> = (0..400)
            .map(|_| {
                (0..60)
                    .map(|_| {
                        state = state
                            .wrapping_mul(6_364_136_223_846_793_005)
                            .wrapping_add(1_442_695_040_888_963_407);
                        b"ACGT"[usize::try_from(state >> 62).unwrap()]
                    })
                    .collect()
            })
            .collect();
        let refs: Vec<&[u8]> = random.iter().map(Vec::as_slice).collect();
        let r = reads_of(&refs);
        assert!(!signal(&r, &Profile::of(&r)).structured(), "random reads have no structure");
    }

    /// An IUPAC code (R, Y, K) in a prefix every read shares is not an anchor: an anchor is bases
    /// only.
    #[test]
    fn an_anchor_is_bases_only() {
        for code in [b'R', b'Y', b'K', b'N'] {
            let shared: Vec<Vec<u8>> = (0..400u32)
                .map(|i| {
                    let bits = format!("{i:016b}");
                    let tail = bits.bytes().map(|b| b"GC"[usize::from(b == b'1')]);
                    b"ACGT"
                        .iter()
                        .copied()
                        .chain([code])
                        .chain(*b"ACGTACGTAC")
                        .chain(tail)
                        .collect()
                })
                .collect();
            let refs: Vec<&[u8]> = shared.iter().map(Vec::as_slice).collect();
            let r = reads_of(&refs);
            let anchor = signal(&r, &Profile::of(&r)).anchor.unwrap();
            assert!(acgt(anchor.0.as_bytes()), "{} anchored {}", char::from(code), anchor.0);
        }
    }
}
