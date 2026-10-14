// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! What kind of library sits behind a matched bead layout. The same 10x beads carry gene
//! expression, immune receptor (VDJ), CRISPR guide and tag or amplicon libraries, and each needs a
//! different pipeline. Each type leaves a mark in the cDNA read:
//!
//! - VDJ: many reads from TCR/BCR constant regions and few from housekeeping genes (receptor
//!   transcripts are enriched; a gene-expression library rich in plasma cells has both),
//! - CRISPR guide: the sgRNA scaffold (standard or optimised variant),
//! - tag or amplicon: almost every read opens with one of a small panel of sequences (antibody or
//!   hashtag tags, or a guide or barcode amplicon read from before its scaffold).
//!
//! Thresholds calibrated on 7 Oct 2026 over 480 public scBaseCount samples
//! (`bench/library_type`), library type from ENA titles: VDJ median 20% constant-region reads at
//! 0.9% housekeeping; gene expression 0.01% (95th percentile 2.7%), housekeeping median 11%; guide
//! libraries >= 24% scaffold reads against 0.3% (95th percentile) for gene expression; tag
//! libraries >= 80% top-200 start cover.

use crate::constants::{
    GUIDE_SHARE_MIN, LIBRARY_READS, SCAFFOLD_KMER_LENGTH, SCAFFOLD_KMER_STEP, TAG_COVER_MIN,
    TAG_HOUSEKEEPING_MAX, TAG_PREFIX_BASES, TAG_STARTS_TOP, VDJ_OVER_HOUSEKEEPING, VDJ_SHARE_MIN,
};
use crate::detect::profile::ratio;
use crate::hash::{FxHashMap, FxHashSet, map_with_capacity};
use crate::input::fastq::Reads;
use crate::signals::kmers::{CdnaHits, PACKED_BASES_MAX, pack_kmer};

const SCAFFOLDS: [&[u8]; 2] = [
    b"GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
    // Optimised (F+E) scaffold of many CROP-seq / Perturb-seq vectors.
    b"GTTTAAGAGCTATGCTGGAAACAGCATAGCAAGTTTAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC",
];

/// What kind of library sits behind a matched bead layout; each needs a different pipeline.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum LibraryType {
    /// Whole-transcriptome gene expression.
    GeneExpression,
    /// Immune receptor (TCR/BCR) enrichment.
    Vdj,
    /// CRISPR guide (sgRNA) capture.
    CrisprGuide,
    /// A tag or amplicon library: almost every read opens with one of a few sequences. Antibody
    /// and hashtag tags read this way, and so do guide or lineage-barcode amplicons whose reads
    /// start before the sgRNA scaffold (a U6 promoter amplicon, `ClonMapper` capture): held-out
    /// Parse runs of both kinds, 8 Oct 2026.
    FeatureTag,
}

impl LibraryType {
    /// A short description for people, e.g. `immune receptor (VDJ: TCR/BCR)`.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            LibraryType::GeneExpression => "gene expression",
            LibraryType::Vdj => "immune receptor (VDJ: TCR/BCR)",
            LibraryType::CrisprGuide => "CRISPR guide capture",
            LibraryType::FeatureTag => "tag or amplicon (antibody, hashtag, guide or barcode)",
        }
    }
}

pub struct LibraryEvidence {
    pub constant_region: f64,
    pub housekeeping: f64,
    pub scaffold: f64,
    pub tag_cover: f64,
}

/// 12-mers of both sgRNA scaffolds, either strand.
fn scaffold_kmers() -> &'static FxHashSet<u64> {
    static SET: std::sync::OnceLock<FxHashSet<u64>> = std::sync::OnceLock::new();
    SET.get_or_init(|| {
        let mut scaffold = FxHashSet::default();
        for s in SCAFFOLDS {
            let reverse_complement: Vec<u8> =
                s.iter().rev().map(|&b| crate::signals::kmers::complement(b)).collect();
            for sequence in [s, &reverse_complement[..]] {
                scaffold.extend(sequence.windows(SCAFFOLD_KMER_LENGTH).filter_map(pack_kmer));
            }
        }
        let windows: usize =
            SCAFFOLDS.iter().map(|s| 2 * (s.len() - SCAFFOLD_KMER_LENGTH + 1)).sum();
        assert!(!scaffold.is_empty(), "the sgRNA scaffolds yield 12-mers");
        assert!(
            scaffold.len() <= windows,
            "{} scaffold 12-mers from {windows} windows",
            scaffold.len()
        );
        scaffold
    })
}

fn windows_in(sequence: &[u8], k: usize, step: usize, set: &FxHashSet<u64>) -> usize {
    assert!((1..=PACKED_BASES_MAX).contains(&k), "a {k}-mer does not pack into a u64");
    assert!(step >= 1, "windows advance by at least one base");
    let found = sequence
        .windows(k)
        .step_by(step)
        .filter(|w| pack_kmer(w).is_some_and(|v| set.contains(&v)))
        .count();
    assert!(found <= sequence.len(), "{found} windows found in a read of {} bases", sequence.len());
    found
}

impl LibraryEvidence {
    /// `hits` is `CdnaHits::of(cdna)`, shared with the strand check: its constant-region and
    /// housekeeping counts cover the same first `LIBRARY_READS` reads as this pass.
    pub fn of(cdna: &Reads, hits: &CdnaHits) -> Self {
        assert_eq!(hits.reads, cdna.len(), "reference hits must come from these reads");
        let scaffold_set = scaffold_kmers();
        let (mut scaffold, mut n) = (0usize, 0usize);
        let mut starts: FxHashMap<u64, usize> = map_with_capacity(LIBRARY_READS);
        for sequence in cdna.iter().take(LIBRARY_READS) {
            n += 1;
            scaffold += usize::from(
                windows_in(sequence, SCAFFOLD_KMER_LENGTH, SCAFFOLD_KMER_STEP, scaffold_set) >= 1,
            );
            if let Some(start) = sequence.get(..TAG_PREFIX_BASES).and_then(pack_kmer) {
                *starts.entry(start).or_insert(0) += 1;
            }
        }
        assert_eq!(n, hits.lib_reads, "library typing and CdnaHits see the same first reads");
        assert!(scaffold <= n, "{scaffold} scaffold reads out of {n}");
        let (constant, housekeeping) = (hits.constant, hits.housekeeping);
        // Paired with CdnaHits::of's postcondition.
        assert!(constant <= n, "{constant} constant-region reads out of {n}");
        assert!(housekeeping <= n, "{housekeeping} housekeeping reads out of {n}");
        let mut counts: Vec<usize> = starts.into_values().collect();
        counts.sort_unstable_by(|a, b| b.cmp(a));
        let top: usize = counts.iter().take(TAG_STARTS_TOP).sum();
        assert!(top <= n, "the top {TAG_STARTS_TOP} read starts cover {top} of {n} reads");
        Self {
            constant_region: ratio(constant, n),
            housekeeping: ratio(housekeeping, n),
            scaffold: ratio(scaffold, n),
            tag_cover: ratio(top, n),
        }
    }

    pub fn call(&self) -> LibraryType {
        let (constant, housekeeping) = (self.constant_region, self.housekeeping);
        assert!((0.0..=1.0).contains(&constant), "constant-region share {constant} outside 0..=1");
        assert!(
            (0.0..=1.0).contains(&housekeeping),
            "housekeeping share {housekeeping} outside 0..=1"
        );
        assert!(
            (0.0..=1.0).contains(&self.scaffold),
            "scaffold share {} outside 0..=1",
            self.scaffold
        );
        assert!(
            (0.0..=1.0).contains(&self.tag_cover),
            "tag cover {} outside 0..=1",
            self.tag_cover
        );
        if self.scaffold >= GUIDE_SHARE_MIN {
            LibraryType::CrisprGuide
        } else if self.tag_cover >= TAG_COVER_MIN && self.housekeeping < TAG_HOUSEKEEPING_MAX {
            LibraryType::FeatureTag
        } else if self.constant_region >= VDJ_SHARE_MIN
            && self.constant_region >= VDJ_OVER_HOUSEKEEPING * self.housekeeping
        {
            LibraryType::Vdj
        } else {
            LibraryType::GeneExpression
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn reads(n: usize, len: usize) -> Reads {
        let mut state = 3u64;
        let mut out = Reads::with_capacity(n, n * len);
        let mut sequence = vec![0u8; len];
        for _ in 0..n {
            for base in &mut sequence {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                *base = b"ACGT"[usize::try_from(state >> 62).unwrap()];
            }
            out.push(&sequence);
        }
        out
    }

    #[test]
    fn scaffold_reads_are_guide_and_random_reads_are_gene_expression() {
        let mut guide = Reads::with_capacity(100, 10_000);
        for _ in 0..100 {
            guide.push(&SCAFFOLDS[1][..60]);
        }
        let of = |r: &Reads| LibraryEvidence::of(r, &CdnaHits::of(r));
        assert_eq!(of(&guide).call(), LibraryType::CrisprGuide);
        assert_eq!(of(&reads(2_000, 90)).call(), LibraryType::GeneExpression);
    }

    /// The pass over reads allocates the same on 1,000 and 10,000 reads (after the one-off build of
    /// the reference sets).
    #[test]
    fn allocates_independently_of_read_count() {
        let one = reads(1, 90);
        std::hint::black_box(LibraryEvidence::of(&one, &CdnaHits::of(&one)));
        let (small, large) = (reads(1_000, 90), reads(10_000, 90));
        let count = |r: &Reads| {
            allocation_counter::measure(|| {
                std::hint::black_box(LibraryEvidence::of(r, &CdnaHits::of(r)));
            })
            .count_total
        };
        assert_eq!(count(&small), count(&large), "library typing allocates per read");
    }
}
