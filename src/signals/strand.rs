// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Which strand the cDNA read is on, without a genome: homopolymer runs and the built-in
//! mitochondrial genomes. A sense cDNA read means a 3' kit, an antisense one a 5' kit.

use crate::constants::{
    HOMOPOLYMER_RUN, STRAND_RATIO_MIN, STRAND_SHARE_MIN, STRAND_STRONG_RATIO_MIN,
    STRAND_STRONG_SHARE_MIN,
};
use crate::detect::profile::ratio;
use crate::input::fastq::Reads;
use crate::signals::kmers::CdnaHits;

/// Which strand the cDNA read is on, from homopolymer runs: a 3' kit's cDNA read is sense and runs
/// into the polyA tail (AAAA); a 5' kit's is antisense and shows the same tail as TTTT. The
/// mitochondrial and nuclear shares come from the built-in references in `data/reference`.
pub struct StrandEvidence {
    pub poly_a: f64,
    pub poly_t: f64,
    /// Share of reads matching the built-in human/mouse chrM forward strand, and its reverse
    /// complement. Most mitochondrial transcripts come from one strand, so a sense cDNA read
    /// (3') matches forward and an antisense one (5') matches the reverse complement.
    pub mito_forward: f64,
    pub mito_reverse: f64,
    /// Share of reads matching ~40 highly expressed human nuclear genes' spliced transcripts in
    /// sense, and antisense (`data/reference/housekeeping_human.fa`). Calibrated on known 10x
    /// data: 3' reads 5-7x sense, 5' reads 24x antisense.
    pub nuclear_sense: f64,
    pub nuclear_antisense: f64,
}

/// A strand call from one signal: the leading side must be 3x the other and at least 0.5% of
/// reads; strong at 10x and 1%.
fn lean(sense: f64, antisense: f64) -> Option<(Strand, bool)> {
    assert!((0.0..=1.0).contains(&sense), "sense share {sense} outside 0..=1");
    assert!((0.0..=1.0).contains(&antisense), "antisense share {antisense} outside 0..=1");
    let (hi, lo, strand) = if sense >= antisense {
        (sense, antisense, Strand::Sense)
    } else {
        (antisense, sense, Strand::Antisense)
    };
    let call = (hi >= STRAND_SHARE_MIN && hi >= STRAND_RATIO_MIN * lo)
        .then_some((strand, hi >= STRAND_STRONG_SHARE_MIN && hi >= STRAND_STRONG_RATIO_MIN * lo));
    if let Some((_, strong)) = call {
        assert!(hi > lo, "a strand call leans to the larger share: {hi} over {lo}");
        if strong {
            assert!(
                hi >= STRAND_STRONG_SHARE_MIN,
                "a strong call stands on STRAND_STRONG_SHARE_MIN, not {hi}"
            );
        }
    }
    call
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Strand {
    Sense,
    Antisense,
}

impl StrandEvidence {
    /// `hits` is `CdnaHits::of(sequences)`, shared with library typing.
    pub fn of(sequences: &Reads, hits: &CdnaHits) -> Self {
        assert_eq!(hits.reads, sequences.len(), "reference hits must come from these reads");
        // Paired with CdnaHits::of's postcondition: no set counts more reads than were tallied.
        for count in
            [hits.mito_forward, hits.mito_reverse, hits.nuclear_sense, hits.nuclear_antisense]
        {
            assert!(count <= hits.reads, "{count} hit reads out of {}", hits.reads);
        }
        let has_run = |sequence: &[u8], base: u8| {
            let mut run = 0;
            sequence.iter().any(|&b| {
                run = if b == base { run + 1 } else { 0 };
                run >= HOMOPOLYMER_RUN
            })
        };
        let share = |f: &dyn Fn(&[u8]) -> bool| {
            ratio(sequences.iter().filter(|s| f(s)).count(), sequences.len())
        };
        let evidence = Self {
            poly_a: share(&|s| has_run(s, b'A')),
            poly_t: share(&|s| has_run(s, b'T')),
            mito_forward: ratio(hits.mito_forward, sequences.len()),
            mito_reverse: ratio(hits.mito_reverse, sequences.len()),
            nuclear_sense: ratio(hits.nuclear_sense, sequences.len()),
            nuclear_antisense: ratio(hits.nuclear_antisense, sequences.len()),
        };
        // A read can show a polyA run and a polyT run, but each share is still a part of the reads.
        assert!(evidence.poly_a <= 1.0, "polyA share {} past every read", evidence.poly_a);
        assert!(evidence.poly_t <= 1.0, "polyT share {} past every read", evidence.poly_t);
        evidence
    }

    pub fn poly_call(&self) -> Option<(Strand, bool)> {
        lean(self.poly_a, self.poly_t)
    }

    pub fn mito_call(&self) -> Option<(Strand, bool)> {
        lean(self.mito_forward, self.mito_reverse)
    }

    pub fn nuclear_call(&self) -> Option<(Strand, bool)> {
        lean(self.nuclear_sense, self.nuclear_antisense)
    }

    /// The combined strand call over the signals that lean: any disagreement gives none; two or
    /// more agreeing give a strong call; one alone keeps its own strength.
    pub fn call(&self) -> Option<(Strand, bool)> {
        let leans: Vec<(Strand, bool)> = [self.poly_call(), self.mito_call(), self.nuclear_call()]
            .into_iter()
            .flatten()
            .collect();
        assert!(leans.len() <= 3, "three signals give at most three leans");
        let (first, rest) = leans.split_first()?;
        if rest.iter().any(|(s, _)| *s != first.0) {
            return None;
        }
        let strong = leans.len() >= 2 || first.1;
        assert!(leans.iter().all(|(s, _)| *s == first.0), "a combined call has every lean agree");
        if !strong {
            assert_eq!(leans.len(), 1, "only a lone weak lean gives a weak combined call");
        }
        Some((first.0, strong))
    }
}
