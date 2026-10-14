// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Weighs the evidence in the first reads against every candidate layout:
//!
//! - barcode segments: hits against the kit's list, or, for an unlisted segment, whether a small
//!   set of sequences covers most reads, as real barcodes do,
//! - fixed linkers at their positions,
//! - UMI randomness: every UMI position should be close to 2 bits of entropy,
//! - the tail after the UMI (polyT for 3', TSO for 5') when the read shows it.
//!
//! Each read is first placed at the allowed start offset that best fits the linkers and barcodes.
//! Every layout is scored with the files as labelled and swapped; the better orientation wins.

pub mod decide;
pub mod profile;
pub mod score;

use crate::constants::SUPPORT_MIN;
use crate::detect::decide::{Verdict, decide};
use crate::detect::profile::Profile;
use crate::detect::score::{Lists, Score, score, vouched};
use crate::input::fastq::Reads;
use crate::layouts::{LAYOUTS, LayoutSet, Mate};
use crate::signals::kmers::CdnaHits;
use crate::signals::library::LibraryEvidence;
use crate::signals::strand::StrandEvidence;
use crate::signals::structure::{StructureSignal, signal};

pub struct Report {
    pub swapped: bool,
    /// Strand evidence from the cDNA read: the file the chosen layout's barcode is not on.
    pub cdna_strand: Option<StrandEvidence>,
    pub r1_profile: Profile,
    pub r2_profile: Option<Profile>,
    pub scores: Vec<Score>,
    pub verdict: Verdict,
    /// Per file (R1, then R2), only when no layout matched: does the file have the positional
    /// structure of a barcode read? Words the hint.
    pub structure: Vec<StructureSignal>,
    /// What kind of library the cDNA read shows, when a layout matched or tied and a cDNA read
    /// exists.
    pub library: Option<LibraryEvidence>,
    /// Layouts not scored because a list they read was not fetched (an optional list).
    pub unchecked: LayoutSet,
}

fn score_all(
    r1: &Reads,
    r2: Option<&Reads>,
    p1: &Profile,
    p2: Option<&Profile>,
    lists: &Lists,
) -> Vec<Score> {
    assert_eq!(p1.reads, r1.len(), "the R1 profile is of the R1 reads");
    assert_eq!(p2.map(|p| p.reads), r2.map(Reads::len), "the R2 profile is of the R2 reads");
    let mut scores: Vec<Score> = LAYOUTS
        .iter()
        .filter(|layout| lists.cover(layout))
        .filter_map(|layout| match layout.read {
            Mate::R1 => Some(score(layout, r1, p1, lists)),
            Mate::R2 => Some(score(layout, r2?, p2?, lists)),
        })
        .collect();
    // A layout the reads do not contradict beats one they do; then support, then tail evidence.
    scores.sort_by(|a, b| {
        a.contradictions
            .is_empty()
            .cmp(&b.contradictions.is_empty())
            .reverse()
            .then(b.support.total_cmp(&a.support))
            .then(a.contradictions.len().cmp(&b.contradictions.len()))
            .then(b.tail_rate.unwrap_or(0.0).total_cmp(&a.tail_rate.unwrap_or(0.0)))
    });
    // Postconditions, paired with `decide`: clean layouts first, and every layout scored when
    // both reads are present.
    assert!(
        scores
            .iter()
            .skip_while(|s| s.contradictions.is_empty())
            .all(|s| !s.contradictions.is_empty()),
        "layouts the reads do not contradict sort first"
    );
    if r2.is_some() {
        let covered = LAYOUTS.iter().filter(|l| lists.cover(l)).count();
        assert_eq!(scores.len(), covered, "with both reads, every layout with its lists is scored");
    }
    scores
}

/// How well an orientation explains the data, counting only what a verdict could stand on: the
/// best clean layout, else the best layout a vendor list vouches for. A contradicted unlisted
/// layout scores nothing, so it can neither win the orientation nor flag a swap on its own.
fn standing(scores: &[Score]) -> f64 {
    let best = |it: &mut dyn Iterator<Item = &Score>| {
        it.map(|s| s.support).fold(None, |m: Option<f64>, x| Some(m.map_or(x, |m| m.max(x))))
    };
    let standing = best(&mut scores.iter().filter(|s| s.contradictions.is_empty()))
        .or_else(|| best(&mut scores.iter().filter(|s| vouched(s))))
        .unwrap_or(0.0);
    assert!((0.0..=1.0).contains(&standing), "standing {standing} outside 0..=1");
    if standing > 0.0 {
        assert!(
            scores
                .iter()
                .any(|s| (s.contradictions.is_empty() || vouched(s)) && s.support >= standing),
            "a standing comes from a clean or list-vouched layout"
        );
    }
    standing
}

/// The cDNA read for a layout whose barcode is on `barcode` (as scored), in the files as given:
/// the other file. Swapped scoring exchanged the files, so the barcode sat in the labelled R2.
fn cdna_read<'a>(
    r1: &'a Reads,
    r2: Option<&'a Reads>,
    barcode: Mate,
    swapped: bool,
) -> Option<&'a Reads> {
    assert!(!swapped || r2.is_some(), "a swap needs two files");
    let barcode_in_r1 = (barcode == Mate::R1) != swapped;
    let cdna = if barcode_in_r1 { r2 } else { Some(r1) };
    assert!(r2.is_some() || barcode_in_r1, "one file holds a barcode read, never cDNA alone");
    cdna
}

/// One pass over the cDNA read's 25-mers, which feeds both the strand check and library typing.
struct CdnaPass<'a> {
    cdna: Option<&'a Reads>,
    hits: Option<CdnaHits>,
    strand: Option<StrandEvidence>,
}

impl<'a> CdnaPass<'a> {
    fn of(cdna: Option<&'a Reads>) -> Self {
        let hits = cdna.map(CdnaHits::of);
        let strand = cdna.zip(hits.as_ref()).map(|(c, h)| StrandEvidence::of(c, h));
        assert!(
            cdna.zip(hits.as_ref()).is_none_or(|(c, h)| h.reads <= c.len()),
            "hits count the cDNA file's reads"
        );
        assert!(
            hits.as_ref()
                .is_none_or(|h| h.constant <= h.lib_reads && h.housekeeping <= h.lib_reads),
            "typed reads are among those library typing looked at"
        );
        Self { cdna, hits, strand }
    }
}

pub fn detect(r1: &Reads, r2: Option<&Reads>, lists: &Lists) -> Report {
    assert!(!r1.is_empty(), "detection needs reads: main refuses an empty R1");
    let unchecked = LAYOUTS
        .iter()
        .enumerate()
        .filter(|(_, layout)| !lists.cover(layout))
        .fold(LayoutSet::EMPTY, |set, (index, _)| set.with(index));
    let p1 = Profile::of(r1);
    let p2 = r2.map(Profile::of);
    let as_labelled = score_all(r1, r2, &p1, p2.as_ref(), lists);

    // Swapped files: score with the labels exchanged and keep whichever explains the data better.
    let exchanged = match (r2, p2.as_ref()) {
        (Some(r2), Some(p2)) => Some(score_all(r2, Some(r1), p2, Some(&p1), lists)),
        _ => None,
    };
    let labelled_standing = standing(&as_labelled);
    let (scores, swapped) = match exchanged {
        Some(s) if standing(&s) >= SUPPORT_MIN && standing(&s) > labelled_standing => (s, true),
        _ => (as_labelled, false),
    };
    assert!(!swapped || r2.is_some(), "a swap needs two files");
    // The cDNA read is the file the layout's barcode is not on. The strand can break a tie, so
    // it is measured before the verdict, on the leading layout's cDNA read; tied layouts share
    // its lists, so its read. Only a LOW match can name another layout, read from the other file.
    let lead = scores.first().map_or(Mate::R1, |s| s.layout.read);
    let mut pass = CdnaPass::of(cdna_read(r1, r2, lead, swapped));
    let verdict = decide(&scores, pass.strand.as_ref());
    let chosen = match &verdict {
        Verdict::Match { best, .. } => scores.get(*best).map_or(lead, |s| s.layout.read),
        Verdict::Ambiguous { .. } | Verdict::NoMatch => lead,
    };
    if chosen != lead {
        pass = CdnaPass::of(cdna_read(r1, r2, chosen, swapped));
    }
    let CdnaPass { cdna, hits, strand: cdna_strand } = pass;
    let structure = if matches!(verdict, Verdict::NoMatch) {
        let mut out = vec![signal(r1, &p1)];
        if let (Some(r2), Some(p2)) = (r2, p2.as_ref()) {
            out.push(signal(r2, p2));
        }
        out
    } else {
        Vec::new()
    };
    let library = match verdict {
        Verdict::NoMatch => None,
        _ => cdna.zip(hits.as_ref()).map(|(c, h)| LibraryEvidence::of(c, h)),
    };
    // The report's parts agree: strand evidence exactly when there is a cDNA read, and structure
    // signals (one per file) only when no layout matched.
    assert_eq!(cdna_strand.is_some(), cdna.is_some(), "strand evidence comes from the cDNA read");
    if matches!(verdict, Verdict::NoMatch) {
        assert_eq!(structure.len(), 1 + usize::from(r2.is_some()), "one structure signal per file");
    } else {
        assert!(structure.is_empty(), "structure signals only word a no-match hint");
    }
    Report {
        swapped,
        cdna_strand,
        r1_profile: p1,
        r2_profile: p2,
        scores,
        verdict,
        structure,
        library,
        unchecked,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cDNA read is the file the barcode is not on.
    #[test]
    fn the_cdna_read_is_the_file_without_the_barcode() {
        let (r1, r2) = (Reads::with_capacity(1, 4), Reads::with_capacity(1, 4));
        let is = |got: Option<&Reads>, want: &Reads| got.is_some_and(|g| std::ptr::eq(g, want));
        assert!(is(cdna_read(&r1, Some(&r2), Mate::R1, false), &r2), "10x as labelled");
        assert!(is(cdna_read(&r1, Some(&r2), Mate::R2, false), &r1), "Parse as labelled");
        assert!(is(cdna_read(&r1, Some(&r2), Mate::R1, true), &r1), "10x swapped");
        assert!(is(cdna_read(&r1, Some(&r2), Mate::R2, true), &r2), "Parse swapped");
        assert!(cdna_read(&r1, None, Mate::R1, false).is_none(), "one file has no cDNA read");
    }
}
