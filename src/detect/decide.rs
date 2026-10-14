// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Turns the ranked layout scores into a verdict: a match with a confidence, an honest ambiguity,
//! or no known layout.

use crate::constants::{
    PRESENT_MIN, RIVAL_SUPPORT_SHARE, SUPPORT_HIGH, SUPPORT_MIN, SUPPORT_TIE_TOLERANCE,
    TAIL_SEEN_BASES_MIN, UNLISTED_LINKER_HIGH,
};
use crate::detect::score::{Score, SegmentEvidence, lists_spread, vouched};
use crate::layouts::Tail;
use crate::signals::strand::{Strand, StrandEvidence};

/// How firmly the reads support a matched layout.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Confidence {
    /// At least 70% support, nothing contradicting it, no rival holding half as much, and evidence
    /// from a barcode list or a linker.
    High,
    /// A clear winner short of HIGH.
    Medium,
    /// A layout a vendor list vouches for (80% of barcodes listed) while something else is off.
    Low,
}

impl Confidence {
    /// The confidence as the text report and `--emit` print it.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Confidence::High => "HIGH",
            Confidence::Medium => "MEDIUM",
            Confidence::Low => "LOW",
        }
    }
}

pub enum Verdict {
    Match {
        best: usize,
        confidence: Confidence,
    },
    /// Equally supported layouts that would be processed differently.
    Ambiguous {
        candidates: Vec<usize>,
        reason: String,
    },
    NoMatch,
}

pub fn decide(scores: &[Score], strand: Option<&StrandEvidence>) -> Verdict {
    let verdict = decide_unchecked(scores, strand);
    // Postconditions: every index names a score, a LOW match stands on a vouching list, and an
    // ambiguity has at least two candidates.
    match &verdict {
        Verdict::Match { best, confidence } => {
            let Some(chosen) = scores.get(*best) else {
                unreachable!("verdict index {best} past {} scores", scores.len())
            };
            assert!(
                *confidence != Confidence::Low || vouched(chosen),
                "a LOW match must be list-vouched"
            );
        }
        Verdict::Ambiguous { candidates, .. } => {
            assert!(
                candidates.len() >= 2,
                "an ambiguity needs two candidates, got {}",
                candidates.len()
            );
            assert!(
                candidates.iter().all(|&i| i < scores.len()),
                "ambiguity index past {} scores",
                scores.len()
            );
        }
        Verdict::NoMatch => {}
    }
    verdict
}

fn decide_unchecked(scores: &[Score], strand: Option<&StrandEvidence>) -> Verdict {
    // Preconditions, paired with `detect::score_all`'s sort: clean layouts come first, and every
    // support is a share. The rules below read `scores[0]` as the best clean layout.
    assert!(
        scores
            .iter()
            .skip_while(|s| s.contradictions.is_empty())
            .all(|s| !s.contradictions.is_empty()),
        "layouts the reads do not contradict sort before those they do"
    );
    assert!(scores.iter().all(|s| (0.0..=1.0).contains(&s.support)), "every support is in 0..=1");
    let Some((top, rest)) = scores.split_first() else { return Verdict::NoMatch };
    if !top.contradictions.is_empty() {
        // Nothing fits cleanly. Claim the best layout a vendor list vouches for, at LOW, with its
        // contradictions shown; an unlisted structure the reads contradict is too easy to fit by
        // accident to claim at all.
        return match scores.iter().position(vouched) {
            Some(best) => Verdict::Match { best, confidence: Confidence::Low },
            None => Verdict::NoMatch,
        };
    }
    if top.support < SUPPORT_MIN {
        return Verdict::NoMatch;
    }

    // Layouts the reads cannot separate: same lists, same linkers, equal support and standing.
    let tied: Vec<(usize, &Score)> = scores
        .iter()
        .enumerate()
        .filter(|(_, s)| {
            s.layout.lists() == top.layout.lists()
                && s.layout.linkers().eq(top.layout.linkers())
                && (s.support - top.support).abs() < SUPPORT_TIE_TOLERANCE
                && s.contradictions.len() == top.contradictions.len()
        })
        .collect();
    assert!(
        tied.first().is_some_and(|&(index, _)| index == 0),
        "the top layout ties with itself, first"
    );
    if tied.len() > 1 {
        return resolve_tie(top, &tied, strand);
    }
    Verdict::Match { best: 0, confidence: confidence_of(top, rest) }
}

/// Separates layouts the barcode evidence cannot: a tail actually seen first, then the cDNA read's
/// strand; otherwise an honest ambiguity with the reason.
fn resolve_tie(top: &Score, tied: &[(usize, &Score)], strand: Option<&StrandEvidence>) -> Verdict {
    assert!(tied.len() >= 2, "a tie needs two layouts, got {}", tied.len());
    assert!(
        tied.iter().all(|(_, s)| s.layout.lists() == top.layout.lists()),
        "tied layouts share the top layout's barcode lists"
    );
    let tails_differ = tied.iter().any(|(_, s)| s.layout.tail() != top.layout.tail());
    // Positive tail evidence decides first: a missing tail is only a note for list-vouched
    // layouts, but a tail actually seen (4+ bases, most reads) still tells 3' from 5'.
    if tails_differ {
        let seen: Vec<usize> = tied
            .iter()
            .filter(|(_, s)| {
                s.tail_seen >= TAIL_SEEN_BASES_MIN && s.tail_rate.is_some_and(|r| r >= PRESENT_MIN)
            })
            .map(|&(i, _)| i)
            .collect();
        if let [best] = seen[..] {
            assert!(tied.iter().any(|&(i, _)| i == best), "the tail picks a tied layout");
            return Verdict::Match { best, confidence: capped(Confidence::High, top) };
        }
    }
    // The cDNA read's strand separates 3' (sense) from 5' (antisense) without a reference.
    if tails_differ && let Some((strand, strong)) = strand.and_then(StrandEvidence::call) {
        let want = match strand {
            Strand::Sense => Tail::PolyT,
            Strand::Antisense => Tail::Tso,
        };
        let fits: Vec<usize> =
            tied.iter().filter(|(_, s)| s.layout.tail() == Some(want)).map(|&(i, _)| i).collect();
        if let [best] = fits[..] {
            assert!(tied.iter().any(|&(i, _)| i == best), "the strand picks a tied layout");
            let confidence = if strong { Confidence::High } else { Confidence::Medium };
            return Verdict::Match { best, confidence: capped(confidence, top) };
        }
    }
    let reason = if tails_differ {
        let seen = if top.tail_seen == 0 {
            "the read ends at the UMI, so polyT (3') and TSO (5') cannot be seen".to_string()
        } else {
            format!(
                "the read shows only {} base(s) after the UMI, the same for polyT and TSO",
                top.tail_seen
            )
        };
        format!(
            "{seen}, and the cDNA read's polyA/polyT runs do not settle the strand; 3' and 5' \
                count opposite strands"
        )
    } else {
        "the reads fit both equally".to_string()
    };
    Verdict::Ambiguous { candidates: tied.iter().map(|&(i, _)| i).collect(), reason }
}

/// A tie-break can raise confidence only as far as the layout's own support allows.
fn capped(confidence: Confidence, top: &Score) -> Confidence {
    assert!((0.0..=1.0).contains(&top.support), "support {} outside 0..=1", top.support);
    let lowered = confidence == Confidence::High && top.support < SUPPORT_HIGH;
    let result = if lowered { Confidence::Medium } else { confidence };
    if result == Confidence::High {
        assert!(top.support >= SUPPORT_HIGH, "HIGH stands on support of at least SUPPORT_HIGH");
    }
    result
}

/// Confidence of a clean, untied top layout against the rest.
fn confidence_of(top: &Score, rest: &[Score]) -> Confidence {
    assert!(top.contradictions.is_empty(), "only a clean top layout is rated here");
    assert!(top.support >= SUPPORT_MIN, "a rated layout clears SUPPORT_MIN");
    let others_hold = rest
        .iter()
        .any(|s| s.contradictions.is_empty() && s.support > RIVAL_SUPPORT_SHARE * top.support);
    // HIGH needs list or linker evidence; concentration and UMI randomness alone cap at MEDIUM.
    // A list counts only when its hits are spread over many sequences (`lists_spread`).
    let listed = lists_spread(&top.evidence);
    let anchored =
        listed || top.evidence.iter().any(|e| matches!(e, SegmentEvidence::Fixed { .. }));
    // Without a list, the linkers are the only direct count of reads from this kit: a mixture
    // carries them in about its share of reads (UNLISTED_LINKER_HIGH).
    let pure = listed
        || top.evidence.iter().all(|e| match e {
            SegmentEvidence::Fixed { rate, .. } => *rate >= UNLISTED_LINKER_HIGH,
            _ => true,
        });
    if top.support >= SUPPORT_HIGH && !others_hold && anchored && pure {
        Confidence::High
    } else {
        Confidence::Medium
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layouts::LAYOUTS;

    /// ScaleBio RNA's evidence as a half-and-half mix shows it: both barcode parts just over their
    /// 50% concentration bar (full marks in support), and the linker in `linker` of reads.
    fn scalebio_with_linker(linker: f64) -> Score {
        assert!((0.0..=1.0).contains(&linker), "a linker rate is a share: {linker}");
        let layout = LAYOUTS.iter().find(|l| l.name == "ScaleBio RNA");
        let Some(layout) = layout else { unreachable!("ScaleBio RNA is a layout") };
        let unlisted = |position, coverage| SegmentEvidence::Unlisted {
            position,
            coverage,
            set: 384,
            cover_min: 0.5,
            eligible: 1000,
        };
        let evidence = vec![
            unlisted(0, 0.53),
            SegmentEvidence::Fixed {
                position: 8,
                sequence: b"TCAGAGC",
                alternatives: 1,
                rate: linker,
                eligible: 1000,
            },
            unlisted(23, 0.50),
        ];
        let support = (1.0 + linker + 1.0) / 3.0;
        assert!(support >= SUPPORT_HIGH, "the mix clears SUPPORT_HIGH on support alone");
        Score {
            layout,
            evidence,
            umi_entropy_min: None,
            tail_seen: 0,
            tail_rate: None,
            offset_shares: Vec::new(),
            support,
            contradictions: Vec::new(),
            notes: Vec::new(),
        }
    }

    /// On a 50/50 mix, a layout with no barcode list passes its barcode-concentration bar in full,
    /// so only the linker shows the mix.
    #[test]
    fn a_listless_layout_needs_its_linkers_in_most_reads_for_high() {
        let mix = confidence_of(&scalebio_with_linker(0.50), &[]);
        assert_eq!(mix, Confidence::Medium, "a 50/50 mix is not HIGH");
        let pure = confidence_of(&scalebio_with_linker(0.98), &[]);
        assert_eq!(pure, Confidence::High, "a pure library still is");
    }
}
