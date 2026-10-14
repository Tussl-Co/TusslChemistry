// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Scores one candidate layout against the sampled reads: one pass tallies what each segment shows,
//! then small steps judge the tally (segment evidence, read length, UMI randomness, the tail) and
//! sum it into a support figure.

use crate::constants::{
    LINKER_BASES_PER_MISMATCH, LIST_TOP_SHARE_MAX, PRESENT_MIN, SEGMENT_HITS_MIN,
    TAIL_BASES_PER_MISMATCH, UMI_ENTROPY_MIN, VOUCHING_HITS,
};
use crate::detect::profile::{BASE_CLASSES, Profile, count_base, dominant, entropy, ratio};
use crate::hash::{FxHashMap, map_with_capacity};
use crate::input::fastq::Reads;
use crate::input::whitelist::{Hit, Whitelist};
use crate::layouts::{Layout, OFFSET_MAX, Segment, TAIL_MAX, Tail};
use crate::lists::{ListId, ListSet};

/// The loaded barcode lists, one slot per `ListId`, indexed by it, and the set of slots filled.
/// An optional list the user has not fetched (`ListId::optional`) is an empty slot; every other
/// slot is filled.
pub struct Lists {
    slots: [Option<Whitelist>; ListId::COUNT],
    loaded: ListSet,
}

impl Lists {
    pub(crate) fn new(loaded: impl IntoIterator<Item = (ListId, Whitelist)>) -> Self {
        let mut lists = Self { slots: std::array::from_fn(|_| None), loaded: ListSet::EMPTY };
        for (id, list) in loaded {
            let Some(slot) = lists.slots.get_mut(id.index()) else {
                unreachable!("list {id:?} has a slot: there is one per ListId")
            };
            assert!(slot.is_none(), "list {id:?} is loaded once");
            *slot = Some(list);
            lists.loaded = lists.loaded.with(id);
        }
        assert!(
            ListId::ALL.iter().all(|&id| id.optional() || lists.loaded.contains(id)),
            "every list that is not optional is loaded"
        );
        assert!(
            ListId::ALL.iter().all(|&id| lists.loaded.contains(id) == lists.slot(id).is_some()),
            "the loaded set names exactly the filled slots"
        );
        lists
    }

    fn slot(&self, id: ListId) -> Option<&Whitelist> {
        let Some(slot) = self.slots.get(id.index()) else {
            unreachable!("list {id:?} has a slot: there is one per ListId")
        };
        // Hot: checked build only.
        debug_assert!(slot.is_some() || id.optional(), "required list {id:?} is always loaded");
        slot.as_ref()
    }

    pub(crate) fn get(&self, id: ListId) -> &Whitelist {
        // Hot: checked build only.
        debug_assert!(self.loaded.contains(id), "only layouts whose lists are loaded are scored");
        let Some(whitelist) = self.slot(id) else {
            unreachable!("list {id:?} not loaded: only layouts whose lists are loaded are scored")
        };
        whitelist
    }

    /// Whether every list `layout` reads is loaded, so the layout can be scored.
    pub(crate) fn cover(&self, layout: &Layout) -> bool {
        let needed = layout.lists();
        let covered = needed.within(self.loaded);
        assert!(covered || needed.iter().any(ListId::optional), "{}: required lists", layout.name);
        assert!(!needed.is_empty() || covered, "{}: no lists, so covered", layout.name);
        covered
    }

    /// The loaded lists, in `ListId` order.
    pub(crate) fn loaded(&self) -> impl Iterator<Item = (ListId, &Whitelist)> {
        assert!(!self.loaded.is_empty(), "the required lists are loaded");
        assert!(self.loaded.len() <= ListId::COUNT, "at most every list is loaded");
        self.loaded.iter().filter_map(|id| self.slot(id).map(|list| (id, list)))
    }
}

fn mismatches(a: &[u8], b: &[u8]) -> usize {
    // Hot: checked build only. Equal lengths, or `zip` would silently compare a prefix.
    debug_assert_eq!(a.len(), b.len(), "compared windows have equal lengths");
    let count = a.iter().zip(b).filter(|(x, y)| x != y).count();
    debug_assert!(count <= a.len(), "at most one mismatch per base");
    count
}

/// Mismatches of a linker window against the closest of the layout's alternatives.
fn fixed_mismatches(window: &[u8], alternatives: &[&[u8]]) -> usize {
    // Hot: checked build only.
    debug_assert!(!alternatives.is_empty(), "a linker has a sequence");
    let best = alternatives.iter().map(|expected| mismatches(window, expected)).min();
    let Some(best) = best else { unreachable!("a linker has a sequence") };
    debug_assert!(best <= window.len(), "mismatches are at most the window");
    best
}

/// Linker tolerance: one mismatch per five bases (a 4-base linker must match exactly).
#[allow(clippy::integer_division)] // floor is the intent: one mismatch per whole five bases
fn fixed_allowed(len: usize) -> usize {
    // Hot: checked build only.
    debug_assert!(len > 0, "a linker has bases");
    let allowed = len / LINKER_BASES_PER_MISMATCH;
    debug_assert!(allowed < len, "a linker never tolerates mismatches at every base");
    allowed
}

/// Tail tolerance: one mismatch per four bases shown.
#[allow(clippy::integer_division)] // floor is the intent: one mismatch per whole four bases
fn tail_allowed(len: usize) -> usize {
    // Hot: checked build only.
    debug_assert!(len >= 2, "a tail is judged on two or more bases");
    let allowed = len / TAIL_BASES_PER_MISMATCH;
    debug_assert!(allowed < len, "a tail never tolerates mismatches at every base");
    allowed
}

/// 2-bit code for packing an unlisted barcode; N and other non-ACGT bases pack as A.
fn two_bit_code(base: u8) -> u32 {
    match base {
        b'C' => 1,
        b'G' => 2,
        b'T' => 3,
        _ => 0,
    }
}

/// Mean of a few per-segment rates (a layout has at most six segments).
#[allow(clippy::cast_precision_loss)] // the count is a handful of segments
fn mean(values: &[f64]) -> f64 {
    assert!(values.iter().all(|v| (0.0..=1.0).contains(v)), "every rate averaged is in 0..=1");
    let mean = values.iter().sum::<f64>() / values.len().max(1) as f64;
    assert!((0.0..=1.0).contains(&mean), "a mean of rates {mean} outside 0..=1");
    mean
}

/// Misfit of `sequence` against the layout's linkers and listed barcodes at `offset`; None when the
/// read is too short for any of them.
fn misfit(layout: &Layout, sequence: &[u8], offset: usize, lists: &Lists) -> Option<usize> {
    // Hot: checked build only.
    debug_assert!(
        layout.offsets.contains(&offset),
        "{} does not allow offset {offset}",
        layout.name
    );
    let mut position = offset;
    let mut cost = 0;
    let mut checked = false;
    for segment in layout.segments {
        let end = position + segment.len();
        match (segment, sequence.get(position..end)) {
            (Segment::Fixed(alternatives), Some(window)) => {
                cost += fixed_mismatches(window, alternatives);
                checked = true;
            }
            (Segment::Barcode { list: id, .. }, Some(window)) => {
                if lists.get(*id).hit(window) == Hit::Miss {
                    cost += 2;
                }
                checked = true;
            }
            _ => {}
        }
        position = end;
    }
    // Hot: checked build only. Each linker base costs at most 1 and each barcode (8+ bases) 2.
    debug_assert!(cost <= 2 * layout.structured_len(), "misfit {cost} past its bound");
    checked.then_some(cost)
}

fn best_offset(layout: &Layout, sequence: &[u8], lists: &Lists) -> usize {
    if layout.offsets.len() == 1 {
        return layout.first_offset();
    }
    let best = layout
        .offsets
        .iter()
        .filter_map(|&o| misfit(layout, sequence, o, lists).map(|c| (c, o)))
        .min()
        .map_or(layout.first_offset(), |(_, o)| o);
    assert!(layout.offsets.contains(&best), "chose offset {best}, not one {} allows", layout.name);
    assert!(best <= OFFSET_MAX, "offset {best} past OFFSET_MAX: the tally has no slot for it");
    best
}

pub enum SegmentEvidence {
    Listed {
        position: usize,
        list: ListId,
        rate: f64,
        eligible: usize,
        /// The share of the hits that the most common hit sequence holds (`LIST_TOP_SHARE_MAX`).
        top_share: f64,
    },
    Unlisted {
        position: usize,
        coverage: f64,
        set: usize,
        cover_min: f64,
        eligible: usize,
    },
    /// A linker: its first sequence, and how many alternatives the layout allows.
    Fixed {
        position: usize,
        sequence: &'static [u8],
        alternatives: usize,
        rate: f64,
        eligible: usize,
    },
}

pub struct Score {
    pub layout: &'static Layout,
    pub evidence: Vec<SegmentEvidence>,
    pub umi_entropy_min: Option<(usize, f64)>,
    /// Tail bases most reads show after the UMI (0 when the read ends there).
    pub tail_seen: usize,
    pub tail_rate: Option<f64>,
    /// Share of reads at each non-zero start offset, for layouts that allow several.
    pub offset_shares: Vec<(usize, f64)>,
    pub support: f64,
    pub contradictions: Vec<String>,
    /// Oddities that do not count against a layout a vendor list already vouches for.
    pub notes: Vec<String>,
}

/// One segment's counts from the pass over the reads.
#[derive(Clone, Default, PartialEq, Eq, Debug)]
struct SegmentTally {
    ok: usize,
    eligible: usize,
    /// How many reads carried each packed sequence: every unlisted barcode, and the hits of a
    /// listed one.
    kmers: FxHashMap<u32, usize>,
}

/// What one pass over the reads saw, per segment, before any judgement.
#[derive(PartialEq, Eq, Debug)]
struct Tally {
    segments: Vec<SegmentTally>,
    umi_counts: Vec<[u32; BASE_CLASSES]>,
    tail_ok: usize,
    tail_eligible: usize,
    /// Reads by how many tail bases they show, `0..=TAIL_MAX`.
    tail_lens: [usize; TAIL_MAX + 1],
    /// Reads by chosen start offset, `0..=OFFSET_MAX`.
    offset_counts: [usize; OFFSET_MAX + 1],
}

fn tally(layout: &Layout, sequences: &Reads, lists: &Lists) -> Tally {
    let mut t = tally_new(layout, sequences.len());
    for sequence in sequences.iter() {
        tally_read(layout, sequence, lists, &mut t);
    }
    for segment_tally in &t.segments {
        assert!(
            segment_tally.ok <= segment_tally.eligible,
            "segment hits {} past eligible {}",
            segment_tally.ok,
            segment_tally.eligible
        );
        assert!(
            segment_tally.eligible <= sequences.len(),
            "segment eligible {} past reads",
            segment_tally.eligible
        );
    }
    assert!(t.tail_ok <= t.tail_eligible, "tail hits past eligible");
    assert!(t.tail_eligible <= sequences.len(), "tail eligible past reads");
    let tail_reads = if layout.tail().is_some() { sequences.len() } else { 0 };
    assert_eq!(
        t.tail_lens.iter().sum::<usize>(),
        tail_reads,
        "every read's tail length counted once"
    );
    assert_eq!(
        t.offset_counts.iter().sum::<usize>(),
        sequences.len(),
        "every read placed at exactly one offset"
    );
    t
}

/// An empty tally for `layout`, sized once for `reads` reads.
fn tally_new(layout: &Layout, reads: usize) -> Tally {
    // Only the first UMI is tallied, so a second would be silently ignored.
    let umis = layout.segments.iter().filter(|s| matches!(s, Segment::Umi(_))).count();
    assert!(umis <= 1, "{} has {umis} UMIs; layouts_are_sound allows one", layout.name);
    let umi_len =
        layout.segments.iter().find_map(|s| if let Segment::Umi(l) = s { Some(*l) } else { None });
    let t = Tally {
        // Sized once: a barcode segment can see at most one new sequence per read.
        segments: layout
            .segments
            .iter()
            .map(|segment| SegmentTally {
                kmers: match segment {
                    Segment::Unlisted { .. } | Segment::Barcode { .. } => map_with_capacity(reads),
                    _ => FxHashMap::default(),
                },
                ..SegmentTally::default()
            })
            .collect(),
        umi_counts: vec![[0u32; BASE_CLASSES]; umi_len.unwrap_or(0)],
        tail_ok: 0,
        tail_eligible: 0,
        tail_lens: [0; TAIL_MAX + 1],
        offset_counts: [0; OFFSET_MAX + 1],
    };
    // `tally_read` zips segments with their tallies; unequal lengths would drop segments.
    assert_eq!(t.segments.len(), layout.segments.len(), "one tally per segment");
    t
}

/// Places one read at its best offset and counts what each segment shows.
fn tally_read(layout: &Layout, sequence: &[u8], lists: &Lists, t: &mut Tally) {
    let offset = best_offset(layout, sequence, lists);
    let Some(slot) = t.offset_counts.get_mut(offset) else {
        unreachable!("offset {offset} is one the layout allows, at most OFFSET_MAX (asserted)")
    };
    *slot += 1;
    let mut position = offset;
    for (segment, segment_tally) in layout.segments.iter().zip(t.segments.iter_mut()) {
        let end = position + segment.len();
        if let Segment::Tail(tail) = segment {
            let (shown, matched) = tally_tail(*tail, sequence, position);
            let Some(slot) = t.tail_lens.get_mut(shown) else {
                unreachable!("tally_tail shows at most TAIL_MAX bases (asserted)")
            };
            *slot += 1;
            if let Some(matched) = matched {
                t.tail_eligible += 1;
                t.tail_ok += usize::from(matched);
            }
        } else if let Some(window) = sequence.get(position..end) {
            tally_segment(*segment, window, segment_tally, &mut t.umi_counts, lists);
        }
        position = end;
    }
    // Hot: checked build only.
    debug_assert_eq!(position, offset + layout.structured_len(), "segments walked end to end");
}

/// Tail bases the read shows from `position` (at most `TAIL_MAX`), and whether they match `tail`
/// within tolerance; None when fewer than two show, too few to judge.
fn tally_tail(tail: Tail, sequence: &[u8], position: usize) -> (usize, Option<bool>) {
    let shown = sequence.len().saturating_sub(position).min(TAIL_MAX);
    // Hot: checked build only.
    debug_assert!(shown <= TAIL_MAX, "a tail window is at most TAIL_MAX bases");
    if shown < 2 {
        return (shown, None);
    }
    let Some(window) = sequence.get(position..position + shown) else {
        unreachable!("shown = min(remaining, TAIL_MAX) keeps the tail window inside the read")
    };
    debug_assert_eq!(window.len(), shown, "the window holds every shown base");
    (shown, Some(mismatches(window, tail.expected(shown)) <= tail_allowed(shown)))
}

/// Counts one segment's window: a linker or listed-barcode hit, an unlisted barcode's sequence,
/// or the UMI's bases per position.
fn tally_segment(
    segment: Segment,
    window: &[u8],
    segment_tally: &mut SegmentTally,
    umi_counts: &mut [[u32; BASE_CLASSES]],
    lists: &Lists,
) {
    // Hot: checked build only.
    debug_assert_eq!(window.len(), segment.len(), "the window spans exactly the segment");
    match segment {
        Segment::Fixed(alternatives) => {
            segment_tally.eligible += 1;
            segment_tally.ok +=
                usize::from(fixed_mismatches(window, alternatives) <= fixed_allowed(window.len()));
        }
        Segment::Barcode { list: id, .. } => {
            segment_tally.eligible += 1;
            if lists.get(id).hit(window) != Hit::Miss {
                segment_tally.ok += 1;
                *segment_tally.kmers.entry(pack(window)).or_insert(0) += 1;
            }
        }
        Segment::Unlisted { .. } => {
            segment_tally.eligible += 1;
            *segment_tally.kmers.entry(pack(window)).or_insert(0) += 1;
        }
        Segment::Umi(_) => {
            debug_assert_eq!(umi_counts.len(), window.len(), "one UMI column per UMI base");
            for (slot, &base) in umi_counts.iter_mut().zip(window) {
                count_base(slot, base);
            }
        }
        Segment::Tail(_) => unreachable!("tally_read handles tails before windows"),
    }
    debug_assert!(
        segment_tally.ok <= segment_tally.eligible,
        "a segment hits at most the reads it saw"
    );
}

/// Whether some listed segment's hits are spread over many sequences, as a real library's are: hit
/// in at least `SEGMENT_HITS_MIN` of reads, with no one sequence holding `LIST_TOP_SHARE_MAX` of
/// the hits. A constant read that happens to sit near a list entry is not.
pub fn lists_spread(evidence: &[SegmentEvidence]) -> bool {
    assert!(
        evidence.iter().all(|e| match e {
            SegmentEvidence::Listed { top_share, .. } => (0.0..=1.0).contains(top_share),
            _ => true,
        }),
        "a top share is a share"
    );
    let spread = evidence.iter().any(|e| {
        matches!(e, SegmentEvidence::Listed { rate, top_share, .. }
            if *rate >= SEGMENT_HITS_MIN && *top_share < LIST_TOP_SHARE_MAX)
    });
    assert!(
        !spread || evidence.iter().any(|e| matches!(e, SegmentEvidence::Listed { .. })),
        "only listed evidence is spread"
    );
    spread
}

/// Whether a layout's lists vouch for it: one hit in `VOUCHING_HITS` of reads, and the hits spread
/// (`lists_spread`).
pub fn lists_vouch(evidence: &[SegmentEvidence]) -> bool {
    let high = evidence
        .iter()
        .any(|e| matches!(e, SegmentEvidence::Listed { rate, .. } if *rate >= VOUCHING_HITS));
    let vouch = high && lists_spread(evidence);
    assert!(!vouch || high, "vouching needs a high hit rate");
    assert!(!vouch || lists_spread(evidence), "vouching needs spread hits");
    vouch
}

/// A window of at most 16 bases as two bits a base (`LEN_MAX`).
fn pack(window: &[u8]) -> u32 {
    // Hot: checked build only.
    debug_assert!(window.len() <= 16, "a barcode window packs into 32 bits");
    let key = window.iter().fold(0u32, |k, &b| (k << 2) | two_bit_code(b));
    debug_assert!(window.len() == 16 || key >> (2 * window.len()) == 0, "the key holds the window");
    key
}

/// The share of a listed segment's hits that its most common hit sequence holds; 0 with no hits.
fn top_share(segment_tally: &SegmentTally) -> f64 {
    let top = segment_tally.kmers.values().copied().max().unwrap_or(0);
    assert!(top <= segment_tally.ok, "one sequence holds at most every hit");
    assert!(
        segment_tally.kmers.values().sum::<usize>() == segment_tally.ok,
        "every hit counted once by sequence"
    );
    ratio(top, segment_tally.ok)
}

/// Each segment's evidence, and a contradiction for every segment the reads do not support.
fn segment_evidence(
    layout: &Layout,
    t: &Tally,
    contradictions: &mut Vec<String>,
) -> Vec<SegmentEvidence> {
    assert_eq!(t.segments.len(), layout.segments.len(), "one tally per segment (tally_new)");
    // Every read placed at one offset (tally): the count every segment's rate is taken over.
    let reads = t.offset_counts.iter().sum::<usize>();
    let contradictions_before = contradictions.len();
    let mut evidence = Vec::new();
    let mut position = 0;
    for (segment, segment_tally) in layout.segments.iter().zip(&t.segments) {
        evidence.extend(segment_evidence_judge(
            segment,
            segment_tally,
            reads,
            position,
            contradictions,
        ));
        position += segment.len();
    }
    // Each judged segment adds at most one contradiction.
    assert!(
        contradictions.len() - contradictions_before <= evidence.len(),
        "more contradictions than judged segments"
    );
    assert_eq!(position, layout.structured_len(), "every segment walked");
    evidence
}

/// The evidence one segment at `position` gives, with a contradiction when the reads do not bear
/// it out; None for segments judged elsewhere (UMI, tail).
///
/// Rates are taken over all `reads`, not only those long enough to reach the segment: a read too
/// short for it does not show it. Otherwise a file mixing short reads of another kit with long
/// reads of this one scores the late segments on this kit's reads alone, and the mix reaches
/// HIGH (chaos, `PIPseq` and ScaleBio mixed with 28-base 10x reads, 9 Oct 2026). A contradiction
/// still needs some read to reach the segment; a read too short for the whole layout is
/// contradicted once, by its length (`score`).
fn segment_evidence_judge(
    segment: &Segment,
    segment_tally: &SegmentTally,
    reads: usize,
    position: usize,
    contradictions: &mut Vec<String>,
) -> Option<SegmentEvidence> {
    assert!(segment_tally.ok <= segment_tally.eligible, "a segment hits at most the reads it saw");
    assert!(segment_tally.eligible <= reads, "a segment sees at most every read");
    let contradictions_before = contradictions.len();
    let eligible = segment_tally.eligible;
    let rate = ratio(segment_tally.ok, reads);
    let evidence = match segment {
        Segment::Barcode { list: id, .. } => {
            if eligible > 0 && rate < SEGMENT_HITS_MIN {
                let label = id.label();
                contradictions
                    .push(format!("barcode at {position}: {:.0}% in {label}", rate * 100.0));
            }
            let top_share = top_share(segment_tally);
            SegmentEvidence::Listed { position, list: *id, rate, eligible, top_share }
        }
        Segment::Unlisted { set, cover_pct_min, .. } => {
            let coverage = segment_evidence_coverage(segment_tally, *set, reads);
            let cover_min = f64::from(*cover_pct_min) / 100.0;
            if eligible > 0 && coverage < cover_min {
                contradictions.push(format!(
                    "position {position} looks random, not a barcode (top {set} sequences \
                     cover {:.0}%, need {cover_pct_min}%)",
                    coverage * 100.0
                ));
            }
            SegmentEvidence::Unlisted { position, coverage, set: *set, cover_min, eligible }
        }
        Segment::Fixed(alternatives) => {
            let Some(&sequence) = alternatives.first() else {
                unreachable!("checked at compile time: a linker has a sequence")
            };
            if eligible > 0 && rate < PRESENT_MIN {
                let others = match alternatives.len() {
                    1 => String::new(),
                    n => format!(" (or one of {} alternatives)", n - 1),
                };
                contradictions.push(format!(
                    "linker {}{others} expected at {position}, found in {:.0}% of reads",
                    String::from_utf8_lossy(sequence),
                    rate * 100.0
                ));
            }
            let alternatives = alternatives.len();
            SegmentEvidence::Fixed { position, sequence, alternatives, rate, eligible }
        }
        _ => return None,
    };
    assert!(contradictions.len() - contradictions_before <= 1, "one segment, one contradiction");
    Some(evidence)
}

/// Share of all `reads` that the `set` commonest sequences of an unlisted segment cover.
fn segment_evidence_coverage(segment_tally: &SegmentTally, set: usize, reads: usize) -> f64 {
    assert!(set > 0, "an unlisted segment names a non-empty top set");
    assert!(segment_tally.eligible <= reads, "a segment sees at most every read");
    // Paired with `tally_segment`: every eligible read added exactly one sequence.
    assert_eq!(
        segment_tally.kmers.values().sum::<usize>(),
        segment_tally.eligible,
        "one sequence per eligible read"
    );
    let mut counts: Vec<usize> = segment_tally.kmers.values().copied().collect();
    counts.sort_unstable_by(|a, b| b.cmp(a));
    let top: usize = counts.iter().take(set).sum();
    ratio(top, reads)
}

/// The least random UMI position, with a contradiction when it is not random enough to be a UMI.
fn umi_check(layout: &Layout, t: &Tally, contradictions: &mut Vec<String>) -> Option<(usize, f64)> {
    let umi_len =
        layout.segments.iter().find_map(|s| if let Segment::Umi(l) = s { Some(*l) } else { None });
    assert_eq!(t.umi_counts.len(), umi_len.unwrap_or(0), "one UMI column per UMI base (tally_new)");
    let umi_start = layout
        .segments
        .iter()
        .take_while(|s| !matches!(s, Segment::Umi(_)))
        .map(|s| s.len())
        .sum::<usize>();
    let least = t
        .umi_counts
        .iter()
        .enumerate()
        .filter_map(|(j, c)| entropy(c).map(|h| (umi_start + j, h, c)))
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .map(|(p, h, c)| {
            if h < UMI_ENTROPY_MIN {
                let (base, frac) = dominant(c);
                contradictions.push(format!(
                    "UMI position {p} is not random ({h:.2} bits, {:.0}% {})",
                    frac * 100.0,
                    base as char
                ));
            }
            (p, h)
        });
    if let Some((p, _)) = least {
        assert!(p >= umi_start, "the least random UMI position {p} precedes the UMI");
        assert!(p < umi_start + t.umi_counts.len(), "the least random UMI position {p} follows it");
    }
    least
}

/// Tail bases most reads show, the share showing the expected tail, and a note or contradiction
/// when it is missing.
fn tail_check(
    layout: &Layout,
    t: &Tally,
    evidence: &[SegmentEvidence],
    contradictions: &mut Vec<String>,
) -> (usize, Option<f64>, Vec<String>) {
    assert!(t.tail_ok <= t.tail_eligible, "tail hits past eligible (tally)");
    // Only lengths some read showed; on a tie the longer length.
    let tail_seen = (0..)
        .zip(t.tail_lens)
        .filter(|&(_, n)| n > 0)
        .max_by_key(|&(k, n)| (n, k))
        .map_or(0, |(k, _)| k);
    let tail_rate = (t.tail_eligible > 0).then(|| ratio(t.tail_ok, t.tail_eligible));
    let mut notes = Vec::new();
    if let (Some(rate), Some(tail)) = (tail_rate, layout.tail())
        && tail_seen >= 2
        && rate < PRESENT_MIN
    {
        let msg = format!(
            "expected {} after the UMI, seen in {:.0}% of reads",
            tail.label(),
            rate * 100.0
        );
        // Trimmed or BAM-rebuilt reads lose their tail but cannot fake a whitelist, so for a
        // list-vouched layout a missing tail is a note; for an unlisted one it is a contradiction.
        if lists_vouch(evidence) { notes.push(msg) } else { contradictions.push(msg) }
    }
    assert!(tail_seen <= TAIL_MAX, "tail bases seen {tail_seen} past TAIL_MAX");
    assert!(notes.len() <= 1, "the missing tail is the only note");
    (tail_seen, tail_rate, notes)
}

/// Support: mean over every checkable segment (listed barcodes, linkers, unlisted-barcode
/// concentration), so each layout answers for its whole structure, not only its easiest parts.
fn support_of(evidence: &[SegmentEvidence]) -> f64 {
    let rates: Vec<f64> = evidence
        .iter()
        .map(|e| match e {
            SegmentEvidence::Listed { rate, eligible, .. }
            | SegmentEvidence::Fixed { rate, eligible, .. } => {
                if *eligible > 0 {
                    *rate
                } else {
                    0.0
                }
            }
            // Concentration is a pass/fail signal, scaled so meeting the bar scores like a hit.
            SegmentEvidence::Unlisted { coverage, cover_min, eligible, .. } => {
                assert!(*cover_min > 0.0, "a zero concentration bar would divide by zero");
                if *eligible > 0 { (coverage / cover_min).min(1.0) } else { 0.0 }
            }
        })
        .collect();
    let support = mean(&rates);
    assert!((0.0..=1.0).contains(&support), "support {support} outside 0..=1");
    support
}

fn offset_shares_of(layout: &Layout, t: &Tally, reads: usize) -> Vec<(usize, f64)> {
    assert_eq!(t.offset_counts.iter().sum::<usize>(), reads, "every read placed once (tally)");
    let mut offset_shares: Vec<(usize, f64)> = if layout.offsets.len() > 1 {
        (0..)
            .zip(t.offset_counts)
            .filter(|&(_, n)| n > 0)
            .map(|(o, n)| (o, ratio(n, reads)))
            .collect()
    } else {
        Vec::new()
    };
    offset_shares.sort_by_key(|&(o, _)| o);
    assert!(
        offset_shares.iter().all(|(o, _)| layout.offsets.contains(o)),
        "reads sit only at offsets {} allows",
        layout.name
    );
    offset_shares
}

pub fn score(
    layout: &'static Layout,
    sequences: &Reads,
    profile: &Profile,
    lists: &Lists,
) -> Score {
    assert_eq!(profile.reads, sequences.len(), "the profile is of these reads");
    let t = tally(layout, sequences, lists);
    // Debug only, because it doubles the scoring work: a tally is a pure function of the set of
    // reads, so the same reads in reverse order must give identical counts. Catches any hidden
    // dependence on read order or on hash-map iteration order.
    #[cfg(debug_assertions)]
    assert_eq!(
        tally(layout, &sequences.reversed(), lists),
        t,
        "{}: tally depends on read order",
        layout.name
    );
    let mut contradictions = Vec::new();
    let evidence = segment_evidence(layout, &t, &mut contradictions);

    let needed = layout.first_offset() + layout.structured_len();
    if profile.len_mode < needed {
        contradictions.push(format!(
            "{:?} is {} bases, too short for this layout (needs {needed})",
            layout.read, profile.len_mode
        ));
    }

    let umi_entropy_min = umi_check(layout, &t, &mut contradictions);
    let (tail_seen, tail_rate, notes) = tail_check(layout, &t, &evidence, &mut contradictions);
    let support = support_of(&evidence);
    let offset_shares = offset_shares_of(layout, &t, sequences.len());
    let score = Score {
        layout,
        evidence,
        umi_entropy_min,
        tail_seen,
        tail_rate,
        offset_shares,
        support,
        contradictions,
        notes,
    };
    // Paired with `tail_check`: a missing tail is only a note for a list-vouched layout.
    if !score.notes.is_empty() {
        assert!(vouched(&score), "{}: a note on a layout no list vouches for", layout.name);
    }
    score
}

pub fn vouched(score: &Score) -> bool {
    assert!(
        score.evidence.len() <= score.layout.segments.len(),
        "at most one evidence per segment"
    );
    let vouched = lists_vouch(&score.evidence);
    if vouched {
        assert!(!score.layout.lists().is_empty(), "only a layout with a barcode list is vouched");
    }
    vouched
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::detect::profile::Profile;
    use crate::input::whitelist::Whitelist;
    use crate::layouts::{LAYOUTS, LayoutSet};
    use crate::signals::kmers::CdnaHits;
    use crate::signals::strand::StrandEvidence;

    /// Deterministic pseudo-random reads: a 10x-like 28 bp barcode read, a cDNA-like 90 bp read.
    fn reads(n: usize, len: usize, seed: u64) -> Reads {
        let mut state = seed;
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

    /// Two barcodes per list, cut to each list's own length, so every list in the registry is
    /// loaded with barcodes it can hold.
    fn lists() -> Lists {
        let patterns: [&[u8]; 2] = [b"ACGTACGTACGTACGT", b"TTTTCCCCGGGGAAAA"];
        Lists::new(ListId::ALL.iter().map(|&id| {
            let length = id.barcode_length();
            let barcodes: Vec<&[u8]> = patterns.iter().map(|p| &p[..length]).collect();
            (id, Whitelist::from_barcodes(&barcodes))
        }))
    }

    /// The hot loops must not allocate per read: allocation counts on 1,000 and 10,000 reads match.
    #[test]
    fn hot_loops_allocate_independently_of_read_count() {
        let lists = lists();
        // Warm-up: the reference k-mer table is built once, on first use, and must not count here.
        let one = reads(1, 90, 1);
        std::hint::black_box(StrandEvidence::of(&one, &CdnaHits::of(&one)));
        for (len, label) in [(28, "barcode read"), (90, "cDNA read")] {
            let (small, large) = (reads(1_000, len, 7), reads(10_000, len, 7));
            for layout in LAYOUTS {
                let count = |r: &Reads| {
                    allocation_counter::measure(|| drop(tally(layout, r, &lists))).count_total
                };
                assert_eq!(
                    count(&small),
                    count(&large),
                    "{}: tally allocates per read on a {label}",
                    layout.name
                );
            }
            let strand = |r: &Reads| {
                allocation_counter::measure(|| {
                    std::hint::black_box(StrandEvidence::of(r, &CdnaHits::of(r)));
                })
                .count_total
            };
            assert_eq!(
                strand(&small),
                strand(&large),
                "strand evidence allocates per read on a {label}"
            );
            let profile =
                |r: &Reads| allocation_counter::measure(|| drop(Profile::of(r))).count_total;
            assert_eq!(profile(&small), profile(&large), "profile allocates per read on a {label}");
        }
    }

    /// Deciding which layouts can be scored, and comparing their lists and linkers, allocates
    /// nothing: the sets are fixed-size and the linkers are iterated in place.
    #[test]
    fn layout_sets_allocate_nothing() {
        let lists = lists();
        let allocations = allocation_counter::measure(|| {
            let mut unchecked = LayoutSet::EMPTY;
            for (index, layout) in LAYOUTS.iter().enumerate() {
                if !lists.cover(layout) {
                    unchecked = unchecked.with(index);
                }
                for other in LAYOUTS {
                    let same =
                        layout.lists() == other.lists() && layout.linkers().eq(other.linkers());
                    std::hint::black_box(same);
                }
            }
            std::hint::black_box(unchecked.iter().count());
        })
        .count_total;
        assert_eq!(allocations, 0, "layout and list sets allocate");
    }

    /// Scoring a layout returns everything it allocated: no leak across calls.
    #[test]
    fn scoring_does_not_leak() {
        let (lists, r) = (lists(), reads(2_000, 28, 11));
        let profile = Profile::of(&r);
        for layout in LAYOUTS {
            let info = allocation_counter::measure(|| drop(score(layout, &r, &profile, &lists)));
            assert_eq!(
                info.count_current, 0,
                "{} leaked {} allocations",
                layout.name, info.count_current
            );
        }
    }
}
