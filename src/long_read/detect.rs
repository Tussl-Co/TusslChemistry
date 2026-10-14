// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Which 10x kit made a long-read single-cell library. Every 10x molecule opens with the Read 1
//! adapter, then a 16-base barcode, then the UMI, then polyT (3' kits) or the template-switch oligo
//! (5' kits). Per read: find the adapter at either end (either orientation), read the barcode after
//! it against each 10x list, and find what follows the UMI. The kit is the one whose list, tail and
//! UMI length the reads agree on. Candidates come from the layout table, so short and long reads
//! can never disagree about what a kit is.

use crate::constants::{
    LONG_ADAPTER_EDITS_MAX, LONG_ADAPTER_SHARE_MIN, LONG_ADAPTER_WINDOW, LONG_BARCODE_SHIFT_MAX,
    LONG_LIST_SEPARATION, LONG_LIST_SHARE_HIGH, LONG_LIST_SHARE_MIN, LONG_POLY_T_RUN,
    LONG_READS_EXAMINED, LONG_TAIL_SEPARATION, LONG_TAIL_SHARE_MIN, LONG_TAIL_WINDOW,
    LONG_TSO_EDITS_MAX, LONG_UMI_TOLERANCE,
};
use crate::detect::decide::Confidence;
use crate::detect::profile::ratio;
use crate::detect::score::Lists;
use crate::input::fastq::Reads;
use crate::input::whitelist::Hit;
use crate::layouts::{LAYOUTS, Layout, Mate, Segment, TSO, Tail};
use crate::lists::ListId;
use crate::long_read::search::Pattern;
use crate::signals::kmers::complement;

/// 10x's partial Read 1 primer, which opens every Chromium gene-expression molecule.
pub const READ1_ADAPTER: &[u8] = b"CTACACGACGCTCTTCCGATCT";
const BARCODE_BASES: usize = 16;
/// The shortest read that can carry a molecule's front end: the adapter, a barcode and the longest
/// UMI. Shorter reads (fragments, the junk some runs open with) cannot show the structure, so they
/// are no evidence against it and are not counted.
pub const USABLE_BASES_MIN: usize = READ1_ADAPTER.len() + BARCODE_BASES + 12;
/// The 10x lists a long read's barcode is tested against, in a fixed order.
pub const TENX_LISTS: [ListId; 4] =
    [ListId::TenxV2, ListId::TenxV3, ListId::TenxV4ThreePrime, ListId::TenxV3FivePrime];
/// UMI lengths measured, 0 up to this; a longer gap means the tail was not where a kit puts it.
const UMI_MEASURED_MAX: usize = 24;

const _: () = assert!(LONG_TAIL_WINDOW > 12 + 13, "the window covers the longest UMI and the TSO");
const _: () = assert!(UMI_MEASURED_MAX < LONG_TAIL_WINDOW);

/// What the long reads show, before any judgement.
#[derive(Debug, PartialEq, Eq)]
pub struct LongEvidence {
    pub reads: usize,
    /// Reads at least `USABLE_BASES_MIN` long, the ones searched.
    pub usable: usize,
    /// Reads where the Read 1 adapter was found, and of those, how many at the read's start.
    pub adapter: usize,
    pub forward: usize,
    /// Per `TENX_LISTS` entry: adapter-bearing reads whose barcode is on that list exactly.
    pub list_hits: [usize; 4],
    /// Per `TENX_LISTS` entry: reads whose barcode is on that list and on no other. The lists share
    /// barcodes (737K and v3 77,142; 3' v4 and 5' v3 253,420), and a shared barcode is no evidence
    /// between them, so lists are told apart on these.
    pub list_exclusive: [usize; 4],
    /// Adapter-bearing reads whose UMI is followed by polyT, and by the TSO.
    pub poly_t: usize,
    pub tso: usize,
    /// Reads by measured UMI length (barcode end to tail start), 0..=`UMI_MEASURED_MAX`.
    pub umi_lengths: [usize; UMI_MEASURED_MAX + 1],
}

/// A 10x kit as a long read shows it: its list, UMI length and tail, from its layout.
#[derive(Clone, Copy)]
pub struct Kit {
    pub layout: &'static Layout,
    pub list: ListId,
    pub umi: usize,
    pub tail: Tail,
}

/// Every layout a long read can show: barcode read R1, 16-base listed barcode, UMI, tail.
pub fn kits() -> Vec<Kit> {
    let kits: Vec<Kit> = LAYOUTS
        .iter()
        .filter(|layout| matches!(layout.read, Mate::R1) && layout.long_read_kit.is_some())
        .filter_map(|layout| match layout.segments {
            [
                Segment::Barcode { list, len: BARCODE_BASES },
                Segment::Umi(umi),
                Segment::Tail(tail),
            ] => Some(Kit { layout, list: *list, umi: *umi, tail: *tail }),
            _ => None,
        })
        .collect();
    assert!(!kits.is_empty(), "the layout table holds 10x kits");
    assert!(kits.iter().all(|k| TENX_LISTS.contains(&k.list)), "every kit's list is tested");
    kits
}

/// Tallies what the first `LONG_READS_EXAMINED` reads show.
pub fn evidence(reads: &Reads, lists: &Lists) -> LongEvidence {
    let adapter = Pattern::new(READ1_ADAPTER);
    let tso = Pattern::new(TSO);
    let mut evidence = LongEvidence {
        reads: 0,
        usable: 0,
        adapter: 0,
        forward: 0,
        list_hits: [0; 4],
        list_exclusive: [0; 4],
        poly_t: 0,
        tso: 0,
        umi_lengths: [0; UMI_MEASURED_MAX + 1],
    };
    let mut reverse = [0u8; LONG_ADAPTER_WINDOW];
    for read in reads.iter().take(LONG_READS_EXAMINED) {
        evidence.reads += 1;
        if read.len() < USABLE_BASES_MIN {
            continue;
        }
        evidence.usable += 1;
        evidence_read(read, &adapter, &tso, lists, &mut reverse, &mut evidence);
    }
    assert!(evidence.usable <= evidence.reads, "usable reads are a subset of reads");
    assert!(evidence.adapter <= evidence.usable, "adapter reads are a subset of usable reads");
    assert!(evidence.forward <= evidence.adapter, "forward reads are a subset of adapter reads");
    assert!(
        evidence.list_hits.iter().all(|&hits| hits <= evidence.adapter),
        "list hits come from adapter-bearing reads"
    );
    assert!(
        evidence.poly_t + evidence.tso <= evidence.adapter,
        "one tail per adapter read at most"
    );
    evidence
}

/// One read's contribution: the adapter at either end, the barcode after it, the tail after the
/// UMI. `reverse` is scratch for the reverse-complemented end, reused across reads.
fn evidence_read(
    read: &[u8],
    adapter: &Pattern,
    tso: &Pattern,
    lists: &Lists,
    reverse: &mut [u8; LONG_ADAPTER_WINDOW],
    evidence: &mut LongEvidence,
) {
    let head = read.get(..LONG_ADAPTER_WINDOW).unwrap_or(read);
    let tail = read.get(read.len().saturating_sub(LONG_ADAPTER_WINDOW)..).unwrap_or(read);
    let Some(scratch) = reverse.get_mut(..tail.len()) else {
        unreachable!("a read's tail window fits the scratch buffer")
    };
    for (slot, &base) in scratch.iter_mut().zip(tail.iter().rev()) {
        *slot = complement(base.to_ascii_uppercase());
    }
    let forward_hit = adapter.search(head, LONG_ADAPTER_EDITS_MAX);
    let reverse_hit = adapter.search(scratch, LONG_ADAPTER_EDITS_MAX);
    // The better match wins; a tie goes to the read's start.
    let (molecule, hit, forward) = match (forward_hit, reverse_hit) {
        (Some(f), Some(r)) if r.distance < f.distance => (&*scratch, r, false),
        (Some(f), _) => (head, f, true),
        (None, Some(r)) => (&*scratch, r, false),
        (None, None) => return,
    };
    assert!(hit.end <= molecule.len(), "the adapter ends inside its window");
    evidence.adapter += 1;
    evidence.forward += usize::from(forward);
    let barcode_end = evidence_barcode(molecule, hit.end, lists, evidence);
    evidence_tail(molecule, barcode_end, tso, evidence);
}

/// Tests the 16 bases after the adapter (start shifted up to `LONG_BARCODE_SHIFT_MAX` either way)
/// against every 10x list; counts each list at most once per read. Returns where the barcode ends,
/// at the shift that matched first (shift 0 when none did).
fn evidence_barcode(
    molecule: &[u8],
    adapter_end: usize,
    lists: &Lists,
    evidence: &mut LongEvidence,
) -> usize {
    let mut barcode_end = None;
    let mut matched = [false; 4];
    for ((slot, hit), &id) in evidence.list_hits.iter_mut().zip(&mut matched).zip(TENX_LISTS.iter())
    {
        let list = lists.get(id);
        let start = (0..=2 * LONG_BARCODE_SHIFT_MAX)
            .map(|step| shift_of(adapter_end, step))
            .filter_map(|start| start.filter(|&s| s + BARCODE_BASES <= molecule.len()))
            .find(|&s| {
                molecule.get(s..s + BARCODE_BASES).is_some_and(|b| list.hit(b) == Hit::Exact)
            });
        if let Some(start) = start {
            *slot += 1;
            *hit = true;
            barcode_end.get_or_insert(start + BARCODE_BASES);
        }
    }
    // On exactly one list: the read counts for that list alone.
    if matched.iter().filter(|&&m| m).count() == 1 {
        let Some(only) = matched.iter().position(|&m| m) else { unreachable!("one list matched") };
        let Some(slot) = evidence.list_exclusive.get_mut(only) else { unreachable!("4 lists") };
        *slot += 1;
    }
    let end = barcode_end.unwrap_or(adapter_end + BARCODE_BASES);
    assert!(end >= BARCODE_BASES, "a barcode end follows a whole barcode");
    assert!(end <= adapter_end + BARCODE_BASES + LONG_BARCODE_SHIFT_MAX, "within the shift bound");
    end
}

/// The barcode's start for search step `step`: 0, -1, +1, -2, +2, ... around the adapter's end.
fn shift_of(adapter_end: usize, step: usize) -> Option<usize> {
    assert!(step <= 2 * LONG_BARCODE_SHIFT_MAX, "steps cover the shift bound");
    let magnitude = step.div_ceil(2);
    let start = if step % 2 == 1 {
        adapter_end.checked_sub(magnitude)
    } else {
        Some(adapter_end + magnitude)
    };
    assert!(
        start.is_none_or(|s| s.abs_diff(adapter_end) <= LONG_BARCODE_SHIFT_MAX),
        "within bound"
    );
    start
}

/// After the barcode: the TSO (5' kits) or a T-run (3' kits), whichever starts first, and the UMI
/// length that puts it there.
fn evidence_tail(molecule: &[u8], barcode_end: usize, tso: &Pattern, evidence: &mut LongEvidence) {
    let Some(window) = molecule.get(barcode_end..) else { return };
    let window = window.get(..LONG_TAIL_WINDOW).unwrap_or(window);
    let tso_start =
        tso.search(window, LONG_TSO_EDITS_MAX).map(|hit| hit.end.saturating_sub(TSO.len()));
    let poly_t_start = window
        .windows(LONG_POLY_T_RUN)
        .position(|run| run.iter().all(|&b| b.eq_ignore_ascii_case(&b'T')));
    let (start, is_tso) = match (tso_start, poly_t_start) {
        (Some(t), Some(p)) if p < t => (p, false),
        (Some(t), _) => (t, true),
        (None, Some(p)) => (p, false),
        (None, None) => return,
    };
    assert!(start <= window.len(), "a tail starts inside its window");
    if start <= UMI_MEASURED_MAX {
        let Some(bin) = evidence.umi_lengths.get_mut(start) else { unreachable!("start checked") };
        *bin += 1;
    }
    if is_tso {
        evidence.tso += 1;
    } else {
        evidence.poly_t += 1;
    }
}

/// What the long-read evidence supports.
#[derive(Debug, PartialEq)]
pub enum LongVerdict {
    Match { kit: Kit, confidence: Confidence, contradictions: Vec<String> },
    Ambiguous { kits: Vec<Kit>, reason: String },
    NoMatch { reason: String },
}

impl std::fmt::Debug for Kit {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Kit({})", self.layout.name)
    }
}

impl PartialEq for Kit {
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.layout, other.layout)
    }
}

/// The most common measured UMI length, if any read showed one.
pub fn umi_mode(evidence: &LongEvidence) -> Option<usize> {
    let (length, count) = evidence
        .umi_lengths
        .iter()
        .enumerate()
        .max_by_key(|&(length, &count)| (count, std::cmp::Reverse(length)))?;
    assert!(length <= UMI_MEASURED_MAX, "the mode is a measured length");
    assert!(*count <= evidence.adapter, "UMIs are measured on adapter-bearing reads");
    (*count > 0).then_some(length)
}

/// Judges the evidence: which list, which tail, which UMI length, and so which kit.
pub fn decide(evidence: &LongEvidence, kits: &[Kit]) -> LongVerdict {
    assert!(!kits.is_empty(), "there are kits to choose from");
    let adapter_share = ratio(evidence.adapter, evidence.usable);
    if adapter_share < LONG_ADAPTER_SHARE_MIN {
        return LongVerdict::NoMatch {
            reason: format!(
                "the 10x Read 1 adapter is in {:.1}% of reads long enough to show it: not 10x \
                 single-cell cDNA, or a kit not yet supported for long reads",
                adapter_share * 100.0
            ),
        };
    }
    let list = match decide_list(evidence) {
        ListCall::One(list) => list,
        ListCall::Unclear(best, next) => {
            let tail = decide_tail(evidence);
            let kits: Vec<Kit> = kits
                .iter()
                .copied()
                .filter(|k| (k.list == best || k.list == next) && tail.is_none_or(|t| k.tail == t))
                .collect();
            let lists = format!(
                "the barcodes fit both the {} and the {} lists",
                best.label(),
                next.label()
            );
            return match kits.as_slice() {
                // The tail leaves one kit: the kit is known even though its list is not clear.
                [kit] => LongVerdict::Match {
                    kit: *kit,
                    confidence: Confidence::Medium,
                    contradictions: vec![lists],
                },
                [] => LongVerdict::NoMatch {
                    reason: format!("{lists}, and the tail fits no kit using them"),
                },
                _ => LongVerdict::Ambiguous { kits, reason: lists },
            };
        }
        ListCall::None => {
            return LongVerdict::NoMatch {
                reason: "the 10x adapter is present, but its barcodes are on no 10x \
                         gene-expression list (Multiome, Visium or a newer kit)"
                    .to_string(),
            };
        }
    };
    let tail = decide_tail(evidence);
    let on_list: Vec<Kit> = kits.iter().copied().filter(|k| k.list == list).collect();
    assert!(!on_list.is_empty(), "every tested list belongs to a kit");
    let fitting: Vec<Kit> =
        on_list.iter().copied().filter(|k| tail.is_none_or(|t| k.tail == t)).collect();
    match fitting.as_slice() {
        [kit] => decide_match(evidence, *kit, tail.is_some()),
        [] => LongVerdict::NoMatch {
            reason: "the barcodes' list and the tail after the UMI fit no known 10x kit"
                .to_string(),
        },
        many => LongVerdict::Ambiguous {
            kits: many.to_vec(),
            reason: "polyT and TSO are not separated after the UMI, so 3' and 5' cannot be told \
                     apart"
                .to_string(),
        },
    }
}

/// What the barcodes say about the list.
enum ListCall {
    /// One list, well clear of every other.
    One(ListId),
    /// Two lists the barcodes cannot separate: the best and the next.
    Unclear(ListId, ListId),
    /// No list clears the bar.
    None,
}

/// The list most adapter-bearing reads match on their own (barcodes on several lists do not
/// count), and whether it stands well clear of the next.
fn decide_list(evidence: &LongEvidence) -> ListCall {
    let shares: Vec<f64> =
        evidence.list_exclusive.iter().map(|&h| ratio(h, evidence.adapter)).collect();
    let share = |index: usize| shares.get(index).copied().unwrap_or(0.0);
    let mut order: Vec<usize> = (0..shares.len()).collect();
    order.sort_by(|&a, &b| share(b).total_cmp(&share(a)).then(a.cmp(&b)));
    let (Some(&best), Some(&next)) = (order.first(), order.get(1)) else {
        unreachable!("four lists are tested")
    };
    let (best_share, next_share) = (share(best), share(next));
    assert!(next_share <= best_share, "the best list holds the most reads");
    assert!(best < TENX_LISTS.len() && next < TENX_LISTS.len(), "indices name lists");
    let (Some(&best_id), Some(&next_id)) = (TENX_LISTS.get(best), TENX_LISTS.get(next)) else {
        unreachable!("indices name lists")
    };
    if best_share < LONG_LIST_SHARE_MIN {
        ListCall::None
    } else if best_share >= LONG_LIST_SEPARATION * next_share {
        ListCall::One(best_id)
    } else {
        ListCall::Unclear(best_id, next_id)
    }
}

/// The tail after the UMI, if one clearly dominates.
fn decide_tail(evidence: &LongEvidence) -> Option<Tail> {
    let (poly_t, tso) =
        (ratio(evidence.poly_t, evidence.adapter), ratio(evidence.tso, evidence.adapter));
    let tail = if poly_t >= LONG_TAIL_SHARE_MIN && poly_t >= LONG_TAIL_SEPARATION * tso {
        Some(Tail::PolyT)
    } else if tso >= LONG_TAIL_SHARE_MIN && tso >= LONG_TAIL_SEPARATION * poly_t {
        Some(Tail::Tso)
    } else {
        None
    };
    if tail == Some(Tail::PolyT) {
        assert!(poly_t > tso, "polyT wins only by leading");
    }
    if tail == Some(Tail::Tso) {
        assert!(tso > poly_t, "the TSO wins only by leading");
    }
    tail
}

/// A match on `kit`: HIGH needs a strong list, a decided tail and a UMI length that agrees.
fn decide_match(evidence: &LongEvidence, kit: Kit, tail_decided: bool) -> LongVerdict {
    let list_index = TENX_LISTS.iter().position(|&id| id == kit.list);
    let Some(list_index) = list_index else { unreachable!("kits use tested lists") };
    let hits = evidence.list_exclusive.get(list_index).copied().unwrap_or(0);
    let share = ratio(hits, evidence.adapter);
    let mut contradictions = Vec::new();
    if let Some(umi) = umi_mode(evidence)
        && umi.abs_diff(kit.umi) > LONG_UMI_TOLERANCE
    {
        contradictions
            .push(format!("the UMI measures {umi} bases; {} has {}", kit.layout.name, kit.umi));
    }
    let confidence = if share >= LONG_LIST_SHARE_HIGH && tail_decided && contradictions.is_empty() {
        Confidence::High
    } else {
        Confidence::Medium
    };
    assert!(share >= LONG_LIST_SHARE_MIN, "a match stands on its list");
    if confidence == Confidence::High {
        assert!(contradictions.is_empty(), "HIGH has nothing against it");
    }
    LongVerdict::Match { kit, confidence, contradictions }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::whitelist::Whitelist;
    use crate::lists::ListId;

    struct Random(u64);

    impl Random {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> usize {
            usize::try_from(self.next() % n).unwrap()
        }
        fn bases(&mut self, length: usize) -> Vec<u8> {
            (0..length).map(|_| b"ACGT"[self.below(4)]).collect()
        }
    }

    type Barcodes = Vec<(ListId, Vec<Vec<u8>>)>;

    /// Every list the layouts name, of the right lengths, each holding `per_list` random barcodes.
    fn lists(random: &mut Random, per_list: usize) -> (Lists, Barcodes) {
        let mut barcodes = Vec::new();
        let lists: Vec<(ListId, Whitelist)> = ListId::ALL
            .iter()
            .map(|&id| {
                let length = id.barcode_length();
                let codes: Vec<Vec<u8>> = (0..per_list).map(|_| random.bases(length)).collect();
                let slices: Vec<&[u8]> = codes.iter().map(Vec::as_slice).collect();
                let list = Whitelist::from_barcodes(&slices);
                barcodes.push((id, codes));
                (id, list)
            })
            .collect();
        (Lists::new(lists), barcodes)
    }

    /// Nanopore-like copying: 3% substitutions, 2% insertions, 2% deletions per base.
    fn with_errors(random: &mut Random, molecule: &[u8]) -> Vec<u8> {
        let mut read = Vec::with_capacity(molecule.len() + 64);
        for &base in molecule {
            match random.below(100) {
                0..=2 => read.push(b"ACGT"[random.below(4)]),
                3..=4 => read.extend([base, b"ACGT"[random.below(4)]]),
                5..=6 => {}
                _ => read.push(base),
            }
        }
        read
    }

    fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
        sequence.iter().rev().map(|&b| complement(b)).collect()
    }

    /// Long reads of `kit`: adapter, a listed barcode, UMI, tail, a transcript; half reversed.
    fn long_reads(random: &mut Random, kit: &Kit, barcodes: &[Vec<u8>], count: usize) -> Reads {
        let mut reads = Reads::with_capacity(count, count * 1_200);
        for _ in 0..count {
            let mut molecule = random.bases(20);
            molecule.extend_from_slice(READ1_ADAPTER);
            molecule.extend_from_slice(&barcodes[random.below(barcodes.len() as u64)]);
            molecule.extend(random.bases(kit.umi));
            match kit.tail {
                Tail::PolyT => molecule.extend([b'T'; 30]),
                Tail::Tso => molecule.extend_from_slice(TSO),
            }
            molecule.extend(random.bases(900));
            let read = with_errors(random, &molecule);
            reads.push(&if random.below(2) == 0 { read } else { reverse_complement(&read) });
        }
        reads
    }

    #[test]
    fn names_every_10x_kit_from_simulated_nanopore_reads() {
        let mut random = Random(5);
        let (lists, barcodes) = lists(&mut random, 400);
        let kits = kits();
        assert_eq!(kits.len(), 5, "3' v2, v3, v4 and 5' v1/v2, v3");
        for kit in &kits {
            let listed = &barcodes.iter().find(|(id, _)| *id == kit.list).unwrap().1;
            let reads = long_reads(&mut random, kit, listed, 2_000);
            let evidence = evidence(&reads, &lists);
            assert!(evidence.adapter * 10 > evidence.reads * 8, "{kit:?}: {evidence:?}");
            match decide(&evidence, &kits) {
                LongVerdict::Match { kit: found, confidence, contradictions } => {
                    assert_eq!(found, *kit, "{evidence:?}");
                    assert_eq!(
                        confidence,
                        Confidence::High,
                        "{kit:?}: {contradictions:?} {evidence:?}"
                    );
                }
                other => panic!("{kit:?}: {other:?} from {evidence:?}"),
            }
        }
    }

    #[test]
    fn random_long_reads_are_no_match() {
        let mut random = Random(9);
        let (lists, _) = lists(&mut random, 400);
        let mut reads = Reads::with_capacity(500, 500 * 1_000);
        for _ in 0..500 {
            reads.push(&random.bases(1_000));
        }
        let evidence = evidence(&reads, &lists);
        assert!(matches!(decide(&evidence, &kits()), LongVerdict::NoMatch { .. }), "{evidence:?}");
    }
}
