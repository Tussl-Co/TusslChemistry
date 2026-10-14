// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Candidate library layouts. A layout names the read that carries the cell barcode, the start
//! offsets that read may have (BD Enhanced Beads prefix 0-3 variable bases), and the ordered
//! segments from that start: barcodes (with a known list, or without one), UMI, fixed linkers and
//! the tail after the UMI.
//!
//! One entry per processing-equivalent layout: kits whose reads and processing are identical share
//! an entry (10x 3' v3 and v3.1 differ in wet-lab chemistry, not in their reads).
//! Geometry strings use simpleaf/piscem's notation, so a result can drive processing directly.

mod bd;
mod mgi;
mod other;
mod pipseq;
mod scalebio;
mod seekgene;
mod singleron;
mod split_parse;
mod tenx;
pub mod validation;

use crate::lists::{ListId, ListSet};

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Tail {
    /// Oligo-dT capture: 3' kits.
    PolyT,
    /// Template-switch oligo: 5' kits.
    Tso,
}

pub const TSO: &[u8] = b"TTTCTTATATGGG";
/// Largest start offset any layout allows (BD Enhanced Beads' 0-3 variable bases). Sizes the
/// per-offset counters; `layouts_are_sound` checks every layout against it.
pub const OFFSET_MAX: usize = 3;
/// Longest tail compared; the TSO's length.
pub const TAIL_MAX: usize = 13;
/// Most segments one layout may have; sizes per-layout tallies.
pub const SEGMENTS_MAX: usize = 10;
/// Longest fixed-length segment: SPLiT-seq's 30-base linkers.
const SEGMENT_BASES_MAX: usize = 32;
/// Most bases from a layout's start offset to the end of its last fixed-length segment (SPLiT-seq
/// spans 94).
const STRUCTURE_BASES_MAX: usize = 128;
const POLY_T: [u8; TAIL_MAX] = [b'T'; TAIL_MAX];
const _: () = assert!(TSO.len() == TAIL_MAX, "both tail patterns must be TAIL_MAX long");

impl Tail {
    /// The first `len` bases this tail should show. Static, so the per-read check allocates
    /// nothing.
    pub fn expected(self, len: usize) -> &'static [u8] {
        assert!(len <= TAIL_MAX, "tail window {len} past TAIL_MAX {TAIL_MAX}");
        let full: &'static [u8] = match self {
            Tail::PolyT => &POLY_T,
            Tail::Tso => TSO,
        };
        // An empty window would match vacuously: the scorer compares only windows of 2+ bases.
        debug_assert!(len >= 1, "a tail window has at least one base"); // Hot: checked build only.
        #[allow(clippy::indexing_slicing)]
        // earned by the assertions above: both patterns are TAIL_MAX long
        &full[..len]
    }

    pub fn label(self) -> &'static str {
        match self {
            Tail::PolyT => "polyT",
            Tail::Tso => "TSO",
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Segment {
    /// A cell barcode from a known list.
    Barcode {
        list: ListId,
        len: usize,
    },
    /// A barcode whose list is not public (Parse round 1) or does not exist (Drop-seq beads).
    /// Real barcodes repeat, so the `set` most common sequences must cover at least
    /// `cover_pct_min` percent of reads; random sequence covers far less.
    Unlisted {
        len: usize,
        set: usize,
        cover_pct_min: u8,
    },
    Umi(usize),
    /// A linker: one of these sequences, all the same length (most kits have one; Singleron's
    /// beads carry one of four at each linker position).
    Fixed(&'static [&'static [u8]]),
    Tail(Tail),
}

impl Segment {
    /// Bases the segment occupies; the tail is open-ended.
    pub const fn len(self) -> usize {
        let len = match self {
            Segment::Barcode { len, .. } | Segment::Unlisted { len, .. } | Segment::Umi(len) => len,
            Segment::Fixed(alternatives) => match alternatives {
                [first, ..] => first.len(),
                [] => 0,
            },
            Segment::Tail(_) => 0,
        };
        // Hot: checked build only. Both are compile-time table checks too (layouts_are_sound).
        debug_assert!(len <= SEGMENT_BASES_MAX, "a segment is at most SEGMENT_BASES_MAX bases");
        if !matches!(self, Segment::Tail(_)) {
            debug_assert!(len >= 1, "a fixed-length segment occupies at least one base");
        }
        len
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Mate {
    R1,
    R2,
}

pub struct Layout {
    pub name: &'static str,
    /// Which read carries the barcode.
    pub read: Mate,
    pub offsets: &'static [usize],
    pub segments: &'static [Segment],
    pub geometry: &'static str,
    /// What processing must know beyond the geometry.
    pub note: &'static str,
    /// The kit as Oxford Nanopore's `wf-single-cell` names it (`--kit`), for 10x kits whose cDNA
    /// is also read whole on long-read platforms; None for kits it does not support.
    pub long_read_kit: Option<&'static str>,
    /// Cell Ranger's `--chemistry` for this layout, as `cellranger count` 10.1.0 accepts it
    /// (checked against its argument parser, bench/emit); None where Cell Ranger has no
    /// `count` chemistry for it (other vendors, and Flex, which needs `cellranger multi`).
    pub cellranger_chemistry: Option<&'static str>,
}

impl Layout {
    /// The first allowed start offset (0 for every layout without a variable prefix).
    pub fn first_offset(&self) -> usize {
        let [first, ..] = self.offsets else {
            unreachable!("checked at compile time: layout {} has a start offset", self.name)
        };
        assert!(*first <= OFFSET_MAX, "{}: offset {first} past OFFSET_MAX", self.name);
        *first
    }

    /// The tail after the UMI, when the layout reads that far. Pairs with `layouts_are_sound`:
    /// a tail is the last segment and there is at most one.
    pub fn tail(&self) -> Option<Tail> {
        let tail = self.segments.iter().find_map(|s| match s {
            Segment::Tail(t) => Some(*t),
            _ => None,
        });
        if let Some(t) = tail {
            let last = self.segments.last();
            assert_eq!(
                last,
                Some(&Segment::Tail(t)),
                "{}: the tail is the last segment",
                self.name
            );
        }
        assert!(self.segments.len() <= SEGMENTS_MAX, "{}: too many segments", self.name);
        tail
    }

    /// Bases from the start offset to the end of the last fixed-length segment.
    pub fn structured_len(&self) -> usize {
        let bases: usize = self.segments.iter().map(|s| s.len()).sum();
        assert!(bases >= 1, "{}: a layout occupies at least one base", self.name);
        assert!(bases <= STRUCTURE_BASES_MAX, "{}: {bases} bases of structure", self.name);
        bases
    }

    /// The barcode lists the layout reads; lists decide which layouts tie. Singleron V1 draws
    /// all three barcodes from one list, which the set holds once.
    pub const fn lists(&self) -> ListSet {
        assert!(self.segments.len() <= SEGMENTS_MAX, "a layout has at most SEGMENTS_MAX segments");
        let (mut set, mut barcodes, mut i) = (ListSet::EMPTY, 0, 0);
        while i < self.segments.len() {
            #[allow(clippy::indexing_slicing)] // earned by the loop bound
            if let Segment::Barcode { list, .. } = self.segments[i] {
                set = set.with(list);
                barcodes += 1;
            }
            i += 1;
        }
        assert!(set.len() <= barcodes, "a set never holds more lists than barcode segments");
        set
    }

    /// The fixed linkers (each its alternatives), in read order; linkers decide which layouts tie.
    pub fn linkers(&self) -> impl Iterator<Item = &'static [&'static [u8]]> + '_ {
        assert!(self.segments.len() <= SEGMENTS_MAX, "{}: too many segments", self.name);
        // An empty linker would match every read; `linker_is_sound` refuses one at compile time.
        assert!(
            self.segments.iter().all(|s| !matches!(s, Segment::Fixed(a) if a.is_empty())),
            "{}: a linker has at least one sequence",
            self.name
        );
        self.segments.iter().filter_map(|s| match s {
            Segment::Fixed(alternatives) => Some(*alternatives),
            _ => None,
        })
    }
}

use Segment::Barcode;

/// A listed barcode segment. Evaluated at compile time over `LAYOUTS`, so a bad length stops the
/// build: the whitelist packs at most 16 bases into a u32.
const fn bc(list: ListId, len: usize) -> Segment {
    assert!(len >= 1, "a barcode has at least one base");
    assert!(
        len <= crate::input::whitelist::LEN_MAX,
        "a listed barcode fits the whitelist's 16 bases"
    );
    assert!(len == list.barcode_length(), "a barcode segment is its list's length");
    Barcode { list, len }
}

/// Every supported layout, in the order candidates are listed when their support ties.
pub const LAYOUTS: &[Layout] = &[
    tenx::THREE_PRIME_V2,
    tenx::FIVE_PRIME_V1_V2,
    tenx::THREE_PRIME_V3,
    tenx::THREE_PRIME_V4,
    tenx::FIVE_PRIME_V3,
    tenx::MULTIOME_GEX,
    tenx::FLEX_V1,
    bd::V1_BEADS,
    bd::ENHANCED_BEADS,
    bd::ENHANCED_V2_V3_BEADS,
    split_parse::SPLIT_SEQ,
    split_parse::PARSE_WT,
    split_parse::PARSE_WT_V3,
    other::INDROPS,
    scalebio::SCALEBIO,
    scalebio::QUANTUMSCALE,
    other::DROP_SEQ,
    singleron::GEXSCOPE_V1,
    singleron::GEXSCOPE_V2,
    singleron::GEXSCOPE_V3,
    mgi::C_SERIES,
    pipseq::V3_V4,
    pipseq::V5,
    seekgene::DD_THREE_PRIME,
    seekgene::DD_FIVE_PRIME,
];

/// Whether the layout's second read is a probe, not cDNA (Flex): no whole-transcriptome counter
/// fits it, only Cell Ranger's `multi` with the probe set.
pub fn probe_based(layout: &Layout) -> bool {
    let probes = layout.name == tenx::FLEX_V1.name;
    assert!(!probes || layout.cellranger_chemistry.is_none(), "a probe kit has no count chemistry");
    assert!(!layout.name.is_empty(), "a layout has a name");
    probes
}

/// Whether the layout's library is always made from nuclei (Multiome): much of a nucleus's RNA is
/// unspliced, so a counter must count intronic reads too, as Cell Ranger ARC does.
pub fn nuclei(layout: &Layout) -> bool {
    let nuclei = layout.name == tenx::MULTIOME_GEX.name;
    assert!(!nuclei || layout.cellranger_chemistry == Some("ARC-v1"), "Multiome counts as ARC");
    assert!(!layout.name.is_empty(), "a layout has a name");
    nuclei
}

/// A set of layouts, one bit per `LAYOUTS` entry: fixed size, copied by value.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct LayoutSet(u32);

const _: () = assert!(LAYOUTS.len() <= u32::BITS as usize, "every layout has a bit in a LayoutSet");

impl LayoutSet {
    pub const EMPTY: LayoutSet = LayoutSet(0);

    /// The set with the layout at `index` in `LAYOUTS` added.
    #[must_use]
    pub fn with(self, index: usize) -> LayoutSet {
        assert!(index < LAYOUTS.len(), "layout {index} is in the table");
        let set = LayoutSet(self.0 | 1 << index);
        assert!(set.0 & self.0 == self.0, "adding a layout keeps the others");
        set
    }

    pub fn is_empty(self) -> bool {
        assert!(self.0 >> LAYOUTS.len() == 0, "a set holds only table layouts");
        let empty = self.0 == 0;
        assert!(empty == (self.0.count_ones() == 0), "an empty set has no layouts");
        empty
    }

    /// The layouts in the set, in table order.
    pub fn iter(self) -> impl Iterator<Item = &'static Layout> {
        assert!(self.0 >> LAYOUTS.len() == 0, "a set holds only table layouts");
        assert!(self.0.count_ones() as usize <= LAYOUTS.len(), "at most every layout");
        LAYOUTS.iter().enumerate().filter(move |(i, _)| self.0 & 1 << i != 0).map(|(_, l)| l)
    }
}

/// Most alternatives one linker may have (Singleron's beads have four).
const LINKER_ALTERNATIVES_MAX: usize = 8;

/// A linker's alternatives: at least one, at most `LINKER_ALTERNATIVES_MAX`, all the same length,
/// since the scorer reads one window per linker.
const fn linker_is_sound(alternatives: &[&[u8]]) {
    assert!(!alternatives.is_empty(), "a linker has at least one sequence");
    assert!(alternatives.len() <= LINKER_ALTERNATIVES_MAX, "a linker has at most 8 alternatives");
    let mut k = 1;
    while k < alternatives.len() {
        #[allow(clippy::indexing_slicing)] // earned by the loop bound and the non-empty check
        let same = alternatives[k].len() == alternatives[0].len();
        assert!(same, "a linker's alternatives are all the same length");
        k += 1;
    }
}

/// Compile-time checks on the layout table, one assertion per property, so a bad entry stops the
/// build with the property it breaks. `first_offset`, `tail` and `Segment::len` assert the same
/// properties at run time.
const fn layouts_are_sound() -> bool {
    let mut i = 0;
    while i < LAYOUTS.len() {
        #[allow(clippy::indexing_slicing)] // earned by the loop bound
        let layout = &LAYOUTS[i];
        assert!(!layout.offsets.is_empty(), "every layout has a start offset");
        assert!(
            layout.segments.len() <= SEGMENTS_MAX,
            "a layout has at most SEGMENTS_MAX segments"
        );
        let mut o = 0;
        while o < layout.offsets.len() {
            #[allow(clippy::indexing_slicing)] // earned by the loop bound
            let offset = layout.offsets[o];
            assert!(offset <= OFFSET_MAX, "every start offset is at most OFFSET_MAX");
            o += 1;
        }
        let (mut j, mut umis, mut bases) = (0, 0, 0);
        while j < layout.segments.len() {
            #[allow(clippy::indexing_slicing)] // earned by the loop bound
            let segment = &layout.segments[j];
            match segment {
                Segment::Umi(_) => umis += 1,
                // The scorer reads the tail after every other segment.
                Segment::Tail(_) => assert!(j + 1 == layout.segments.len(), "a tail is last"),
                Segment::Fixed(alternatives) => linker_is_sound(alternatives),
                _ => {}
            }
            let len = segment.len();
            // A zero-length barcode, UMI or linker would check nothing.
            if !matches!(segment, Segment::Tail(_)) {
                assert!(len >= 1, "a barcode, UMI or linker occupies at least one base");
            }
            assert!(len <= SEGMENT_BASES_MAX, "a segment is at most SEGMENT_BASES_MAX bases");
            bases += len;
            j += 1;
        }
        // The scorer keeps one set of UMI position counts per layout.
        assert!(umis <= 1, "a layout has at most one UMI");
        assert!(bases >= 1, "a layout occupies at least one base");
        assert!(bases <= STRUCTURE_BASES_MAX, "a layout's structure fits STRUCTURE_BASES_MAX");
        i += 1;
    }
    true
}

const _: () = assert!(layouts_are_sound(), "the layout table is unsound");

#[cfg(test)]
mod tests {
    use super::*;

    /// The table check also runs here, so it still guards the table if the compile-time
    /// assertion that calls it is ever lost.
    #[test]
    fn the_table_is_sound() {
        assert!(layouts_are_sound());
    }

    /// Every layout answers its table questions without tripping an invariant: the lists it
    /// uses (each named once), its linkers, its tail and its first offset.
    #[test]
    fn every_layout_describes_itself() {
        for layout in LAYOUTS {
            assert!(layout.lists().len() <= SEGMENTS_MAX, "{}", layout.name);
            assert!(layout.linkers().count() <= SEGMENTS_MAX, "{}", layout.name);
            let _ = layout.tail();
            assert!(layout.first_offset() <= OFFSET_MAX, "{}", layout.name);
        }
        assert_eq!(
            LAYOUTS.iter().find(|l| l.name == "Singleron GEXSCOPE V1").map(Layout::lists),
            Some(ListSet::EMPTY.with(ListId::SingleronV1)),
            "one list at three positions is named once"
        );
    }
}
