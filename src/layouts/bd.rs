// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! BD Rhapsody whole-transcriptome beads: the original V1 beads, the Enhanced Beads, and the
//! Enhanced V2/V3 beads (WTA Next, OMICS-One).

use super::ListId;
use super::bc;
use super::{Layout, Mate, Segment, Tail};
use Segment::{Fixed, Umi};

/// The Enhanced and Enhanced V2/V3 beads' linkers between the three cell label sections.
const ENHANCED_L1: &[u8] = b"GTGA";
const ENHANCED_L2: &[u8] = b"GACA";

pub const V1_BEADS: Layout = Layout {
    name: "BD Rhapsody (original V1 beads)",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        bc(ListId::BdCls1, 9),
        Fixed(&[b"ACTGGCCTGCGA"]),
        bc(ListId::BdCls2, 9),
        Fixed(&[b"GGTAGCGGTGACA"]),
        bc(ListId::BdCls3, 9),
        Umi(8),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[9]f[ACTGGCCTGCGA]b[9]f[GGTAGCGGTGACA]b[9]u[8]x:}2{r:}",
    note: "cell = CLS1+CLS2+CLS3",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const ENHANCED_BEADS: Layout = Layout {
    name: "BD Rhapsody (Enhanced Beads)",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        bc(ListId::BdCls1, 9),
        Fixed(&[ENHANCED_L1]),
        bc(ListId::BdCls2, 9),
        Fixed(&[ENHANCED_L2]),
        bc(ListId::BdCls3, 9),
        Umi(8),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{x[0-3]b[9]f[GTGA]b[9]f[GACA]b[9]u[8]x:}2{r:}",
    note: "cell = CLS1+CLS2+CLS3; R1 starts with 0-3 variable bases (none, A, GT, TCA)",
    long_read_kit: None,
    cellranger_chemistry: None,
};

/// Enhanced V2/V3 beads (BD's pipeline docs and cell-label file): the Enhanced structure with 384
/// sequences per cell label section; V2 and V3 differ only in the 5' capture oligo, which a
/// whole-transcriptome R1 never shows. BD's WTA Next demos hit these lists at 88-92% and the
/// 96-sequence lists at 0%.
pub const ENHANCED_V2_V3_BEADS: Layout = Layout {
    name: "BD Rhapsody (Enhanced V2/V3 beads)",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        bc(ListId::BdV2Cls1, 9),
        Fixed(&[ENHANCED_L1]),
        bc(ListId::BdV2Cls2, 9),
        Fixed(&[ENHANCED_L2]),
        bc(ListId::BdV2Cls3, 9),
        Umi(8),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{x[0-3]b[9]f[GTGA]b[9]f[GACA]b[9]u[8]x:}2{r:}",
    note: "cell = CLS1+CLS2+CLS3 from the 384-sequence lists; R1 starts with 0-3 variable bases \
           set by CLS1; BD's pipeline names these Enhanced V2/V3 beads",
    long_read_kit: None,
    cellranger_chemistry: None,
};
