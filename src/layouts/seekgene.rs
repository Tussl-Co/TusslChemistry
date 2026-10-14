// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `SeekGene` `SeekOne` single-cell RNA, from `SeekGene`'s own pipeline (`SeekSoul` Tools, MIT
//! licence, `chemistry.py` at a9a508e), confirmed on one design run each (9 Oct 2026,
//! bench/seekgene/RESEARCH.md):
//! - DD 3' (`DDV2`): R1 = 17-base barcode, 12-base UMI, polyT (polyT from base 29 in SRR27533148).
//! - DD 5' (`DD5V1`): the same barcode and UMI, then the template-switch oligo TTTCTTATATGGG (at
//!   bases 29-41 in SRR23931570).
//! - MM (also `DDV1`): a 1-4 base anchor ending in A, then three 8-base barcodes from one 479-entry
//!   list split by linkers ACTGGTGA and GGTAGTGACA, an 8-base UMI and polyT (the linkers in 94% of
//!   SRR22784277's reads; anchors of 1-4 bases in 99.6%).

use super::ListId;
use super::bc;
use super::{Layout, Mate, Segment, Tail};
use Segment::{Umi, Unlisted};

/// The 17-base barcode's first 7 bases come from 960 sequences. They are checked by
/// concentration, not a list: allowing one substitution, 960 7-mers cover nearly every 7-mer, so a
/// list would match any read (`src/lists/seekgene.rs`). The kit's reads put almost all their
/// 7-mers on the 960 (97% and 95% in the two design runs), so the bar is higher than the 50% other
/// unlisted parts use: at 50%, a half-and-half mix with another kit's reads (52% here, the list
/// part 45%) still averaged HIGH support (bench/chaos, `mix_with_other_chemistry`).
const P3CB_PART_1: usize = 960;
const P3CB_PART_1_COVER_PCT_MIN: u8 = 80;

pub const DD_THREE_PRIME: Layout = Layout {
    name: "SeekGene SeekOne DD 3'",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        Unlisted { len: 7, set: P3CB_PART_1, cover_pct_min: P3CB_PART_1_COVER_PCT_MIN },
        bc(ListId::SeekgeneP3cb2, 10),
        Umi(12),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[17]u[12]x:}2{r:}",
    note: "cell = the 17-base barcode; SeekGene's SeekSoul Tools (--chemistry DDV2) processes it",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const DD_FIVE_PRIME: Layout = Layout {
    name: "SeekGene SeekOne DD 5'",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        Unlisted { len: 7, set: P3CB_PART_1, cover_pct_min: P3CB_PART_1_COVER_PCT_MIN },
        bc(ListId::SeekgeneP3cb2, 10),
        Umi(12),
        Segment::Tail(Tail::Tso),
    ],
    geometry: "1{b[17]u[12]x:}2{r:}",
    note: "cell = the 17-base barcode; SeekGene's SeekSoul Tools (--chemistry DD5V1) processes it",
    long_read_kit: None,
    cellranger_chemistry: None,
};
