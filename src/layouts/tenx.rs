// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! 10x Genomics Chromium 3' and 5' gene expression, v2 to v4 (GEM-X), the gene expression library
//! of Single Cell Multiome ATAC + Gene Expression, and Flex (fixed RNA, probe-based).

use super::ListId;
use super::bc;
use super::{Layout, Mate, Segment, Tail};
use Segment::Umi;

pub const THREE_PRIME_V2: Layout = Layout {
    name: "10x 3' v2",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxV2, 16), Umi(10), Segment::Tail(Tail::PolyT)],
    geometry: "1{b[16]u[10]x:}2{r:}",
    note: "",
    long_read_kit: Some("3prime:v2"),
    cellranger_chemistry: Some("SC3Pv2"),
};

pub const FIVE_PRIME_V1_V2: Layout = Layout {
    name: "10x 5' v1/v2",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxV2, 16), Umi(10), Segment::Tail(Tail::Tso)],
    geometry: "1{b[16]u[10]x:}2{r:}",
    note: "",
    long_read_kit: Some("5prime:v2"),
    cellranger_chemistry: Some("SC5P-R2"),
};

pub const THREE_PRIME_V3: Layout = Layout {
    name: "10x 3' v3/v3.1",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxV3, 16), Umi(12), Segment::Tail(Tail::PolyT)],
    geometry: "1{b[16]u[12]x:}2{r:}",
    note: "",
    long_read_kit: Some("3prime:v3"),
    cellranger_chemistry: Some("SC3Pv3"),
};

pub const THREE_PRIME_V4: Layout = Layout {
    name: "10x 3' v4 (GEM-X)",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxV4ThreePrime, 16), Umi(12), Segment::Tail(Tail::PolyT)],
    geometry: "1{b[16]u[12]x:}2{r:}",
    note: "",
    long_read_kit: Some("3prime:v4"),
    cellranger_chemistry: Some("SC3Pv4"),
};

/// TSO position from 10x's published 5' GEM-X read structure. Short reads: 4 of 4 vendor runs at
/// HIGH (`bench/vendor_gold`). Long reads: 6 of 6 development runs with a stated kit at HIGH
/// (`bench/long_read`); the 2 held-out runs called it have no stated kit.
pub const FIVE_PRIME_V3: Layout = Layout {
    name: "10x 5' v3 (GEM-X)",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxV3FivePrime, 16), Umi(12), Segment::Tail(Tail::Tso)],
    geometry: "1{b[16]u[12]x:}2{r:}",
    note: "",
    long_read_kit: Some("5prime:v3"),
    cellranger_chemistry: Some("SC5P-R2-v3"),
};

/// Multiome's gene expression library reads like 3' v3 but draws its barcodes from its own list
/// (Cell Ranger ARC's 737K-arc-v1), which is how the two are told apart. Long reads of it are
/// not yet checked, so the long-read detector does not offer it.
pub const MULTIOME_GEX: Layout = Layout {
    name: "10x Multiome GEX",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxArcGex, 16), Umi(12), Segment::Tail(Tail::PolyT)],
    geometry: "1{b[16]u[12]x:}2{r:}",
    note: "the paired ATAC library's barcodes map to these through Cell Ranger ARC's list",
    long_read_kit: None,
    cellranger_chemistry: Some("ARC-v1"),
};

/// Flex v1 (Fixed RNA Profiling, Next GEM and GEM-X): R1 is a barcode from its own list and a
/// 12-base UMI; R2 reads a ligated probe pair, with an 8-base probe barcode at bases 69 to 76 that
/// names the sample. Its list ships only with Cell Ranger, so it is optional.
pub const FLEX_V1: Layout = Layout {
    name: "10x Flex v1",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::TenxFlexV1, 16), Umi(12)],
    geometry: "1{b[16]u[12]x:}2{r:}",
    note: "probe-based: R2 reads probes, not cDNA; quantify with Cell Ranger multi and the \
           probe set (a whole-transcriptome index does not fit)",
    long_read_kit: None,
    cellranger_chemistry: None,
};
