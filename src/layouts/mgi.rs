// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! MGI `DNBelab` C-series single-cell RNA (C4 and the high-throughput kit), from MGI's own
//! read-structure config (MIT licence): R1 is two 10-base barcodes from one table, then a 10-base
//! UMI; R2 is cDNA. Both kits read alike, so they share one entry. Checked 8 Oct 2026 on public
//! runs: a 30-base R1 whose barcodes hit the table at 72-84% allowing one substitution; the misses
//! are diffuse, not a second table.

use super::bc;
use super::{Layout, ListId, Mate, Segment};
use Segment::Umi;

pub const C_SERIES: Layout = Layout {
    name: "MGI DNBelab C4 / C-series HT",
    read: Mate::R1,
    offsets: &[0],
    segments: &[bc(ListId::MgiCSeries, 10), bc(ListId::MgiCSeries, 10), Umi(10)],
    geometry: "1{b[10]b[10]u[10]x:}2{r:}",
    note: "cell = barcode 1 + barcode 2, each from one 1,536-barcode table; MGI's dnbc4tools \
           processes it; archives often hold the companion oligo library beside the cDNA",
    long_read_kit: None,
    cellranger_chemistry: None,
};
