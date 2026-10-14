// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Kits whose barcodes have no public list, so barcodes are judged by concentration: inDrops and
//! Drop-seq.

use super::{Layout, Mate, Segment, Tail};
use Segment::{Fixed, Umi, Unlisted};

/// Read off public runs with scripts/discover.py (7 Oct 2026): BC1 is 8-11 bases, so the W1 linker
/// starts at 8-11; the variable start is modelled as offsets 0-3 on BC1's last 8 bases.
pub const INDROPS: Layout = Layout {
    name: "inDrops",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        Unlisted { len: 8, set: 384, cover_pct_min: 50 },
        Fixed(&[b"GAGTGATTGCTTGTGACGCCTT"]),
        Unlisted { len: 8, set: 384, cover_pct_min: 50 },
        Umi(6),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[8-11]f[GAGTGATTGCTTGTGACGCCTT]b[8]u[6]x:}2{r:}",
    note: "cell = BC1+BC2; BC1 is 8-11 bases (Klein lab gel bead barcodes)",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const DROP_SEQ: Layout = Layout {
    name: "Drop-seq",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        Unlisted { len: 12, set: 1000, cover_pct_min: 10 },
        Umi(8),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[12]u[8]x:}2{r:}",
    note: "no barcode list: bead barcodes are random, so cells are called by read count \
           (knee); beads often carry a synthesis error near the barcode's end",
    long_read_kit: None,
    cellranger_chemistry: None,
};
