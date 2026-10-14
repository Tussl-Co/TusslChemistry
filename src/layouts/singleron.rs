// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Singleron GEXSCOPE kits, from Singleron's own pipeline (`CeleScope`, `chemistry_dict.py`, MIT
//! licence): three barcodes split by linkers, a fixed C, a 12-base UMI and polyT on R1. V1 and V2
//! linkers come in four alternatives each (one per bead set); V3 shortens them to 6 bases and
//! starts 0 to 3 bases in. Checked 7 Oct 2026 on public runs: the base before the UMI is C in 83%
//! to 99% of reads, and polyT is often degraded in V1/V2 runs, so a missing tail is only a note.
//! `GEXSCOPE` `MicroBead` (12-base barcode, 8-base UMI, no list or linker) is structurally
//! Drop-seq, and is reported as Drop-seq.

use super::bc;
use super::{Layout, ListId, Mate, Segment, Tail};
use Segment::{Fixed, Umi};

/// `GEXSCOPE` V1 and V2 linker alternatives (`CeleScope` `linker1.txt`, `linker2.txt`).
const LINKER1: &[&[u8]] =
    &[b"ATCCACGTGCTTGAGA", b"TCGGTGACAGCCATAT", b"CGAACATGTAGGTCTC", b"GATTGTCACTAACGCG"];
const LINKER2: &[&[u8]] =
    &[b"TCAGCATGCGGCTACG", b"CGTAGTCAGAAGCTGA", b"GACTACGTATTAGCAT", b"ATGCTGACTCCTAGTC"];
/// The fixed base between the third barcode and the UMI ("L1" in `CeleScope`'s patterns).
const BEFORE_UMI: &[&[u8]] = &[b"C"];

pub const GEXSCOPE_V1: Layout = Layout {
    name: "Singleron GEXSCOPE V1",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        bc(ListId::SingleronV1, 8),
        Fixed(LINKER1),
        bc(ListId::SingleronV1, 8),
        Fixed(LINKER2),
        bc(ListId::SingleronV1, 8),
        Fixed(BEFORE_UMI),
        Umi(12),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[8]x[16]b[8]x[16]b[8]x[1]u[12]x:}2{r:}",
    note: "cell = bc1+bc2+bc3, one 96-barcode list at all three; linkers vary by bead set \
           (CeleScope chemistry GEXSCOPE-V1 handles them)",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const GEXSCOPE_V2: Layout = Layout {
    name: "Singleron GEXSCOPE V2",
    read: Mate::R1,
    offsets: &[0],
    segments: &[
        bc(ListId::SingleronBc1, 9),
        Fixed(LINKER1),
        bc(ListId::SingleronBc2, 9),
        Fixed(LINKER2),
        bc(ListId::SingleronBc3, 9),
        Fixed(BEFORE_UMI),
        Umi(12),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[9]x[16]b[9]x[16]b[9]x[1]u[12]x:}2{r:}",
    note: "cell = bc1+bc2+bc3; linkers vary by bead set (CeleScope chemistry GEXSCOPE-V2)",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const GEXSCOPE_V3: Layout = Layout {
    name: "Singleron GEXSCOPE V3",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        bc(ListId::SingleronBc1, 9),
        Fixed(&[b"ACGATG"]),
        bc(ListId::SingleronBc2, 9),
        Fixed(&[b"CATAGT"]),
        bc(ListId::SingleronBc3, 9),
        Fixed(BEFORE_UMI),
        Umi(12),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{x[0-3]b[9]f[ACGATG]b[9]f[CATAGT]b[9]x[1]u[12]x:}2{r:}",
    note: "cell = bc1+bc2+bc3; R1 starts with 0-3 variable bases (CeleScope chemistry \
           GEXSCOPE-V3)",
    long_read_kit: None,
    cellranger_chemistry: None,
};
