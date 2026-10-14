// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! ScaleBio single-cell RNA kits: the original RNA kit (v1.0/v1.1), and `QuantumScale`.

use super::bc;
use super::{Layout, ListId, Mate, Segment, Tail};
use Segment::{Fixed, Umi, Unlisted};

/// Read off public runs with scripts/discover.py (7 Oct 2026): ligation barcode 8-9 bases, then
/// TCAGAGC, UMI 8, RT barcode 10, polyT.
pub const SCALEBIO: Layout = Layout {
    name: "ScaleBio RNA",
    read: Mate::R1,
    offsets: &[0, 1],
    segments: &[
        Unlisted { len: 8, set: 384, cover_pct_min: 50 },
        Fixed(&[b"TCAGAGC"]),
        Umi(8),
        Unlisted { len: 10, set: 384, cover_pct_min: 50 },
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{b[8-9]f[TCAGAGC]u[8]b[10]x:}2{r:}",
    note: "cell = ligation barcode + RT barcode (+ PCR index); the ligation barcode is 8-9 \
        bases",
    long_read_kit: None,
    cellranger_chemistry: None,
};

/// `QuantumScale` (`ScaleRna` `libQuantumV1.0`): R1 is cDNA; R2 is a 17-base read of an 8-base
/// UMI then the 8-base RT barcode; the three bead barcodes are in I1 (bases 0, 12 and 24, split by
/// CTGT and GGTT) and a PCR barcode in I2. Checked 8 Oct 2026 on ScaleBio's test data: RT 98% at
/// R2 8-16. R2 alone names the kit, which matters because archives usually drop the index reads;
/// the bead barcodes in I1 are not checked here.
pub const QUANTUMSCALE: Layout = Layout {
    name: "ScaleBio QuantumScale RNA",
    read: Mate::R2,
    offsets: &[0],
    segments: &[Umi(8), bc(ListId::ScaleQuantumRt, 8)],
    geometry: "1{r:}2{u[8]b[8]x:}",
    note: "cell = RT barcode (R2 8-16) + three bead barcodes in I1 + PCR barcode in I2; give the \
           index reads to ScaleBio's pipeline, which needs them",
    long_read_kit: None,
    cellranger_chemistry: None,
};
