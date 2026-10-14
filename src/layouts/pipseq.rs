// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Fluent `PIPseq`, sold by Illumina as Single Cell 3' RNA Prep. The barcode read opens with 0 to 3
//! phasing bases, then four barcode parts split by three short linkers, a 12-base UMI and polyT
//! (structure from `PIPseeker`'s own output, decoded by Xi Chen, 2023). Read off one public design
//! run per version (8 Oct 2026): v3 and v4 share it, the linkers in about 97% of reads across the
//! four phases, each part a few hundred sequences. The part lists ship inside `PIPseeker`, so the
//! parts are judged by concentration.
//!
//! `PIPseq` V, sold by Illumina as Single Cell 3' RNA Prep, keeps the parts and linkers and
//! replaces the UMI with a 3-base binning index and an ACG spacer before polyT (`DRAGEN`'s
//! documentation and Fontanez et al. 2024, bioRxiv 10.1101/2024.10.04.616561; read off its design
//! run SRR33161101, 9 Oct 2026: the three linkers in 96% of reads, ACG in 93% of those). On
//! Illumina's 45-base Read 1 the spacer shows only in reads with no phasing bases, so such runs
//! are not called: a recall gap, not a wrong call. Its polyT is not required: in the design run
//! only 28% of reads show it cleanly after the spacer (a long homopolymer late in the read), and
//! the ACG spacer already tells it from v3/v4, whose reads hold a random UMI there.

use super::{Layout, Mate, Segment, Tail};
use Segment::{Fixed, Umi, Unlisted};

/// Each part is one of a few hundred sequences: 384 of them must cover most reads.
const PART_SET: usize = 384;
const PART_COVER_PCT_MIN: u8 = 50;

pub const V3_V4: Layout = Layout {
    name: "Illumina Single Cell 3' / PIPseq v3-v4",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        Unlisted { len: 8, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"ATG"]),
        Unlisted { len: 6, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"GAG"]),
        Unlisted { len: 6, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"TCGAG"]),
        Unlisted { len: 8, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Umi(12),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{x[0-3]b[8]f[ATG]b[6]f[GAG]b[6]f[TCGAG]b[8]u[12]x:}2{r:}",
    note: "cell = the four parts; PIPseeker (or DRAGEN for Illumina's kit) maps them to its \
           barcode list; R1 opens with 0-3 phasing bases",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const V5: Layout = Layout {
    name: "Illumina Single Cell 3' RNA Prep / PIPseq V",
    read: Mate::R1,
    offsets: &[0, 1, 2, 3],
    segments: &[
        Unlisted { len: 8, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"ATG"]),
        Unlisted { len: 6, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"GAG"]),
        Unlisted { len: 6, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Fixed(&[b"TCGAG"]),
        Unlisted { len: 8, set: PART_SET, cover_pct_min: PART_COVER_PCT_MIN },
        Umi(3),
        Fixed(&[b"ACG"]),
    ],
    geometry: "1{x[0-3]b[8]f[ATG]b[6]f[GAG]b[6]f[TCGAG]b[8]u[3]f[ACG]x:}2{r:}",
    note: "cell = the four parts; DRAGEN (--scrna-umi-position 39_41 for the binning index) or \
           PIPseeker processes it; no UMI: molecules are counted by R2 start within 64 bins",
    long_read_kit: None,
    cellranger_chemistry: None,
};
