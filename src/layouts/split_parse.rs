// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Split-pool combinatorial barcoding: SPLiT-seq (Rosenberg et al. 2018) and Parse Evercode WT,
//! in its original layout (WT v1, v2 and early v3 kits) and the short-linker v3 layout.

use super::ListId;
use super::bc;
use super::{Layout, Mate, Segment, Tail};
use Segment::{Fixed, Umi, Unlisted};

/// SPLiT-seq / Parse linker between the round 3 and round 2 barcodes.
const SPLIT_L1: &[u8] = b"GTGGCCGATGTTTCGCATCGGCGTACGACT";
/// SPLiT-seq (Rosenberg et al. 2018) linker between round 2 and round 1.
const SPLIT_L2: &[u8] = b"ATCCACGTGCTTGAGAGGCCAGAGCATTCG";
/// Parse Evercode's shorter round 2 to round 1 linker, read off public Parse runs (7 Oct 2026):
/// shares SPLiT-seq's first 15 bases, then diverges; R2 is 86 bases instead of 94.
const PARSE_L2: &[u8] = b"ATCCACGTGCTTGAGACTGTGG";
/// Parse Evercode WT v3's short linkers, read off public v3 runs (8 Oct 2026): round 3 to round 2,
/// and round 2 to round 1. The first's GGGG often reads badly on two-colour sequencers, which see
/// G as no signal; the linker tolerance (two mismatches in twelve bases) absorbs it.
const PARSE_V3_L1: &[u8] = b"ATGAGGGGTCAG";
const PARSE_V3_L2: &[u8] = b"TCCAACCACCTC";

pub const SPLIT_SEQ: Layout = Layout {
    name: "SPLiT-seq (Rosenberg 2018)",
    read: Mate::R2,
    offsets: &[0],
    segments: &[
        Umi(10),
        bc(ListId::SplitRound3, 8),
        Fixed(&[SPLIT_L1]),
        bc(ListId::SplitRound2, 8),
        Fixed(&[SPLIT_L2]),
        bc(ListId::SplitRound1, 8),
    ],
    geometry: concat!(
        "1{r:}2{u[10]b[8]f[GTGGCCGATGTTTCGCATCGGCGTACGACT]b[8]",
        "f[ATCCACGTGCTTGAGAGGCCAGAGCATTCG]b[8]x:}",
    ),
    note: "cDNA is R1; round 1 pairs an oligo-dT and a random-hexamer barcode per well, \
           merge them; process each sublibrary (i7) separately",
    long_read_kit: None,
    cellranger_chemistry: None,
};

pub const PARSE_WT: Layout = Layout {
    name: "Parse Evercode WT",
    read: Mate::R2,
    offsets: &[0],
    segments: &[
        Umi(10),
        bc(ListId::SplitRound3, 8),
        Fixed(&[SPLIT_L1]),
        bc(ListId::SplitRound2, 8),
        Fixed(&[PARSE_L2]),
        Unlisted { len: 8, set: 384, cover_pct_min: 50 },
    ],
    geometry: concat!(
        "1{r:}2{u[10]b[8]f[GTGGCCGATGTTTCGCATCGGCGTACGACT]b[8]",
        "f[ATCCACGTGCTTGAGACTGTGG]b[8]x:}",
    ),
    note: "cDNA is R1; round 1 barcodes are Parse's own (not public), merge its \
           oligo-dT/hexamer pairs; process each sublibrary separately",
    long_read_kit: None,
    cellranger_chemistry: None,
};

/// Parse Evercode WT v3 as sold from 2024: shorter linkers, Parse's v3 tables for all three rounds,
/// polyT after round 1 when R2 reads that far (archives often trim R2 to its 58 structured bases).
pub const PARSE_WT_V3: Layout = Layout {
    name: "Parse Evercode WT v3",
    read: Mate::R2,
    offsets: &[0],
    segments: &[
        Umi(10),
        bc(ListId::ParseV3Rounds23, 8),
        Fixed(&[PARSE_V3_L1]),
        bc(ListId::ParseV3Rounds23, 8),
        Fixed(&[PARSE_V3_L2]),
        bc(ListId::ParseV3Round1, 8),
        Segment::Tail(Tail::PolyT),
    ],
    geometry: "1{r:}2{u[10]b[8]f[ATGAGGGGTCAG]b[8]f[TCCAACCACCTC]b[8]x:}",
    note: "cDNA is R1; cell = round 1 + round 2 + round 3; round 1 pairs an oligo-dT and a \
           random-hexamer barcode per well, merge them (Parse's split-pipe does)",
    long_read_kit: None,
    cellranger_chemistry: None,
};
