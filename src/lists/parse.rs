// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Parse Biosciences Evercode WT v3 barcode tables. Parse ships them inside its own pipeline
//! (split-pipe, behind a login); they are fetched from a public mirror in an IGVF lab's pipeline
//! repository, pinned to a commit. Rounds 2 and 3 share one 96-barcode table; round 1 (the RT
//! barcode, oligo-dT and random-hexamer wells) is the largest v3 table, 208 barcodes, which holds
//! the smaller kits' tables (26 and 102). Checked 8 Oct 2026 on public v3 runs: rounds 2/3 hit
//! 94-98%, round 1 93-96%.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const PARSE_TERMS: &str =
    "Parse Biosciences' list, under Parse's terms; fetched, never redistributed";

pub const V3_ROUNDS_2_3: ListEntry = ListEntry {
    id: ListId::ParseV3Rounds23,
    file: "parse_v3_rounds23.txt",
    label: "Parse WT v3 rounds 2-3",
    length: 8,
    count: 96,
    sha256: "42422e1c25b76ec7f934f7e6959f136cf5032689632aa3848b09911ad0fdd81f",
    recipe: Recipe::CsvColumn2,
    source: (Host::ParseMirror, "bc_data_R3_v3.csv"),
    fallback: None,
    terms: PARSE_TERMS,
};

pub const V3_ROUND_1: ListEntry = ListEntry {
    id: ListId::ParseV3Round1,
    file: "parse_v3_round1.txt",
    label: "Parse WT v3 round 1",
    length: 8,
    count: 208,
    sha256: "b1e3e31671b6326c7240a6b4a45c03e5f7745d7438db6aeab079781c21a46cf0",
    recipe: Recipe::CsvColumn2,
    source: (Host::ParseMirror, "bc_data_n208_R1_v3_3.csv"),
    fallback: None,
    terms: PARSE_TERMS,
};
