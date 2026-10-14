// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! `SeekGene` `SeekOne` barcode lists, from `SeekGene`'s `SeekSoul` Tools (MIT licence), pinned to
//! commit a9a508e (the 2023 copy of the 17-base list is corrupt: 959 entries, a byte-order mark,
//! CRLF). The DD kits' 17-base list is exactly every pairing of 960 7-base parts and 1,536
//! 10-base parts (checked 9 Oct 2026). Only the 10-base part is held as a list: allowing one
//! substitution, 960 of the 16,384 7-mers cover nearly all of them, so a 7-base list would match
//! any read (9 Oct 2026: unrelated tag libraries were vouched for as `SeekGene` at LOW); the
//! layout checks the 7-base part by concentration instead. The MM kit's 479-entry list is not
//! registered while its layout is not claimed (`src/layouts/seekgene.rs`). `SeekArcTools`' copy
//! of the 17-base file is byte-identical.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const SEEKGENE_TERMS: &str = "SeekGene's list, MIT licence (SeekSoul Tools)";
const P3CB: &str = "P3CBGB/P3CB.barcode.txt.gz";

pub const P3CB_PART_2: ListEntry = ListEntry {
    id: ListId::SeekgeneP3cb2,
    file: "seekgene_p3cb_part2.txt",
    label: "SeekGene DD barcode, bases 8-17",
    length: 10,
    count: 1536,
    sha256: "93e422e68c1f5eaa1f6b4ede79ba5a1e55769c4d6e43d3581bd462e25e63f0db",
    recipe: Recipe::Slice(7, 10),
    source: (Host::SeekSoul, P3CB),
    fallback: Some((Host::SeekArc, P3CB)),
    terms: SEEKGENE_TERMS,
};
