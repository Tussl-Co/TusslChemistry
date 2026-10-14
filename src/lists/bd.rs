// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! BD Rhapsody's cell label lists: 96 per CLS for the original V1 and Enhanced beads (97 in the
//! mirrored files), and 384 per CLS for Enhanced V2/V3 beads, read from BD's own cell-label file
//! (`rhapsody_cell_label.py.txt`, BD's public bucket), whose tuples `B384_cell_key1`-3 hold them.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const BD_TERMS: &str = "BD's list, under BD's terms; fetched, never redistributed";

pub const CLS1: ListEntry = ListEntry {
    id: ListId::BdCls1,
    file: "bd_cls1.txt",
    label: "BD CLS1",
    length: 9,
    count: 97,
    sha256: "71bc753f8bcbbe7e46e26a55f93df8da23efb7ccd6a6c7033e26d082e28fe0bc",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "BD/BD_CLS1.txt"),
    fallback: None,
    terms: BD_TERMS,
};

pub const CLS2: ListEntry = ListEntry {
    id: ListId::BdCls2,
    file: "bd_cls2.txt",
    label: "BD CLS2",
    length: 9,
    count: 97,
    sha256: "e97b234f42aded8a435c629823971c8cd3e29147132bc70dd5d316654075cef5",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "BD/BD_CLS2.txt"),
    fallback: None,
    terms: BD_TERMS,
};

pub const CLS3: ListEntry = ListEntry {
    id: ListId::BdCls3,
    file: "bd_cls3.txt",
    label: "BD CLS3",
    length: 9,
    count: 97,
    sha256: "2a27c1093694a5d58cedec670e4f1c4ef97da6a54793913187c8eb2cc072386b",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "BD/BD_CLS3.txt"),
    fallback: None,
    terms: BD_TERMS,
};

pub const V2_CLS1: ListEntry = ListEntry {
    id: ListId::BdV2Cls1,
    file: "bd_v2_cls1.txt",
    label: "BD Enhanced V2/V3 CLS1",
    length: 9,
    count: 384,
    sha256: "d23b8d87221530ba282a888345e3044cc0f33907a44e22f274f3802fbc2a60b3",
    recipe: Recipe::Quoted(Some("B384_cell_key1")),
    source: (Host::BdPublic, "rhapsody_cell_label.py.txt"),
    fallback: None,
    terms: BD_TERMS,
};

pub const V2_CLS2: ListEntry = ListEntry {
    id: ListId::BdV2Cls2,
    file: "bd_v2_cls2.txt",
    label: "BD Enhanced V2/V3 CLS2",
    length: 9,
    count: 384,
    sha256: "9569286961eefeb0227affe598f4f054963978e5462b9296b79237ebd231fc29",
    recipe: Recipe::Quoted(Some("B384_cell_key2")),
    source: (Host::BdPublic, "rhapsody_cell_label.py.txt"),
    fallback: None,
    terms: BD_TERMS,
};

pub const V2_CLS3: ListEntry = ListEntry {
    id: ListId::BdV2Cls3,
    file: "bd_v2_cls3.txt",
    label: "BD Enhanced V2/V3 CLS3",
    length: 9,
    count: 384,
    sha256: "d33fdcdeb4f4b404d20eb979ce294dfc615d72a00ef5b6d8ffa48643a4b0d8a1",
    recipe: Recipe::Quoted(Some("B384_cell_key3")),
    source: (Host::BdPublic, "rhapsody_cell_label.py.txt"),
    fallback: None,
    terms: BD_TERMS,
};
