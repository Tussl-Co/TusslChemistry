// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! ScaleBio `QuantumScale`'s RT barcodes, from ScaleBio's open RNA pipeline (`ScaleRna`), pinned
//! to a commit. ScaleBio's licence permits use only with ScaleBio products, so the list is fetched
//! by the user and never shipped.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

pub const QUANTUM_RT: ListEntry = ListEntry {
    id: ListId::ScaleQuantumRt,
    file: "scale_quantum_rt.txt",
    label: "ScaleBio QuantumScale RT barcodes",
    length: 8,
    count: 96,
    sha256: "3eb70409468fdb8f80d8bb3d107d55bae36bc59fec9a63af165ef795c8913106",
    recipe: Recipe::TsvColumn2,
    source: (Host::ScaleRna, "references/quantum_rt.txt"),
    fallback: None,
    terms: "ScaleBio's list, under ScaleBio's licence (use with ScaleBio products); \
            fetched, never redistributed",
};
