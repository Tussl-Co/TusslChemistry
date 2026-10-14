// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! MGI `DNBelab` C-series cell barcodes: one 1,536-barcode table read at both barcode positions,
//! from MGI's own read-structure config (MIT licence), pinned to a commit. The high-throughput
//! kit reads the same positions with the same table (checked 8 Oct 2026 on public runs).

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

pub const C_SERIES: ListEntry = ListEntry {
    id: ListId::MgiCSeries,
    file: "mgi_c_series_bc.txt",
    label: "MGI DNBelab C-series barcodes",
    length: 10,
    count: 1536,
    sha256: "9b67b81a65056a8b1c6cecccd8d388cf214bf87942df59466530b0b82aa30315",
    recipe: Recipe::Quoted(None),
    source: (Host::MgiGithub, "config/DNBelabC4_scRNA_readStructure.json"),
    fallback: None,
    terms: "MGI's list, MIT licence (DNBelab C Series scRNA analysis software)",
};
