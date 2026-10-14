// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Singleron GEXSCOPE's barcode lists, from Singleron's `CeleScope` pipeline (MIT licence).

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const SINGLERON_TERMS: &str = "Singleron's list, MIT licence (CeleScope)";

pub const V1: ListEntry = ListEntry {
    id: ListId::SingleronV1,
    file: "singleron_v1_bc.txt",
    label: "Singleron V1 barcodes",
    length: 8,
    count: 96,
    sha256: "346ecc6bf50e72c91e0ec115ca22918abbfa8d3dae7be0dab9e602a1564de199",
    recipe: Recipe::Lines,
    source: (Host::Celescope, "GEXSCOPE-V1/bc.txt"),
    fallback: None,
    terms: SINGLERON_TERMS,
};

pub const BC1: ListEntry = ListEntry {
    id: ListId::SingleronBc1,
    file: "singleron_bc1.txt",
    label: "Singleron barcode 1",
    length: 9,
    count: 96,
    sha256: "d8c09192245ea1e80cddc46cc2235cd9e4b1a936a4db6ab384d097ecc8307d19",
    recipe: Recipe::Lines,
    source: (Host::Celescope, "GEXSCOPE-V2/bc1.txt"),
    fallback: None,
    terms: SINGLERON_TERMS,
};

pub const BC2: ListEntry = ListEntry {
    id: ListId::SingleronBc2,
    file: "singleron_bc2.txt",
    label: "Singleron barcode 2",
    length: 9,
    count: 96,
    sha256: "0529f87f57cc0b9623be4a497fbb97322214b4a44a572942306d42a23f8967f2",
    recipe: Recipe::Lines,
    source: (Host::Celescope, "GEXSCOPE-V2/bc2.txt"),
    fallback: None,
    terms: SINGLERON_TERMS,
};

pub const BC3: ListEntry = ListEntry {
    id: ListId::SingleronBc3,
    file: "singleron_bc3.txt",
    label: "Singleron barcode 3",
    length: 9,
    count: 192,
    sha256: "fb643607c4897464bdf5ea4c181ac5cd1f57c4eee16272477ae1e8a916a4adb2",
    recipe: Recipe::Lines,
    source: (Host::Celescope, "GEXSCOPE-V2/bc3.txt"),
    fallback: None,
    terms: SINGLERON_TERMS,
};
