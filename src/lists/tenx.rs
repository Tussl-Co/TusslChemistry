// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! 10x Genomics' barcode lists. 10x's lists are 10x's, under 10x's terms: fetched, never
//! shipped. The Flex list ships only inside Cell Ranger.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const TENX_TERMS: &str = "10x Genomics' list, under 10x's terms; fetched, never redistributed";

pub const V2: ListEntry = ListEntry {
    id: ListId::TenxV2,
    file: "10x_v2_permit.txt",
    label: "10x 737K-august-2016",
    length: 16,
    count: 737_280,
    sha256: "b0dda4b114d8fc7ea8def94fb472e6205c1831a72e3ac245ed13cc332d9f20d1",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "10X-Genomics/737K-august-2016.txt.gz"),
    fallback: Some((Host::SimpleafBox, "jbs2wszgbj7k4ic2hass9ts6nhqkwq1p")),
    terms: TENX_TERMS,
};

pub const V3: ListEntry = ListEntry {
    id: ListId::TenxV3,
    file: "10x_v3_permit.txt",
    label: "10x 3M-february-2018",
    length: 16,
    count: 6_794_880,
    sha256: "843a6f7038db8cb3c06f3dc21cc69d04139ffa10689780518d7a9e42dc2e819b",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "10X-Genomics/3M-february-2018.txt.gz"),
    fallback: Some((Host::SimpleafBox, "vc9zd4qyjj581gvtolw5kj638wmg4f3s")),
    terms: TENX_TERMS,
};

pub const V4_THREE_PRIME: ListEntry = ListEntry {
    id: ListId::TenxV4ThreePrime,
    file: "10x_v4_3p_permit.txt",
    label: "10x 3M-3pgex-may-2023",
    length: 16,
    count: 7_372_800,
    sha256: "62df04844ed0d1281b302f862ec0d005d3298fc169f507f68dd6e8d4cb5ccccd",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "10X-Genomics/3M-3pgex-may-2023.txt.gz"),
    fallback: Some((Host::SimpleafBox, "4a2q2lxngogaq4vkniptpu2gts5a38fb.txt")),
    terms: TENX_TERMS,
};

pub const V3_FIVE_PRIME: ListEntry = ListEntry {
    id: ListId::TenxV3FivePrime,
    file: "10x_v3_5p_permit.txt",
    label: "10x 3M-5pgex-jan-2023",
    length: 16,
    count: 3_686_400,
    sha256: "b402db8acbcd3a076bf50991538d6c08464ddb137e0b81916a9fa92cf222dd0d",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "10X-Genomics/3M-5pgex-jan-2023.txt.gz"),
    fallback: Some((Host::SimpleafBox, "cbpv1c4zi6ty81nvcgy3pyta2oj7vea1.txt")),
    terms: TENX_TERMS,
};

pub const ARC_GEX: ListEntry = ListEntry {
    id: ListId::TenxArcGex,
    file: "10x_arc_v1_gex_permit.txt",
    label: "10x 737K-arc-v1 (gene expression)",
    length: 16,
    count: 736_320,
    sha256: "23ea010e4854aa5097f213719a86479e47ff3b33b6541478cbe0c84c6e3818b5",
    recipe: Recipe::Lines,
    source: (Host::Teichlab, "10X-Genomics/gex_737K-arc-v1.txt.gz"),
    fallback: None,
    terms: TENX_TERMS,
};

pub const FLEX_V1: ListEntry = ListEntry {
    id: ListId::TenxFlexV1,
    file: "10x_flex_v1_permit.txt",
    label: "10x 737K-fixed-rna-profiling",
    length: 16,
    count: 737_280,
    sha256: "74c2c9d186f4ffb3ea86fcb35d0f7a0de91bc0aee3b88e35d84bad89bbf7f25f",
    recipe: Recipe::Lines,
    source: (Host::CellRanger, "737K-fixed-rna-profiling.txt.gz"),
    fallback: None,
    terms: TENX_TERMS,
};
