// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The registry: each list's entry type, the table of every entry in `ListId` order (assembled
//! from the vendor modules), its compile-time checks, and the test that renders `data/lists.tsv`.

use super::{ListId, bd, mgi, parse, scalebio, seekgene, singleron, split_seq, tenx};

/// How the fetch script turns a download into a list of one barcode per line. Either way a gzip
/// download is decompressed first, carriage returns are dropped and empty lines skipped.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Recipe {
    /// One barcode per line.
    Lines,
    /// A CSV with a header row and the barcode in column 2, sorted (byte order) and deduplicated.
    CsvColumn2,
    /// A tab-separated file without a header and the barcode in column 2, sorted and deduplicated.
    TsvColumn2,
    /// Every double-quoted ACGT string in a source file (a vendor's JSON config or code), or with a
    /// name only those in the block that opens with `name =`; sorted and deduplicated.
    Quoted(Option<&'static str>),
    /// The distinct substrings at `[start, start + len)` of every line, sorted and deduplicated:
    /// one part of a list whose barcodes are every pairing of two parts (`SeekGene`'s 17-base list
    /// is 960 7-base parts by 1,536 10-base parts, past the 16 bases a packed list holds).
    Slice(usize, usize),
}

/// Where lists are fetched from. GitHub sources are pinned to a commit, so a list cannot change
/// under the checksum; the checksum still guards every download.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Host {
    /// The Teichmann lab's `scg_lib_structs`. It has no licence of its own: it mirrors lists that
    /// belong to their vendors, whose terms apply.
    Teichlab,
    /// Singleron's `CeleScope` pipeline (MIT licence).
    Celescope,
    /// simpleaf's mirror of 10x's lists (COMBINE-lab/simpleaf, `resources/permit_list_info.json`).
    SimpleafBox,
    /// A Cell Ranger install, which the user downloads from 10x under 10x's licence; the fetch
    /// script copies from it when given `--cellranger <dir>`. The path is under the install.
    CellRanger,
    /// Parse's barcode tables as mirrored in `mortazavilab/parse_pipeline` (IGVF processing; no
    /// licence of its own: the tables are Parse's, under Parse's terms).
    ParseMirror,
    /// MGI's `DNBelab` C-series analysis software on GitHub (MIT licence).
    MgiGithub,
    /// BD's public bucket for its sequence analysis pipeline's cell-label files (no licence stated:
    /// BD's lists, under BD's terms).
    BdPublic,
    /// ScaleBio's open RNA pipeline on GitHub (`ScaleRna`; ScaleBio's licence, use with ScaleBio
    /// products only).
    ScaleRna,
    /// `SeekGene`'s `SeekSoul` Tools pipeline on GitHub (MIT licence).
    SeekSoul,
    /// `SeekGene`'s `SeekArcTools` (MIT licence), whose barcode files are byte-identical copies.
    SeekArc,
}

impl Host {
    #[cfg_attr(not(test), expect(dead_code, reason = "read where data/lists.tsv is rendered"))]
    pub const fn url(self) -> &'static str {
        match self {
            Host::Teichlab => {
                "https://raw.githubusercontent.com/Teichlab/scg_lib_structs/\
                 3b61a870094889ece95a4daee3314b5186e849a5/data"
            }
            Host::Celescope => {
                "https://raw.githubusercontent.com/singleron-RD/CeleScope/\
                 851c7b891a747a31d4271f51571b09d02a64f633/celescope/data/chemistry"
            }
            Host::SimpleafBox => "https://umd.box.com/shared/static",
            Host::CellRanger => "cellranger:lib/python/cellranger/barcodes",
            Host::BdPublic => "https://bd-rhapsody-public.s3.amazonaws.com/CellLabel",
            Host::ScaleRna => {
                "https://raw.githubusercontent.com/ScaleBio/ScaleRna/\
                 32fac729093130ab70c47a987c10daf55063b99b"
            }
            Host::MgiGithub => {
                "https://raw.githubusercontent.com/MGI-tech-bioinformatics/\
                 DNBelab_C_Series_scRNA-analysis-software/2b5c226c189ae24c1dc3ebab242e8cdd65a112f9"
            }
            Host::SeekSoul => {
                "https://raw.githubusercontent.com/seekgenebio/seeksoultools/\
                 a9a508ea75e9dc41d7821f47f7475ffbfc0ab74a/src/seeksoultools/utils/barcode"
            }
            Host::SeekArc => {
                "https://raw.githubusercontent.com/seekgene/SeekArcTools/\
                 a49e76b0f443b1fd9568381ee56806d19873f51f/src/seekarctools/utils/barcode"
            }
            Host::ParseMirror => {
                "https://raw.githubusercontent.com/mortazavilab/parse_pipeline/\
                 b156bbd63fb21acfb8e1c65c51f14e0ab49a6b78/barcodes"
            }
        }
    }
}

pub struct ListEntry {
    pub id: ListId,
    pub file: &'static str,
    pub label: &'static str,
    pub length: usize,
    /// Barcodes in the list file.
    pub count: usize,
    /// sha256 of the list file exactly as the fetch script writes it.
    pub sha256: &'static str,
    #[cfg_attr(not(test), expect(dead_code, reason = "read where data/lists.tsv is rendered"))]
    pub recipe: Recipe,
    /// The host and the path under it.
    pub source: (Host, &'static str),
    /// A second source, tried when the first fails or its download does not check.
    pub fallback: Option<(Host, &'static str)>,
    /// Who the list belongs to, and on what terms it may be used.
    #[cfg_attr(not(test), expect(dead_code, reason = "read where data/lists.tsv is rendered"))]
    pub terms: &'static str,
}

/// Every list, in `ListId` order. Counts and checksums are of lists fetched on 8 Oct 2026; for
/// 10x's lists both sources gave byte-identical files, and so did Cell Ranger 10.1.0's copy of
/// 737K-arc-v1. The Flex list's checksum is of Cell Ranger 10.1.0's copy. SPLiT-seq's three round
/// files are identical upstream: the same 96 sequences serve all three rounds.
pub const REGISTRY: [ListEntry; 24] = [
    tenx::V2,
    tenx::V3,
    tenx::V4_THREE_PRIME,
    tenx::V3_FIVE_PRIME,
    bd::CLS1,
    bd::CLS2,
    bd::CLS3,
    split_seq::ROUND1,
    split_seq::ROUND2,
    split_seq::ROUND3,
    singleron::V1,
    singleron::BC1,
    singleron::BC2,
    singleron::BC3,
    tenx::ARC_GEX,
    tenx::FLEX_V1,
    parse::V3_ROUNDS_2_3,
    parse::V3_ROUND_1,
    mgi::C_SERIES,
    bd::V2_CLS1,
    bd::V2_CLS2,
    bd::V2_CLS3,
    scalebio::QUANTUM_RT,
    seekgene::P3CB_PART_2,
];

/// sha256 digests are written as 64 lower-case hex digits.
const SHA256_HEX_DIGITS: usize = 64;

/// Compile-time checks on the registry, one property per assertion.
const fn registry_is_sound() -> bool {
    let mut i = 0;
    while i < REGISTRY.len() {
        #[allow(clippy::indexing_slicing)] // earned by the loop bound
        let entry = &REGISTRY[i];
        assert!(entry.id.index() == i, "the registry is in list id order");
        assert!(entry.length >= 1, "a barcode has at least one base");
        assert!(entry.length <= crate::input::whitelist::LEN_MAX, "a barcode packs into a u32");
        assert!(entry.count >= 1, "a list has at least one barcode");
        assert!(!entry.file.is_empty(), "a list has a file name");
        assert!(entry.sha256.len() == SHA256_HEX_DIGITS, "a sha256 is 64 hex digits");
        let digits = entry.sha256.as_bytes();
        let mut d = 0;
        while d < digits.len() {
            #[allow(clippy::indexing_slicing)] // earned by the loop bound
            let digit = digits[d];
            assert!(digit.is_ascii_digit() || (digit >= b'a' && digit <= b'f'), "lower-case hex");
            d += 1;
        }
        i += 1;
    }
    assert!(i == ListId::ALL.len(), "one entry per list id");
    true
}

const _: () = assert!(registry_is_sound(), "the list registry is unsound");

#[cfg(test)]
mod tests {
    use super::*;

    const REGISTRY_FILE: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/data/lists.tsv");

    fn url((host, path): (Host, &str)) -> String {
        format!("{}/{path}", host.url())
    }

    /// The registry as the fetch script reads it: a header, then one tab-separated row per list.
    fn registry_tsv() -> String {
        use std::fmt::Write;
        let mut out = String::from(
            "# Rendered from src/lists/registry.rs by its tests (TCHEM_BLESS=1 cargo test \
             registry); do not edit.\n\
             file\tlength\tcount\tsha256\trecipe\tsource\tfallback\tterms\n",
        );
        for entry in &REGISTRY {
            let recipe = match entry.recipe {
                Recipe::Lines => "lines".to_string(),
                Recipe::CsvColumn2 => "csv-column-2".to_string(),
                Recipe::TsvColumn2 => "tsv-column-2".to_string(),
                Recipe::Quoted(None) => "quoted".to_string(),
                Recipe::Quoted(Some(name)) => format!("quoted:{name}"),
                Recipe::Slice(start, len) => format!("slice:{start}:{len}"),
            };
            // "-", not an empty field: shell `read` merges consecutive tabs.
            let fallback = entry.fallback.map_or_else(|| "-".to_string(), url);
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{recipe}\t{}\t{fallback}\t{}",
                entry.file,
                entry.length,
                entry.count,
                entry.sha256,
                url(entry.source),
                entry.terms
            )
            .unwrap();
        }
        out
    }

    #[test]
    fn the_registry_is_sound() {
        assert!(registry_is_sound());
        for id in ListId::ALL {
            assert_eq!(id.entry().id, id);
        }
    }

    /// `data/lists.tsv` is the registry rendered; regenerate with `TCHEM_BLESS=1`.
    #[test]
    fn the_registry_file_is_current() {
        let want = registry_tsv();
        if std::env::var_os("TCHEM_BLESS").is_some() {
            std::fs::write(REGISTRY_FILE, &want).unwrap();
        }
        let have = std::fs::read_to_string(REGISTRY_FILE).unwrap_or_default();
        assert!(have == want, "data/lists.tsv is stale; rerun with TCHEM_BLESS=1");
    }
}
