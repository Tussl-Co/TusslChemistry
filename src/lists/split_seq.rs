// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! The SPLiT-seq round barcodes (Rosenberg et al. 2018, Science, table S12). The three round
//! files are identical upstream: the same 96 sequences serve all three rounds.

use super::ListId;
use super::registry::{Host, ListEntry, Recipe};

const SPLIT_TERMS: &str = "published by Rosenberg et al. 2018, Science, table S12";

pub const ROUND1: ListEntry = ListEntry {
    id: ListId::SplitRound1,
    file: "splitseq_round1.txt",
    label: "SPLiT-seq round 1",
    length: 8,
    count: 96,
    sha256: "8c100ab73ceb375eeeca5994cf811899a6c058bf653388f291500f0464ef5096",
    recipe: Recipe::CsvColumn2,
    source: (Host::Teichlab, "SPLiT-seq/SPLiT-seq_Round1_bc.csv"),
    fallback: None,
    terms: SPLIT_TERMS,
};

pub const ROUND2: ListEntry = ListEntry {
    id: ListId::SplitRound2,
    file: "splitseq_round2.txt",
    label: "SPLiT-seq round 2",
    length: 8,
    count: 96,
    sha256: "8c100ab73ceb375eeeca5994cf811899a6c058bf653388f291500f0464ef5096",
    recipe: Recipe::CsvColumn2,
    source: (Host::Teichlab, "SPLiT-seq/SPLiT-seq_Round2_bc.csv"),
    fallback: None,
    terms: SPLIT_TERMS,
};

pub const ROUND3: ListEntry = ListEntry {
    id: ListId::SplitRound3,
    file: "splitseq_round3.txt",
    label: "SPLiT-seq round 3",
    length: 8,
    count: 96,
    sha256: "8c100ab73ceb375eeeca5994cf811899a6c058bf653388f291500f0464ef5096",
    recipe: Recipe::CsvColumn2,
    source: (Host::Teichlab, "SPLiT-seq/SPLiT-seq_Round3_bc.csv"),
    fallback: None,
    terms: SPLIT_TERMS,
};
