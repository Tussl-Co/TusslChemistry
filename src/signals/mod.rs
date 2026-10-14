// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Evidence from the cDNA read that the barcode layout cannot give: which strand it is on (3' or
//! 5' kits), what kind of library it is, and, when no layout fits, what structure the reads show.
//! All of it is reference-light: k-mers of a few built-in sequences, no genome index.

pub mod instrument;
pub mod kmers;
pub mod library;
pub mod strand;
pub mod structure;
