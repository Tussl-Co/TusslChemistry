// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Any bytes, as a file or stream tchem is given: FASTQ (plain, gzip, BGZF), BAM or an SRA archive,
//! sniffed from the first bytes. A malformed input must be a typed refusal, never a panic.

#![no_main]

use libfuzzer_sys::fuzz_target;
use std::path::Path;
use tussl_chemistry::Reads;

fuzz_target!(|data: &[u8]| {
    let _ = Reads::from_reader(data, Path::new("fuzz"), 1000);
});
