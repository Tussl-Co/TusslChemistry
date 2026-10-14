// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Any bytes, split into two inputs at the first 0xFF byte, read as R1 and R2 and run through the
//! whole of detection and every output: the text report, JSON and each pipeline's command. Reads
//! that parse can be anything, so detection must report, never panic. The barcode lists are the
//! real ones, loaded once from TCHEM_FUZZ_WHITELISTS.

#![no_main]

use libfuzzer_sys::fuzz_target;
use std::path::Path;
use std::sync::OnceLock;
use tussl_chemistry::{Pipeline, Reads, Whitelists, detect};

static LISTS: OnceLock<Whitelists> = OnceLock::new();

fn lists() -> &'static Whitelists {
    LISTS.get_or_init(|| {
        let dir = std::env::var("TCHEM_FUZZ_WHITELISTS").expect("set TCHEM_FUZZ_WHITELISTS");
        Whitelists::load(Path::new(&dir)).expect("the barcode lists load")
    })
}

fuzz_target!(|data: &[u8]| {
    let (a, b) = match data.iter().position(|&x| x == 0xFF) {
        Some(i) => (&data[..i], Some(&data[i + 1..])),
        None => (data, None),
    };
    let Ok(r1) = Reads::from_reader(a, Path::new("R1"), 2000) else { return };
    let r2 = match b {
        Some(b) => match Reads::from_reader(b, Path::new("R2"), 2000) {
            Ok(r) => Some(r),
            Err(_) => return,
        },
        None => None,
    };
    let Ok(report) = detect(&r1, r2.as_ref(), lists()) else { return };
    let _ = report.to_text(true);
    let r2_path = r2.as_ref().map(|_| Path::new("R2"));
    let _ = report.to_json(Path::new("R1"), r2_path);
    for pipeline in Pipeline::ALL {
        let _ = report.emit(pipeline, Path::new("R1"), r2_path, true);
    }
});
