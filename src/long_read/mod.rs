// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Long-read single-cell libraries: 10x cDNA read whole on Oxford Nanopore or PacBio. One read
//! holds the adapter, barcode, UMI, polyT or TSO and the transcript, in either orientation, and
//! indel errors shift every position, so the layout is found by searching, not by position.

pub mod detect;
pub mod search;

use crate::constants::{LONG_READ_SHARE_MIN, LONG_READS_EXAMINED, SHORT_READ_BASES_MAX};
use crate::detect::profile::Profile;
use crate::detect::profile::ratio;
use crate::detect::score::Lists;
use crate::input::fastq::Reads;
use detect::{LongEvidence, LongVerdict};

/// Everything a long-read detection found: the sample's read lengths, the evidence, the verdict.
pub struct LongReport {
    pub profile: Profile,
    pub evidence: LongEvidence,
    pub verdict: LongVerdict,
}

/// Whether a sample is long-read data: at least `LONG_READ_SHARE_MIN` of its first reads are
/// longer than any short read.
pub fn is_long(reads: &Reads) -> bool {
    let examined = reads.len().min(LONG_READS_EXAMINED);
    let long = reads.iter().take(LONG_READS_EXAMINED).filter(|r| r.len() > SHORT_READ_BASES_MAX);
    let long = long.count();
    assert!(long <= examined, "long reads are a subset of the reads examined");
    assert!(examined <= LONG_READS_EXAMINED, "routing examines a bounded sample");
    examined > 0 && ratio(long, examined) >= LONG_READ_SHARE_MIN
}

/// Detects the 10x kit behind a sample of long reads.
pub fn detect(reads: &Reads, profile: Profile, lists: &Lists) -> LongReport {
    assert!(is_long(reads), "the long-read detector is for long reads");
    assert_eq!(profile.reads, reads.len(), "the profile describes these reads");
    let evidence = detect::evidence(reads, lists);
    let verdict = detect::decide(&evidence, &detect::kits());
    LongReport { profile, evidence, verdict }
}
