// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! What the tool reads: the first reads of each FASTQ or BAM file, and the vendors' barcode lists.
//! Both are untrusted, so every fault is a typed refusal, never a panic.

pub mod bam;
pub mod fastq;
pub mod whitelist;
