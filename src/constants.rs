// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Every tunable number in the tool, in one place: the limits on what it reads, and every
//! threshold behind a judgement, each with where its value comes from. Retuning the tool means
//! editing this file; the parity oracle then shows exactly which verdicts moved.
//!
//! Facts about formats and data stay next to the code that uses them: the 2-bit base code, the
//! packed whitelist header, the JSON number layout, the kits' linker and tail sequences
//! (`layouts/`), the adapter sequences (`signals/structure.rs`) and the table bounds the layout
//! checks enforce. Those are not choices, so they are not here.

// ---------------------------------------------------------------------------------------------
// What the tool reads: input limits.
// ---------------------------------------------------------------------------------------------

/// Most reads one sample holds. A million covers every layout check with room to spare, and bounds
/// the sample's memory alongside `SAMPLE_BASES_MAX`.
pub const SAMPLE_READS_MAX: usize = 1_000_000;
/// A sensible sample size: enough reads for every check, about a second of work.
pub const SAMPLE_READS_DEFAULT: usize = 100_000;
/// Most bases kept across all sampled reads (512 MiB). 100k short reads use under 30 MB; the cap
/// only bites on long reads, where sampling stops early rather than exhausting memory.
pub const SAMPLE_BASES_MAX: usize = 1 << 29;
/// Longest FASTQ line accepted, in bytes. Long-read platforms reach tens of kilobases per read; a
/// line past a mebibyte is not a read this tool can use, and the cap bounds memory per line.
pub const FASTQ_LINE_BYTES_MAX: usize = 1 << 20;
/// Most barcodes one list may hold. The largest 10x list (3M-3pgex-may-2023) has 7,372,800; twice
/// that bounds memory at about 120 MB per list while leaving room for future kits.
pub const BARCODES_MAX: usize = 1 << 24;
/// Longest whitelist line accepted: a barcode, optional surrounding whitespace, and the newline.
pub const WHITELIST_LINE_BYTES_MAX: usize = 64;

const _: () = assert!(SAMPLE_READS_DEFAULT >= 1);
const _: () = assert!(SAMPLE_READS_DEFAULT <= SAMPLE_READS_MAX);
const _: () = assert!(FASTQ_LINE_BYTES_MAX < SAMPLE_BASES_MAX, "a sample holds at least one line");

// ---------------------------------------------------------------------------------------------
// Scoring a layout against the reads.
// ---------------------------------------------------------------------------------------------

/// Positions profiled per read: the base-composition sketch and the fixed-column count.
pub const PROFILE_POSITIONS: usize = 96;
/// A UMI position below this entropy (bits) is not random enough to be a UMI.
pub const UMI_ENTROPY_MIN: f64 = 1.5;
/// Share of reads that must show a linker or the expected tail for it to count as present.
pub const PRESENT_MIN: f64 = 0.5;
/// Share of reads a listed barcode segment must hit.
pub const SEGMENT_HITS_MIN: f64 = 0.3;
/// A listed barcode hit rate this high vouches for a layout even when something else is off.
pub const VOUCHING_HITS: f64 = 0.8;
/// A list's hits come from many cells, so no one sequence holds this share of them. A constant
/// sequence (a tag or amplicon read, one adapter) puts nearly all its hits on one: a BD tag
/// library's amplicon read sat one substitution from a `SeekGene` 10-mer and hit that list in 97%
/// of reads, which vouched for `SeekGene` DD and won the orientation from BD (coverage at
/// c157518, 9 Oct 2026). Real libraries spread their hits over thousands of barcodes, so the bar
/// has a wide margin on both sides. A layout's lists vouch for it, or anchor HIGH, only when at
/// least one of its listed segments is spread under this bar; a segment that is legitimately
/// narrow (a round-1 barcode from a few wells) keeps its rate.
pub const LIST_TOP_SHARE_MAX: f64 = 0.5;
/// Linker tolerance: one mismatch allowed per this many bases (a 4-base linker must match
/// exactly).
pub const LINKER_BASES_PER_MISMATCH: usize = 5;
/// Tail tolerance: one mismatch allowed per this many tail bases shown.
pub const TAIL_BASES_PER_MISMATCH: usize = 4;

const _: () =
    assert!(SEGMENT_HITS_MIN < VOUCHING_HITS, "vouching needs more than a passing hit rate");
const _: () = assert!(PRESENT_MIN > 0.0 && PRESENT_MIN < 1.0);
const _: () = assert!(LIST_TOP_SHARE_MAX > 0.0 && LIST_TOP_SHARE_MAX < 1.0);

// ---------------------------------------------------------------------------------------------
// Deciding between layouts.
// ---------------------------------------------------------------------------------------------

/// Below this support no layout is claimed.
pub const SUPPORT_MIN: f64 = 0.30;
/// HIGH needs at least this support, however the layout won: a mixture of two libraries supports
/// each at about half and must not be called HIGH even when a tie-break is decisive.
pub const SUPPORT_HIGH: f64 = 0.7;
/// HIGH for a layout with no barcode list needs each of its linkers in at least this share of
/// reads. Support scores barcode concentration as pass/fail against each kit's own bar, so a
/// 50/50 mix of two libraries passed those segments in full and reached HIGH (chaos: ScaleBio RNA
/// mixed with 10x reads, 9 Oct 2026). A linker is in a read or not, so its rate is the share of
/// reads from this kit: 92-100% in clean design samples (`PIPseq` v3/v4, ScaleBio RNA, inDrops),
/// 50% in the mixes. Concentration is not used: a real inDrops run covers its second barcode in
/// 67% of reads (wild ERR10146550) with its linker in 93%. A listed layout's hit rate already
/// counts every read, so it needs no such bar.
pub const UNLISTED_LINKER_HIGH: f64 = 0.7;
/// Layouts whose support differs by less than this are tied (the reads cannot separate them).
pub const SUPPORT_TIE_TOLERANCE: f64 = 0.02;
/// A clean rival holding more than this share of the top layout's support denies HIGH.
pub const RIVAL_SUPPORT_SHARE: f64 = 0.5;
/// Tail bases a read must show for the tail to break a 3'/5' tie (a polyT and a TSO differ from
/// the first base, and four bases make a chance match unlikely).
pub const TAIL_SEEN_BASES_MIN: usize = 4;

const _: () = assert!(SUPPORT_MIN < SUPPORT_HIGH, "HIGH asks more than claiming at all");

// ---------------------------------------------------------------------------------------------
// The cDNA read's strand (3' or 5' kits), without a reference genome.
// ---------------------------------------------------------------------------------------------

/// A homopolymer run this long marks a polyA (sense) or polyT (antisense) tail in the cDNA read.
pub const HOMOPOLYMER_RUN: usize = 15;
/// k-mer length and window step of every reference set in the shared table: 25-mers sampled every
/// 5 bases.
pub const KMER_LENGTH: usize = 25;
pub const KMER_STEP: usize = 5;
/// A read counts for a reference set when at least this many of its windows are in it.
pub const KMER_WINDOWS_MIN: usize = 2;
/// A k-mer with a run of this many identical bases is low-complexity (A-rich transcript ends
/// match polyT stretches on the other strand), so it never votes.
pub const LOW_COMPLEXITY_RUN: usize = 8;
/// A strand signal leans when its leading side holds at least this share of reads and
/// `STRAND_RATIO_MIN` times the other side; it is strong at `STRAND_STRONG_*`. Calibrated on
/// known 10x data (3' reads 5-7x sense, 5' reads 24x antisense).
pub const STRAND_SHARE_MIN: f64 = 0.005;
pub const STRAND_RATIO_MIN: f64 = 3.0;
pub const STRAND_STRONG_SHARE_MIN: f64 = 0.01;
pub const STRAND_STRONG_RATIO_MIN: f64 = 10.0;

const _: () = assert!(KMER_LENGTH <= 32, "a k-mer packs into a u64");
const _: () = assert!(KMER_STEP >= 1 && KMER_STEP <= KMER_LENGTH, "windows overlap");
const _: () = assert!(LOW_COMPLEXITY_RUN <= KMER_LENGTH);
const _: () = assert!(STRAND_SHARE_MIN <= STRAND_STRONG_SHARE_MIN);
const _: () = assert!(STRAND_RATIO_MIN <= STRAND_STRONG_RATIO_MIN);

// ---------------------------------------------------------------------------------------------
// Library type. Calibrated on 7 Oct 2026 over 480 public scBaseCount samples
// (`bench/library_type`), library type from ENA titles: VDJ median 20% constant-region reads at
// 0.9% housekeeping; gene expression 0.01% (95th percentile 2.7%), housekeeping median 11%; guide
// libraries >= 24% scaffold reads against 0.3% (95th percentile) for gene expression; tag
// libraries >= 80% top-200 start cover.
// ---------------------------------------------------------------------------------------------

/// Reads examined; enough for shares in the low percent.
pub const LIBRARY_READS: usize = 50_000;
/// VDJ: constant-region reads at least this share, and at least `VDJ_OVER_HOUSEKEEPING` times the
/// housekeeping share (3x: VDJ libraries sit far above it, gene expression near 1x).
pub const VDJ_SHARE_MIN: f64 = 0.05;
pub const VDJ_OVER_HOUSEKEEPING: f64 = 3.0;
/// Guide: sgRNA scaffold in at least this share of reads (gene expression from Perturb-seq, whose
/// guide transcripts are polyadenylated, reaches 6%).
pub const GUIDE_SHARE_MIN: f64 = 0.15;
/// Scaffold k-mers: 12-mers every 3 bases, short enough to survive sequencing errors.
pub const SCAFFOLD_KMER_LENGTH: usize = 12;
pub const SCAFFOLD_KMER_STEP: usize = 3;
/// Tag: the commonest `TAG_STARTS_TOP` read starts (first `TAG_PREFIX_BASES` bases) cover at least
/// `TAG_COVER_MIN` of reads, with under `TAG_HOUSEKEEPING_MAX` housekeeping reads.
pub const TAG_PREFIX_BASES: usize = 15;
pub const TAG_STARTS_TOP: usize = 200;
pub const TAG_COVER_MIN: f64 = 0.80;
pub const TAG_HOUSEKEEPING_MAX: f64 = 0.02;

const _: () = assert!(LIBRARY_READS <= SAMPLE_READS_MAX);
const _: () =
    assert!(SCAFFOLD_KMER_LENGTH <= 32 && TAG_PREFIX_BASES <= 32, "k-mers pack into a u64");

// ---------------------------------------------------------------------------------------------
// When no layout matches: what structure the reads show, worded as a hint.
// ---------------------------------------------------------------------------------------------

/// Anchor k-mers: 12-mers, long enough to be specific, short enough to fit every linker.
pub const ANCHOR_KMER_LENGTH: usize = 12;
/// Reads scanned for anchors; enough to see a linker in a few percent of reads.
pub const ANCHOR_READS: usize = 20_000;
/// Bases scanned per read: barcode-read structure lives in the first ~100 bases, and the cap keeps
/// long reads (kilobases) from dominating the cost.
pub const ANCHOR_SPAN: usize = 150;
/// An anchor in at least this share of reads means a fixed sequence at some position.
pub const ANCHOR_SHARE_MIN: f64 = 0.15;
/// Reads a file needs before its fingerprint carries a sketch or an anchor. Both are consensus
/// over the reads, and over few reads a consensus is the reads themselves: with 1 or 2 reads the
/// sketch spelled them out base for base (chaos probe, 8 Oct 2026). Over 1,000 reads a column
/// leads only where hundreds of molecules agree, which a kit's fixed sequence does and one
/// person's transcript does not. Far under the default sample of 100,000.
pub const FINGERPRINT_READS_MIN: usize = 1_000;
/// A run of this many identical bases makes an anchor low-complexity (polyA/polyT tails).
pub const ANCHOR_RUN_MAX: usize = 7;
/// At least this many profile columns dominated by one base mean fixed positions.
pub const FIXED_COLUMNS_MIN: usize = 8;
/// A column is dominated (fixed) when one base holds more than this share; the sketch prints it
/// upper case.
pub const FIXED_BASE_SHARE_MIN: f64 = 0.7;
/// The text hint calls a file's reads long when their most common length is over this. Routing
/// to the long-read detector uses `SHORT_READ_BASES_MAX` and `LONG_READ_SHARE_MIN`.
pub const LONG_READ_BASES_MIN: usize = 400;
/// A file whose reads are at least this long is cDNA-length, too long for a barcode read alone.
pub const CDNA_LENGTH_MIN: usize = 50;

const _: () = assert!(ANCHOR_KMER_LENGTH <= ANCHOR_SPAN);
const _: () = assert!(ANCHOR_RUN_MAX < ANCHOR_KMER_LENGTH);
const _: () = assert!(CDNA_LENGTH_MIN < LONG_READ_BASES_MIN);

// ---------------------------------------------------------------------------------------------
// The report for people.
// ---------------------------------------------------------------------------------------------

/// Candidates listed without `--all`.
pub const CANDIDATES_SHOWN: usize = 3;
/// Linker bases shown in an evidence line before `..`.
pub const LINKER_PREVIEW_BASES: usize = 8;
/// The profile sketch: '.' marks a random position (entropy above this many bits), lower case a
/// base that leads with more than `SKETCH_LEAD_SHARE_MIN` (upper case is `FIXED_BASE_SHARE_MIN`).
pub const SKETCH_RANDOM_BITS_MIN: f64 = 1.8;
pub const SKETCH_LEAD_SHARE_MIN: f64 = 0.4;

const _: () = assert!(SKETCH_LEAD_SHARE_MIN < FIXED_BASE_SHARE_MIN);

// ---------------------------------------------------------------------------------------------
// Long-read single-cell (10x cDNA read whole on Oxford Nanopore or PacBio). Set from the read
// structure and nanopore error rates, then checked on `bench/long_read` and its held-out test
// (docs/validation.md, "Long reads").
// ---------------------------------------------------------------------------------------------

/// The longest read Illumina produces (`MiSeq`, 2x301): no short read is longer.
pub const SHORT_READ_BASES_MAX: usize = 301;
/// A single file is long-read data when at least this share of its reads is longer than any short
/// read. Not the most common length: nanopore runs carry many short fragments (measured 7 Oct 2026:
/// two 10x 3' v3 nanopore runs had a most common length under 400 yet 62% and 68% of reads over
/// 301; Illumina has none).
pub const LONG_READ_SHARE_MIN: f64 = 0.2;
/// Reads the long-read detector examines.
pub const LONG_READS_EXAMINED: usize = 20_000;
/// Bases searched for the 10x Read 1 adapter at each end of a read: the adapter sits within the
/// first few dozen bases of a full-length molecule, more after a PacBio array adapter.
pub const LONG_ADAPTER_WINDOW: usize = 200;
/// Edits allowed in the 22-base Read 1 adapter (about 18%, above modern nanopore error rates).
pub const LONG_ADAPTER_EDITS_MAX: usize = 4;
/// Bases the barcode's start may sit either side of the adapter's end, for indel errors.
pub const LONG_BARCODE_SHIFT_MAX: usize = 2;
/// Edits allowed in the 13-base template-switch oligo that follows a 5' kit's UMI.
pub const LONG_TSO_EDITS_MAX: usize = 2;
/// A run of this many Ts marks the polyT that follows a 3' kit's UMI.
pub const LONG_POLY_T_RUN: usize = 6;
/// Bases after the barcode searched for the tail: the longest UMI (12), the TSO (13) and slack.
pub const LONG_TAIL_WINDOW: usize = 36;
/// A 10x claim needs the adapter in at least this share of reads.
pub const LONG_ADAPTER_SHARE_MIN: f64 = 0.1;
/// A barcode list is claimed when at least this share of adapter-bearing reads carries one of its
/// barcodes exactly, and at least `LONG_LIST_SEPARATION` times the next list's share; HIGH needs
/// `LONG_LIST_SHARE_HIGH`. Set against chance: a random 16-mer exactly matches the largest 10x list
/// (7.4M of 4^16 barcodes) at one of 5 shifts about 0.9% of the time, so the claim needs 5x chance
/// and HIGH 25x. A real barcode survives intact with probability (1 - error rate)^16: about 85% at
/// 1% errors (modern nanopore and PacBio HiFi), 31% at 7% (older nanopore chemistry).
pub const LONG_LIST_SHARE_MIN: f64 = 0.05;
pub const LONG_LIST_SHARE_HIGH: f64 = 0.25;
/// The separation is set by how much the lists share, not by error rate (errors lower every
/// list's share alike): 6.9% of the 5' v3 (GEM-X) list's barcodes are also on the 3' v4 list, and
/// real GEM-X runs show the other list at 7 to 9% of the true list's share (measured 7 Oct 2026 on
/// Oxford Nanopore 293T 3' v4 and Jurkat 5' v3 open data).
pub const LONG_LIST_SEPARATION: f64 = 5.0;
/// The tail (polyT or TSO) is decided when one is seen in at least this share of adapter-bearing
/// reads and `LONG_TAIL_SEPARATION` times the other.
pub const LONG_TAIL_SHARE_MIN: f64 = 0.2;
pub const LONG_TAIL_SEPARATION: f64 = 3.0;
/// The measured UMI length agrees with a kit's when within this many bases (indel errors).
pub const LONG_UMI_TOLERANCE: usize = 1;

const _: () = assert!(LONG_READS_EXAMINED <= SAMPLE_READS_MAX);
const _: () = assert!(LONG_LIST_SHARE_MIN < LONG_LIST_SHARE_HIGH);
const _: () = assert!(LONG_ADAPTER_EDITS_MAX < 22, "the adapter is 22 bases");
const _: () = assert!(LONG_TSO_EDITS_MAX < 13, "the TSO is 13 bases");

// ---------------------------------------------------------------------------------------------
// Sample-index reads (I1/I2), set aside when a sample's whole file set is given.
// ---------------------------------------------------------------------------------------------

/// Index reads are 6 to 10 bases (dual 10x indices are 10); a barcode read is at least 20 (Drop-
/// seq). A file whose reads are this short or shorter may be index reads.
pub const INDEX_READ_BASES_MAX: usize = 12;
/// Reads examined for the diversity test.
pub const INDEX_READS_EXAMINED: usize = 10_000;
/// One sample's index is one oligo, or four for a 10x single-index set: index reads are almost
/// all one of a handful of sequences, where barcode reads are nearly all distinct.
pub const INDEX_TOP_SEQUENCES: usize = 8;
pub const INDEX_COVER_MIN: f64 = 0.5;

const _: () = assert!(INDEX_READ_BASES_MAX < 20, "shorter than the shortest barcode read");
const _: () = assert!(INDEX_READS_EXAMINED <= SAMPLE_READS_MAX);
