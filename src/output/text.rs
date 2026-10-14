// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Prints a detection report: the read profiles and strand evidence, the top candidates with their
//! evidence, and the verdict.

use crate::ReadKind;
use crate::constants::{
    CANDIDATES_SHOWN, CDNA_LENGTH_MIN, FIXED_COLUMNS_MIN, LINKER_PREVIEW_BASES,
    LONG_READ_BASES_MIN, PROFILE_POSITIONS,
};
use crate::detect::Report;
use crate::detect::decide::Verdict;
use crate::detect::profile::{Profile, ratio};
use crate::detect::score::{Score, SegmentEvidence};
use crate::layouts::validation::{self, Level, Validation};
use crate::long_read::LongReport;
use crate::long_read::detect::{LongEvidence, LongVerdict, TENX_LISTS, umi_mode};
use crate::signals::instrument::{Evidence, Instrument};
use crate::signals::strand::Strand;
use crate::signals::structure::StructureSignal;
use std::fmt::Write as _;

/// Room for the whole report up front; a report is about 1.5 KB.
const OUTPUT_CAPACITY_BYTES: usize = 4096;

/// The report for people: read profiles, strand and library evidence, the top candidates (every
/// one with `all`) and the verdict, one line per `\n`.
pub fn render(report: &Report, all: bool) -> String {
    assert!(!report.scores.is_empty(), "every layout is scored, so a report has candidates");
    if report.swapped {
        assert!(report.r2_profile.is_some(), "a swap needs two files");
    }
    let mut out = String::with_capacity(OUTPUT_CAPACITY_BYTES);
    let out_ref = &mut out;
    write_header(out_ref, report);
    write_candidates(out_ref, report, all);
    write_verdict(out_ref, report);
    assert!(out.contains("VERDICT  "), "every report ends in a verdict");
    assert!(out.ends_with('\n'), "every report line is terminated");
    out
}

/// Appends one formatted line; a `String` sink cannot fail.
fn line(out: &mut String, text: std::fmt::Arguments) {
    let before = out.len();
    let Ok(()) = out.write_fmt(text) else { unreachable!("writing to a String cannot fail") };
    out.push('\n');
    assert!(out.len() > before, "a line adds at least its newline");
    assert!(out.ends_with('\n'), "a line ends the output with a newline");
}

fn write_profile(out: &mut String, label: &str, p: &Profile) {
    // Paired with Profile::of: a profile covers reads, and its lengths are ordered.
    assert!(p.reads > 0, "{label} profile covers no reads");
    assert!(p.len_min <= p.len_mode, "{label} shortest read longer than the mode");
    assert!(p.len_mode <= p.len_max, "{label} mode longer than the longest read");
    line(
        out,
        format_args!(
            "{label}  {} reads, length {} (min {}, max {})",
            p.reads, p.len_mode, p.len_min, p.len_max
        ),
    );
    line(out, format_args!("    {}", p.sketch()));
}

fn write_header(out: &mut String, report: &Report) {
    // Strand and library type both come from the cDNA read in one pass, so library typing never
    // has evidence the strand check lacked.
    if report.library.is_some() {
        assert!(report.cdna_strand.is_some(), "library typed without strand evidence");
    }
    if report.cdna_strand.is_some() {
        assert!(report.r2_profile.is_some(), "strand evidence needs a cDNA read, so two files");
    }
    write_profile(out, "R1", &report.r1_profile);
    if let Some(p2) = &report.r2_profile {
        write_profile(out, "R2", p2);
    }
    line(
        out,
        format_args!(
            "    ('.' random, upper case = one base over 70%, lower case = leads, '-' mixed)"
        ),
    );
    if let Some(strand) = &report.cdna_strand {
        let call = match strand.call() {
            Some((Strand::Sense, _)) => "sense (3' kits)",
            Some((Strand::Antisense, _)) => "antisense (5' kits)",
            None => "unclear",
        };
        line(
            out,
            format_args!(
                "    cDNA read strand: polyA runs {:.2}% vs polyT {:.2}%, chrM sense {:.2}% vs \
                antisense {:.2}%, nuclear genes sense {:.2}% vs antisense {:.2}% -> {call}",
                strand.poly_a * 100.0,
                strand.poly_t * 100.0,
                strand.mito_forward * 100.0,
                strand.mito_reverse * 100.0,
                strand.nuclear_sense * 100.0,
                strand.nuclear_antisense * 100.0
            ),
        );
    }
    if let Some(lib) = &report.library {
        line(
            out,
            format_args!(
                "    library: constant-region reads {:.1}%, housekeeping {:.1}%, sgRNA scaffold \
                {:.1}%, top tag starts {:.0}% -> {}",
                lib.constant_region * 100.0,
                lib.housekeeping * 100.0,
                lib.scaffold * 100.0,
                lib.tag_cover * 100.0,
                lib.call().label()
            ),
        );
    }
    if report.swapped {
        line(
            out,
            format_args!(
                "!!  R1 and R2 were SWAPPED: the files fit far better with their labels exchanged."
            ),
        );
    }
}

/// One line of evidence for a candidate: every segment, the UMI, the tail and start offsets.
fn evidence_line(s: &Score) -> String {
    assert!(
        s.evidence.len() <= s.layout.segments.len(),
        "{}: more evidence than segments",
        s.layout.name
    );
    if !s.offset_shares.is_empty() {
        assert!(s.layout.offsets.len() > 1, "{}: offset shares for one offset", s.layout.name);
    }
    let mut parts: Vec<String> = s
        .evidence
        .iter()
        .map(|e| match e {
            SegmentEvidence::Listed { position, list, rate, .. } => {
                format!("{}@{position} {:.0}%", list.label(), rate * 100.0)
            }
            SegmentEvidence::Unlisted { position, coverage, set, .. } => {
                format!("unlisted barcode@{position} top-{set} cover {:.0}%", coverage * 100.0)
            }
            SegmentEvidence::Fixed { position, sequence, alternatives, rate, .. } => {
                let shown: String =
                    String::from_utf8_lossy(sequence).chars().take(LINKER_PREVIEW_BASES).collect();
                let ell = if sequence.len() > LINKER_PREVIEW_BASES { ".." } else { "" };
                let more = if *alternatives > 1 {
                    format!(" (+{} alt)", alternatives - 1)
                } else {
                    String::new()
                };
                format!("{shown}{ell}{more}@{position} {:.0}%", rate * 100.0)
            }
        })
        .collect();
    if let Some((position, h)) = s.umi_entropy_min {
        parts.push(format!("UMI min {h:.2} bits@{position}"));
    }
    if let (Some(tail), Some(rate)) = (s.layout.tail(), s.tail_rate) {
        parts.push(format!("{} {:.0}% over {} base(s)", tail.label(), rate * 100.0, s.tail_seen));
    }
    if !s.offset_shares.is_empty() {
        let shares: Vec<String> =
            s.offset_shares.iter().map(|(o, f)| format!("+{o}:{:.0}%", f * 100.0)).collect();
        parts.push(format!("start offsets {}", shares.join(" ")));
    }
    let line = parts.join(" | ");
    assert!(!line.contains('\n'), "{}: evidence is one line", s.layout.name);
    line
}

fn write_candidates(out: &mut String, report: &Report, all: bool) {
    let shown = if all { report.scores.len() } else { CANDIDATES_SHOWN.min(report.scores.len()) };
    assert!(shown <= report.scores.len(), "shows {shown} of {} candidates", report.scores.len());
    if !all {
        assert!(shown <= CANDIDATES_SHOWN, "without --all at most three candidates, not {shown}");
    }
    line(
        out,
        format_args!(
            "\nCandidates (support = mean over barcode, linker and unlisted-barcode evidence):"
        ),
    );
    for s in report.scores.iter().take(shown) {
        line(
            out,
            format_args!(
                "  {:<32} on {:?}  support {:>5.1}%",
                s.layout.name,
                s.layout.read,
                s.support * 100.0
            ),
        );
        line(out, format_args!("      {}", evidence_line(s)));
        for c in &s.contradictions {
            line(out, format_args!("      x {c}"));
        }
        for n in &s.notes {
            line(out, format_args!("      ~ {n}"));
        }
    }
    if shown < report.scores.len() {
        line(out, format_args!("  ({} more; --all shows them)", report.scores.len() - shown));
    }
    if !report.unchecked.is_empty() {
        write_unchecked(out, report);
    }
}

/// One line naming the layouts not checked, written straight into the report.
fn write_unchecked(out: &mut String, report: &Report) {
    assert!(!report.unchecked.is_empty(), "only layouts not checked are named");
    let before = out.len();
    out.push_str("  not checked, list not fetched: ");
    for (i, layout) in report.unchecked.iter().enumerate() {
        if i > 0 {
            out.push_str(", ");
        }
        out.push_str(layout.name);
    }
    line(out, format_args!(" (scripts/fetch_whitelists.sh --cellranger)"));
    assert!(out.len() > before, "the line is written");
}

fn write_verdict(out: &mut String, report: &Report) {
    line(out, format_args!(""));
    match &report.verdict {
        Verdict::Match { best, confidence } => {
            let Some(s) = report.scores.get(*best) else {
                unreachable!("decide checks every verdict index names a score")
            };
            let confidence = confidence.label();
            line(out, format_args!("VERDICT  {}  (confidence {confidence})", s.layout.name));
            write_library_warning(out, report);
            line(out, format_args!("         geometry {}", s.layout.geometry));
            if !s.layout.note.is_empty() {
                line(out, format_args!("         note: {}", s.layout.note));
            }
            let record = validation::of(s.layout.name, ReadKind::Short);
            // The record vouches for HIGH calls; below HIGH it describes the layout, not this call.
            let scope = if confidence == "HIGH" { "" } else { " (for HIGH calls of this layout)" };
            line(out, format_args!("         validation{scope}: {}", validation_sentence(record)));
            write_family_validation(out, s.layout.name, ReadKind::Short, scope);
            for c in &s.contradictions {
                line(out, format_args!("         but: {c}"));
            }
        }
        Verdict::Ambiguous { candidates, reason } => {
            let names: Vec<_> = candidates
                .iter()
                .filter_map(|&i| report.scores.get(i))
                .map(|s| s.layout.name)
                .collect();
            // Paired with decide's postconditions.
            assert_eq!(names.len(), candidates.len(), "every candidate index names a score");
            assert!(names.len() >= 2, "an ambiguity names two layouts at least");
            line(out, format_args!("VERDICT  AMBIGUOUS between {}", names.join(" and ")));
            write_library_warning(out, report);
            line(out, format_args!("         {reason}"));
        }
        Verdict::NoMatch => {
            // Two ways to no match (decide): every candidate contradicted with no list agreeing,
            // or the best clean one under SUPPORT_MIN. They call for different next steps.
            if every_candidate_contradicted(report) {
                line(
                    out,
                    format_args!(
                        "VERDICT  no known layout (the reads contradict every candidate, and no \
                         vendor barcode list agrees with one)"
                    ),
                );
            } else {
                line(out, format_args!("VERDICT  no known layout (best support under 30%)"));
            }
            // The most common real-world cause: archives (ENA's two-file export especially) drop
            // the short barcode read and keep only cDNA.
            let hint = no_match_hint(report);
            for hint_line in &hint {
                line(out, format_args!("{hint_line}"));
            }
            // The link goes with the unknown-kit hint only: a dropped barcode read is not a kit.
            if hint.iter().any(|l| l == ASK_FOR_KIT) {
                line(out, format_args!("         {}", kit_request_url(report)));
            }
        }
    }
}

/// The unknown-kit hint's line that introduces the link, which follows it and only it.
const ASK_FOR_KIT: &str = concat!(
    "         open the link below: the form opens with the fingerprint filled in, ",
    "which holds no reads"
);

/// The unknown-kit form, which a link can fill in (GitHub issue forms take a field's id as a query
/// parameter: see .github/ISSUE_TEMPLATE/unknown-kit.yml).
const KIT_REQUEST_FORM: &str =
    "https://github.com/Tussl-Co/TusslChemistry/issues/new?template=unknown-kit.yml";
/// Browsers and GitHub take URLs well past this; a longer one drops the fingerprint and says so.
const URL_BYTES_MAX: usize = 8000;

/// Bytes an escaped byte takes: `%` and two hex digits.
const PERCENT_ENCODED_BYTES: usize = 3;

/// `text` percent-encoded for a URL query: unreserved characters stay, every other byte is `%XX`.
pub fn percent_encode(text: &str) -> String {
    let unreserved = |byte: u8| byte.is_ascii_alphanumeric() || b"-._~".contains(&byte);
    let mut out = String::with_capacity(text.len() * PERCENT_ENCODED_BYTES);
    for byte in text.bytes() {
        if unreserved(byte) {
            out.push(char::from(byte));
        } else {
            let _ = write!(out, "%{byte:02X}");
        }
    }
    let escaped = text.bytes().filter(|&b| !unreserved(b)).count();
    assert_eq!(out.len(), text.len() + escaped * (PERCENT_ENCODED_BYTES - 1), "each escape is %XX");
    assert_eq!(out.matches('%').count(), escaped, "every reserved byte is escaped, once");
    out
}

/// A link that opens the unknown-kit form with this run's fingerprint and the tool's version in.
fn kit_request_url(report: &Report) -> String {
    let fingerprint = crate::output::json::compact(&crate::output::json::fingerprint_value(report));
    let base = format!(
        "{KIT_REQUEST_FORM}&title={}&version={}",
        percent_encode("Kit request: "),
        percent_encode(env!("CARGO_PKG_VERSION"))
    );
    let full = format!("{base}&fingerprint={}", percent_encode(&fingerprint));
    let fits = full.len() <= URL_BYTES_MAX;
    let link = if fits { &full } else { &base };
    // The link is one URL token a terminal can open: every part went through percent_encode.
    assert!(link.len() <= URL_BYTES_MAX, "the link fits a browser, {} bytes", link.len());
    assert!(link.is_ascii() && !link.contains([' ', '\n']), "a link is one URL token");
    if fits {
        format!("ask for this kit: {link}")
    } else {
        format!("ask for this kit: {link}  (paste the fingerprint from --format json)")
    }
}

/// The report's first line: the instrument the reads' source records, and what an archive states
/// for the run, side by side when both are known (they disagree often enough to show both).
pub fn instrument_line(i: &Instrument, archive_model: Option<&str>) -> String {
    let named = match (i.platform, i.model.as_deref()) {
        (Some(p), Some(model)) => format!("{} {model}", p.label()),
        (Some(p), None) => p.label().to_owned(),
        (None, _) => "unknown".to_owned(),
    };
    let basis = match (i.evidence, i.id.as_deref()) {
        (Some(Evidence::ReadGroup), _) => " (from the BAM read group)".to_owned(),
        (Some(_), Some(id)) => format!(" ({id}, from the read names)"),
        (Some(_), None) => " (from the read names)".to_owned(),
        (None, _) => " (the read names carry none)".to_owned(),
    };
    let archive = archive_model.map_or(String::new(), |m| format!("; the archive records {m}"));
    let out = format!("Instrument  {named}{basis}{archive}\n");
    assert!(out.ends_with('\n') && out.matches('\n').count() == 1, "one line");
    assert!(i.platform.is_some() || !basis.contains("from the"), "an unknown names no source");
    out
}

/// The family record under the layout's own, for a layout in a kit family
/// (`validation::FAMILIES`): how often a HIGH call naming any of its layouts is that family.
fn write_family_validation(out: &mut String, layout: &str, kind: ReadKind, scope: &str) {
    let Some((family, record)) = validation::of_family(layout, kind) else { return };
    assert!(!family.is_empty(), "a family has a name");
    assert_eq!(record.family_only_calls, 0, "a family record has no family-only calls");
    line(out, format_args!("         {family} family{scope}: {}", validation_sentence(record)));
}

/// How far the layout's HIGH calls have been checked on held-out runs, in one sentence: the
/// level, the counts behind it, and what was set aside (family-only and disputed calls).
pub fn validation_sentence(v: Validation) -> String {
    let level = v.level();
    let mut sentence = match (level, v.error_bound()) {
        (Level::Unvalidated, _) => {
            "unvalidated: no held-out HIGH call of this exact layout yet".to_owned()
        }
        (Level::Contradicted, _) => format!(
            "contradicted: {} of {} held-out HIGH calls wrong",
            v.wrong,
            v.held_out_high_calls + v.wrong
        ),
        (_, Some(bound)) => format!(
            "{}: {} held-out HIGH calls from {} {}, none wrong (error under {:.1}% at 95%)",
            level.label(),
            v.held_out_high_calls,
            v.studies,
            if v.studies == 1 { "study" } else { "studies" },
            bound * 100.0
        ),
        (_, None) => unreachable!("a checked record without a wrong call has a bound: {v:?}"),
    };
    if v.family_only_calls > 0 {
        let more = if v.held_out_high_calls == 0 { "held-out HIGH calls" } else { "more" };
        let _ = write!(sentence, "; {} {more} right at kit-family level", v.family_only_calls);
    }
    if v.read_label_calls > 0 {
        let _ = write!(
            sentence,
            "; {} more agree only with labels read off the reads (not counted)",
            v.read_label_calls
        );
    }
    if v.disputed > 0 {
        let _ = write!(sentence, "; {} disputed (label names another kit)", v.disputed);
    }
    assert!(sentence.starts_with(level.label()), "the sentence opens with the level");
    assert!(!sentence.contains('\n'), "one line");
    sentence
}

/// A bead layout says nothing about what was captured: flag libraries a gene-expression mapper
/// must not take. Printed inside the verdict, so a pipeline reading it cannot miss it.
fn write_library_warning(out: &mut String, report: &Report) {
    let Some(lib) = &report.library else { return };
    // Library typing runs only once a layout matched or tied (detect's contract).
    assert!(!matches!(report.verdict, Verdict::NoMatch), "a library typed with no layout");
    assert!(report.cdna_strand.is_some(), "library typed without the cDNA pass behind it");
    let kind = lib.call();
    if kind != crate::signals::library::LibraryType::GeneExpression {
        line(
            out,
            format_args!(
                "         library: {}, NOT gene expression: use a pipeline for this library type",
                kind.label()
            ),
        );
    }
}

/// Why nothing matched, judged first from the shape of the input and then from structure in the
/// reads. No read statistic alone separates a cDNA read that opens with library sequence from the
/// barcode read of an unknown kit (measured 7 Oct 2026: some real cDNA files score like barcode
/// reads), so the number of files decides between "barcode read missing" and "unsupported kit".
pub fn no_match_hint(report: &Report) -> Vec<String> {
    assert!(matches!(report.verdict, Verdict::NoMatch), "a hint explains a no-match, nothing else");
    let profiles: Vec<&Profile> =
        std::iter::once(&report.r1_profile).chain(report.r2_profile.as_ref()).collect();
    assert!(
        report.structure.len() <= profiles.len(),
        "{} structure signals for {} files",
        report.structure.len(),
        profiles.len()
    );
    let strongest = report
        .structure
        .iter()
        .zip(["R1", "R2"])
        .filter(|(s, _)| s.structured())
        .max_by(|a, b| a.0.strength().total_cmp(&b.0.strength()));
    let out = if profiles.iter().any(|p| p.len_mode > LONG_READ_BASES_MIN) {
        no_match_hint_long_reads(&profiles)
    } else if report.r2_profile.is_none() && report.r1_profile.len_mode >= CDNA_LENGTH_MIN {
        no_match_hint_one_file(report, strongest)
    } else if let Some((s, file)) = strongest {
        no_match_hint_unknown_kit(s, file, every_candidate_contradicted(report))
    } else if profiles.iter().all(|p| p.len_mode >= CDNA_LENGTH_MIN) {
        no_match_hint_all_cdna(&profiles)
    } else {
        Vec::new()
    };
    // The JSON output trims this indent away, and the text report lines up under VERDICT.
    assert!(out.iter().all(|l| l.starts_with("         ")), "hint lines sit under the verdict");
    assert!(out.len() <= 5, "a hint is at most five lines, not {}", out.len());
    out
}

/// One line on a file's structure: its fixed columns and, when there is one, its anchor.
fn no_match_hint_describe(s: &StructureSignal) -> String {
    assert!(s.fixed_columns <= PROFILE_POSITIONS, "{} fixed columns", s.fixed_columns);
    if let Some((_, share)) = &s.anchor {
        assert!((0.0..=1.0).contains(share), "anchor share {share} outside 0..=1");
    }
    let anchor = s
        .anchor
        .as_ref()
        .map_or(String::new(), |(k, share)| format!(", {k} in {:.0}% of reads", share * 100.0));
    format!("{} fixed columns{anchor}", s.fixed_columns)
}

fn no_match_hint_long_reads(profiles: &[&Profile]) -> Vec<String> {
    assert!((1..=2).contains(&profiles.len()), "one or two files, not {}", profiles.len());
    assert!(profiles.iter().any(|p| p.len_mode > LONG_READ_BASES_MIN), "a file of long reads");
    let long = profiles.iter().filter(|p| p.len_mode > LONG_READ_BASES_MIN).count();
    // Two files of long reads, or one, go to the long-read detector: the files reaching this
    // hint mix long and short reads, or start with short reads that hide the long ones.
    let first = if profiles.len() == 2 && long == 1 {
        format!(
            "         hint: one file has long reads (over {LONG_READ_BASES_MIN} bases) and \
             the other short reads."
        )
    } else {
        format!(
            "         hint: reads are over {LONG_READ_BASES_MIN} bases, but the first reads \
             of the sample are short."
        )
    };
    vec![first, "         Give tchem the long-read files alone.".to_string()]
}

fn no_match_hint_one_file(
    report: &Report,
    strongest: Option<(&StructureSignal, &str)>,
) -> Vec<String> {
    assert!(report.r2_profile.is_none(), "the one-file hint is for one file");
    assert!(report.r1_profile.len_mode >= CDNA_LENGTH_MIN, "the one file is cDNA-length");
    let mut out = vec![
        "         hint: only one file, and it is cDNA-length. A single-cell library has \
            a barcode read too,"
            .to_string(),
        "         so it is probably missing: ENA's two-file export often drops it. Fetch \
            the SRA copy with"
            .to_string(),
        "         technical reads (fasterq-dump --include-technical).".to_string(),
    ];
    if let Some((s, file)) = strongest {
        assert_eq!(file, "R1", "one file has structure only in R1");
        out.push(format!(
            "         The file does show fixed-position structure ({}); if it is \
                itself a barcode",
            no_match_hint_describe(s)
        ));
        out.push("         read, the kit is not yet supported.".to_string());
    }
    out
}

/// Whether the no-match came from every candidate being contradicted, not from low support: the
/// scores sort clean layouts first, so a contradicted first score means none is clean.
fn every_candidate_contradicted(report: &Report) -> bool {
    assert!(matches!(report.verdict, Verdict::NoMatch), "asked of a no-match only");
    let contradicted = report.scores.first().is_some_and(|s| !s.contradictions.is_empty());
    assert!(
        !contradicted || report.scores.iter().all(|s| !s.contradictions.is_empty()),
        "a contradicted first score means every score is contradicted (detect's sort)"
    );
    contradicted
}

fn no_match_hint_unknown_kit(s: &StructureSignal, file: &str, contradicted: bool) -> Vec<String> {
    assert!(s.structured(), "the unknown-kit hint names a structured file");
    assert!(file == "R1" || file == "R2", "a hint names R1 or R2, not {file}");
    // When a known layout fits but for a contradiction (a polyT a few bases late, say), "not yet
    // supported" would send the user the wrong way: point at the contradictions instead.
    let (why, then) = if contradicted {
        (
            "         The reads contradict every candidate (see each candidate's contradictions). \
             It can be",
            "         a variant of a known kit, trimmed reads, or a kit not yet supported. To ask \
             for the kit,",
        )
    } else {
        (
            "         No known layout fits it: probably a single-cell kit not yet supported (or a \
             non-single-cell",
            "         library such as amplicons). scripts/discover.py can describe it. To ask \
             for the kit,",
        )
    };
    // A file is structured by its fixed columns or by one sequence in many reads, at any position
    // (a linker, an adapter or a primer): only the first is a fixed-position layout.
    let first = match (&s.anchor, s.fixed_columns >= FIXED_COLUMNS_MIN) {
        (_, true) => format!(
            "         hint: {file} has fixed-position structure ({}), as a barcode read does.",
            no_match_hint_describe(s)
        ),
        (Some((kmer, share)), false) => format!(
            "         hint: {file} carries one sequence in {:.0}% of reads ({kmer}): a linker, \
             an adapter or a primer.",
            share * 100.0
        ),
        (None, false) => unreachable!("a structured file without fixed columns has an anchor"),
    };
    assert!(
        !first.contains("fixed-position") || s.fixed_columns >= FIXED_COLUMNS_MIN,
        "a hint calls a file fixed-position only with {FIXED_COLUMNS_MIN} fixed columns: {first}"
    );
    vec![
        first,
        why.to_string(),
        then.to_string(),
        ASK_FOR_KIT.to_string(),
        "         (CONTRIBUTING.md; it is also \"fingerprint\" in --format json).".to_string(),
    ]
}

fn no_match_hint_all_cdna(profiles: &[&Profile]) -> Vec<String> {
    assert!((1..=2).contains(&profiles.len()), "one or two files, not {}", profiles.len());
    assert!(profiles.iter().all(|p| p.len_mode >= CDNA_LENGTH_MIN), "every file is cDNA-length");
    vec![
        "         hint: every file looks like cDNA (no short barcode read). If this is \
            single-cell"
            .to_string(),
        "         data, the barcode read was probably dropped by the archive; ENA's \
            two-file export"
            .to_string(),
        "         often omits it. Fetch the SRA copy with technical reads (fasterq-dump \
            --include-technical)."
            .to_string(),
    ]
}

/// The report for people on a long-read sample: read lengths, what the adapter search found, and
/// the verdict with the setting a long-read pipeline needs.
pub fn render_long(report: &LongReport) -> String {
    let e = &report.evidence;
    assert!(e.reads > 0, "a long report examined reads");
    assert!(e.adapter <= e.reads, "adapter reads are a subset of reads");
    let p = &report.profile;
    let mut out = String::with_capacity(OUTPUT_CAPACITY_BYTES);
    let out_ref = &mut out;
    line(
        out_ref,
        format_args!(
            "Long reads  {} sampled, length {} (min {}, max {}); {} examined",
            p.reads, p.len_mode, p.len_min, p.len_max, e.reads
        ),
    );
    write_long_evidence(out_ref, e);
    line(out_ref, format_args!(""));
    write_long_verdict(out_ref, &report.verdict);
    assert!(out.contains("VERDICT  "), "every report ends in a verdict");
    out
}

fn write_long_evidence(out: &mut String, e: &LongEvidence) {
    assert!(e.forward <= e.adapter, "forward reads are a subset of adapter reads");
    assert!(e.poly_t + e.tso <= e.adapter, "one tail per adapter read at most");
    let share = ratio(e.adapter, e.usable) * 100.0;
    if e.adapter == 0 {
        line(
            out,
            format_args!(
                "    10x Read 1 adapter in {share:.1}% of usable reads ({} of {} long enough to \
                 show it)",
                e.usable, e.reads
            ),
        );
    } else {
        line(
            out,
            format_args!(
                "    10x Read 1 adapter in {share:.1}% of usable reads ({:.1}% at the read start, \
                 the rest reverse-complemented at the end)",
                ratio(e.forward, e.adapter) * 100.0
            ),
        );
    }
    let lists: Vec<String> = TENX_LISTS
        .iter()
        .zip(e.list_hits)
        .zip(e.list_exclusive)
        .map(|((id, hits), alone)| {
            let (all, only) = (ratio(hits, e.adapter) * 100.0, ratio(alone, e.adapter) * 100.0);
            format!("{} {all:.1}% ({only:.1}% on it alone)", id.label())
        })
        .collect();
    line(out, format_args!("    barcodes listed exactly: {}", lists.join(" | ")));
    let umi = umi_mode(e).map_or_else(|| "not measured".to_string(), |u| format!("{u} bases"));
    line(
        out,
        format_args!(
            "    after the UMI: polyT {:.1}%, TSO {:.1}%; UMI length {umi}",
            ratio(e.poly_t, e.adapter) * 100.0,
            ratio(e.tso, e.adapter) * 100.0
        ),
    );
}

fn write_long_verdict(out: &mut String, verdict: &LongVerdict) {
    let before = out.len();
    match verdict {
        LongVerdict::Match { kit, confidence, contradictions } => {
            let confidence = confidence.label();
            line(
                out,
                format_args!("VERDICT  {}, long reads  (confidence {confidence})", kit.layout.name),
            );
            if let Some(kit) = kit.layout.long_read_kit {
                line(out, format_args!("         wf-single-cell --kit {kit}"));
            }
            let record = validation::of(kit.layout.name, ReadKind::Long);
            // The record vouches for HIGH calls; below HIGH it describes the layout, not this call.
            let scope = if confidence == "HIGH" { "" } else { " (for HIGH calls of this layout)" };
            line(out, format_args!("         validation{scope}: {}", validation_sentence(record)));
            write_family_validation(out, kit.layout.name, ReadKind::Long, scope);
            for c in contradictions {
                line(out, format_args!("         but: {c}"));
            }
        }
        LongVerdict::Ambiguous { kits, reason } => {
            let names: Vec<&str> = kits.iter().map(|k| k.layout.name).collect();
            line(
                out,
                format_args!("VERDICT  AMBIGUOUS between {}, long reads", names.join(" and ")),
            );
            line(out, format_args!("         {reason}"));
        }
        LongVerdict::NoMatch { reason } => {
            line(out, format_args!("VERDICT  no known layout, long reads"));
            line(out, format_args!("         hint: {reason}"));
        }
    }
    assert!(out.len() > before, "a verdict writes at least one line");
    assert!(out.ends_with('\n'), "the verdict's lines are terminated");
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Unreserved characters stay; everything else, multi-byte UTF-8 included, is %XX per byte.
    #[test]
    fn percent_encoding() {
        assert_eq!(percent_encode("Kit request: "), "Kit%20request%3A%20");
        assert_eq!(percent_encode("a-b_c.d~e"), "a-b_c.d~e");
        assert_eq!(percent_encode("{\"x\":\"é&=\"}"), "%7B%22x%22%3A%22%C3%A9%26%3D%22%7D");
    }

    /// A file structured only by one common sequence is not called fixed-position: the sequence
    /// can be a primer at any position (user test, 12 Oct 2026).
    #[test]
    fn a_common_sequence_alone_is_not_a_fixed_layout() {
        let anchor = Some(("TTACGTATTGCT".to_owned(), 0.51));
        let common = StructureSignal { fixed_columns: 2, anchor: anchor.clone() };
        let hint = no_match_hint_unknown_kit(&common, "R2", true);
        assert!(hint[0].contains("one sequence in 51% of reads (TTACGTATTGCT)"), "{hint:?}");
        assert!(!hint[0].contains("fixed-position"), "{hint:?}");
        let fixed = StructureSignal { fixed_columns: FIXED_COLUMNS_MIN, anchor };
        let hint = no_match_hint_unknown_kit(&fixed, "R1", true);
        assert!(hint[0].contains("R1 has fixed-position structure (8 fixed columns"), "{hint:?}");
    }
}
