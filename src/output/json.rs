// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Machine-readable output: one JSON object a pipeline can act on without a person reading the
//! report. It names the layout and its geometry, the barcode list files, which input file is the
//! barcode read and which the cDNA (after swap detection), the strand to give the mapper, and the
//! library type, with every signal behind the call.
//!
//! Written by hand: a small value tree and a pretty printer with two-space indents, `"key": value`
//! pairs, `[]`/`{}` when empty, strings escaped per RFC 8259, counts as integers, and shares in
//! their shortest round-trip form. The layout matches `serde_json`'s pretty printer byte for byte:
//! over every k/m share a run can produce (m up to a million reads) and 5 million random doubles,
//! the only differences are ties between two equally short digit strings above 1e13, which a share
//! never reaches.

use crate::ReadKind;
use crate::constants::{ANCHOR_SHARE_MIN, FINGERPRINT_READS_MIN};
use crate::detect::Report;
use crate::detect::decide::{Confidence, Verdict};
use crate::detect::profile::{Profile, ratio};
use crate::detect::score::Score;
use crate::layouts::validation::{self, Level, Validation};
use crate::layouts::{Layout, Mate, Tail};
use crate::long_read::LongReport;
use crate::long_read::detect::{LongVerdict, TENX_LISTS, umi_mode};
use crate::signals::instrument::{Evidence, Instrument, Platform};
use crate::signals::library::{LibraryEvidence, LibraryType};
use crate::signals::strand::{Strand, StrandEvidence};
use crate::signals::structure::StructureSignal;
use std::fmt::{Display, Write as _};
use std::path::Path;

/// The output contract's version, separate from the tool's: pipelines key on this, not on the
/// tool version. A new key is a minor bump; renaming, removing or retyping a key is a major
/// bump. `schema/report.v1.json` describes it, and the gate checks every recorded report
/// against that file.
pub const SCHEMA_VERSION: &str = "1.1";
/// The `--emit --format json` object's own contract version (schema/emit.v1.json), separate from
/// the report's: the two outputs change independently.
pub const EMIT_SCHEMA_VERSION: &str = "1.0";
/// The fingerprint's own format version: bumped when what it holds changes.
const FINGERPRINT_VERSION: u64 = 1;
/// Deepest nesting the output uses (report, fingerprint, files, file, anchor); the printer
/// refuses deeper trees.
const DEPTH_MAX: usize = 5;
/// Room for the whole object up front; a report's JSON is about 1.5 KB.
const OUTPUT_CAPACITY_BYTES: usize = 4096;
/// An f64's shortest round-trip form has at most 17 significant digits.
const F64_DIGITS_MAX: i32 = 17;
/// `serde_json`'s layout: plain decimals while the point falls within 16 digits of the start and
/// at most 4 zeros after `0.`; exponent form outside.
const PLAIN_POINT_MAX: i32 = 16;
const PLAIN_LEADING_ZEROS_MAX: i32 = 4;

#[derive(Debug, PartialEq)]
pub enum Value {
    Null,
    Bool(bool),
    Integer(u64),
    /// A fraction of reads, 0 to 1.
    Share(f64),
    String(String),
    Array(Vec<Value>),
    Object(Vec<(&'static str, Value)>),
}

/// Why the output could not be written.
#[derive(Debug)]
pub enum JsonRefusal {
    /// JSON strings are Unicode; a path that is not cannot be named faithfully.
    NonUnicodePath(std::path::PathBuf),
}

impl std::fmt::Display for JsonRefusal {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonUnicodePath(path) => {
                write!(f, "path contains invalid UTF-8 characters: {}", path.display())
            }
        }
    }
}

impl std::error::Error for JsonRefusal {}

fn text(value: &str) -> Value {
    Value::String(value.to_owned())
}

/// A file named in the output. Every path here names a file the run opened or joined onto the
/// whitelist directory, so it is never empty.
fn path(value: &Path) -> Result<Value, JsonRefusal> {
    assert!(!value.as_os_str().is_empty(), "a named file has a non-empty path");
    let Some(unicode) = value.to_str() else {
        return Err(JsonRefusal::NonUnicodePath(value.to_path_buf()));
    };
    assert_eq!(unicode.len(), value.as_os_str().len(), "a Unicode path keeps its bytes");
    Ok(text(unicode))
}

/// The keys that open every object the tool prints: which tool, which version of it, and which
/// version of that object's contract (`version_key`). One place, so the report's two shapes and
/// the emitted command cannot drift apart.
fn opened(
    version_key: &'static str,
    contract: &str,
    fields: impl IntoIterator<Item = (&'static str, Value)>,
) -> Value {
    let mut object = vec![
        ("tool", text("tussl-chemistry")),
        ("version", text(env!("CARGO_PKG_VERSION"))),
        (version_key, text(contract)),
    ];
    let opening = object.len();
    object.extend(fields);
    assert!(object.len() > opening, "an object carries more than its opening keys");
    assert!(
        object
            .iter()
            .skip(opening)
            .all(|(key, _)| !matches!(*key, "tool" | "version") && *key != version_key),
        "the opening keys appear once, first"
    );
    Value::Object(object)
}

/// A report's opening keys, then `fields`.
fn report_object(fields: impl IntoIterator<Item = (&'static str, Value)>) -> Value {
    let report = opened("schema_version", SCHEMA_VERSION, fields);
    let Value::Object(fields) = &report else { unreachable!("opened builds an object") };
    // A key written twice would be valid JSON that most parsers read as the last one only.
    let unique =
        fields.iter().enumerate().all(|(i, (k, _))| fields.iter().take(i).all(|(j, _)| j != k));
    assert!(unique, "a report's keys are distinct");
    assert!(fields.iter().any(|(k, _)| *k == "verdict"), "a report states its verdict");
    report
}

/// The `--emit --format json` object: the command as an argument list (a placeholder appears as
/// `${NAME}`, to be replaced by the caller), what each placeholder must hold, notes, and the shell
/// form. Not a report: it has its own schema (schema/emit.v1.json).
pub fn emit_value(e: &crate::Emitted) -> Value {
    let argv = e.argv.iter().map(|arg| match arg {
        crate::Arg::Text(t) => text(t),
        crate::Arg::Var(p) => Value::String(format!("${{{}}}", p.name)),
    });
    let placeholders = e
        .placeholders
        .iter()
        .map(|p| Value::Object(vec![("name", text(p.name)), ("meaning", text(p.meaning))]));
    let value = opened(
        "emit_schema_version",
        EMIT_SCHEMA_VERSION,
        [
            ("pipeline", text(e.pipeline.name())),
            ("layout", text(e.layout)),
            ("confidence", text(&e.confidence.to_lowercase())),
            ("validation", text(&e.validation)),
            ("argv", Value::Array(argv.collect())),
            ("placeholders", Value::Array(placeholders.collect())),
            ("notes", Value::Array(e.notes.iter().map(|n| text(n)).collect())),
            ("shell", Value::String(crate::render_shell(e))),
        ],
    );
    assert!(!e.argv.is_empty(), "a command has a program");
    assert!(
        e.placeholders
            .iter()
            .all(|p| e.argv.iter().any(|a| matches!(a, crate::Arg::Var(q) if q.name == p.name))),
        "every declared placeholder is used in argv"
    );
    value
}

/// The JSON object for `report`; `r1`/`r2` are the files as labelled on the command line.
pub fn build(
    report: &Report,
    r1: &Path,
    r2: Option<&Path>,
    whitelists: &Path,
) -> Result<Value, JsonRefusal> {
    let (verdict, confidence, chosen, candidates) = verdict_parts(report);
    assert!(chosen.is_none() || verdict == "match", "only a match names a layout");
    assert!(candidates.is_empty() || verdict == "ambiguous", "only an ambiguity lists candidates");
    let layout = match chosen {
        Some(score) => layout_value(score.layout, ReadKind::Short, whitelists)?,
        None => Value::Null,
    };
    let files = match chosen {
        Some(score) => files_value(report, score, r1, r2)?,
        None => Value::Null,
    };
    let cdna_strand = if verdict == "match" { mapper_strand(report, chosen) } else { None };
    let hint = hint_sentence(report);
    if hint.is_some() {
        assert_eq!(verdict, "no_match", "only a no-match carries a hint");
    }
    let no_strings: &[String] = &[];
    let strings = |values: &[String]| Value::Array(values.iter().map(|s| text(s)).collect());
    Ok(report_object([
        ("verdict", text(verdict)),
        ("confidence", confidence.map_or(Value::Null, text)),
        ("layout", layout),
        ("candidates", Value::Array(candidates.into_iter().map(text).collect())),
        ("files", files),
        ("cdna_strand", cdna_strand.map_or(Value::Null, text)),
        ("library", report.library.as_ref().map_or(Value::Null, library_value)),
        ("swapped", Value::Bool(report.swapped)),
        ("strand_evidence", report.cdna_strand.as_ref().map_or(Value::Null, strand_value)),
        ("contradictions", strings(chosen.map_or(no_strings, |s| &s.contradictions))),
        ("notes", strings(chosen.map_or(no_strings, |s| &s.notes))),
        ("hint", hint.map_or(Value::Null, Value::String)),
        ("not_checked", Value::Array(report.unchecked.iter().map(|l| text(l.name)).collect())),
        ("reads_sampled", Value::Integer(report.r1_profile.reads as u64)),
        ("read_kind", text("short")),
        ("long_read_evidence", Value::Null),
        ("fingerprint", fingerprint_value(report)),
    ]))
}

/// A no-match's fingerprint: what each file shows, to share when asking for a kit
/// the tool does not know. It holds no read: per position only a consensus over every read
/// sampled, and only for files with a barcode read's structure (a cDNA file gets its lengths);
/// a k-mer only when it starts at least `ANCHOR_SHARE_MIN` of reads, which a linker does and a
/// person's transcript does not; and neither from fewer than `FINGERPRINT_READS_MIN` reads,
/// whose consensus is the reads. Null unless the verdict is a no-match.
pub(crate) fn fingerprint_value(report: &Report) -> Value {
    if report.structure.is_empty() {
        assert!(!matches!(report.verdict, Verdict::NoMatch), "a no-match has structure signals");
        return Value::Null;
    }
    let profiles = std::iter::once(&report.r1_profile).chain(report.r2_profile.as_ref());
    let files: Vec<Value> = ["R1", "R2"]
        .into_iter()
        .zip(profiles)
        .zip(&report.structure)
        .map(|((label, profile), structure)| fingerprint_file(label, profile, structure))
        .collect();
    assert_eq!(files.len(), report.structure.len(), "one entry per file");
    Value::Object(vec![
        ("fingerprint_version", Value::Integer(FINGERPRINT_VERSION)),
        ("files", Value::Array(files)),
    ])
}

fn fingerprint_file(label: &str, p: &Profile, s: &StructureSignal) -> Value {
    let structured = s.structured();
    // Consensus over few reads is the reads: below the floor, neither sketch nor anchor.
    let enough = p.reads >= FINGERPRINT_READS_MIN;
    let anchor = s.anchor.as_ref().filter(|(_, share)| enough && *share >= ANCHOR_SHARE_MIN);
    if let Some((kmer, share)) = anchor {
        // The privacy rule, checked where it is applied: only a k-mer shared by many reads.
        assert!(*share >= ANCHOR_SHARE_MIN, "{label}: anchor at {share} below the share floor");
        assert!(p.reads >= FINGERPRINT_READS_MIN, "{label}: anchor over {} reads", p.reads);
        assert!(kmer.bytes().all(|b| b"ACGT".contains(&b)), "{label}: anchor {kmer}");
    }
    let sketch = (structured && enough).then(|| p.sketch());
    Value::Object(vec![
        ("file", text(label)),
        ("reads", Value::Integer(p.reads as u64)),
        ("length_min", Value::Integer(p.len_min as u64)),
        ("length_mode", Value::Integer(p.len_mode as u64)),
        ("length_max", Value::Integer(p.len_max as u64)),
        ("structured", Value::Bool(structured)),
        ("fixed_columns", Value::Integer(s.fixed_columns as u64)),
        ("sketch", sketch.as_deref().map_or(Value::Null, text)),
        (
            "anchor",
            anchor.map_or(Value::Null, |(kmer, share)| {
                Value::Object(vec![("kmer", text(kmer)), ("share", Value::Share(*share))])
            }),
        ),
    ])
}

/// `report` with its `instrument` object placed after `read_kind`: what the reads' source records
/// (platform, model, the instrument or flowcell id, and what that rests on), and the model an
/// archive states for the run when the reads came from one.
pub fn with_instrument(report: Value, i: &Instrument, archive_model: Option<&str>) -> Value {
    let Value::Object(mut fields) = report else { unreachable!("a report is an object") };
    let platform = i.platform.map(|p| match p {
        Platform::Illumina => "illumina",
        Platform::Element => "element",
        Platform::Mgi => "mgi",
        Platform::Ultima => "ultima",
        Platform::PacBio => "pacbio",
        Platform::OxfordNanopore => "oxford_nanopore",
    });
    let from = i.evidence.map(|e| match e {
        Evidence::ReadName => "read_name",
        Evidence::ReadGroup => "read_group",
    });
    let optional = |v: Option<&str>| v.map_or(Value::Null, text);
    let value = Value::Object(vec![
        ("platform", optional(platform)),
        ("model", optional(i.model.as_deref())),
        ("id", optional(i.id.as_deref())),
        ("from", optional(from)),
        ("archive_model", optional(archive_model)),
    ]);
    let Some(at) = fields.iter().position(|(key, _)| *key == "read_kind") else {
        unreachable!("every report has read_kind")
    };
    fields.insert(at + 1, ("instrument", value));
    assert!(fields.iter().filter(|(k, _)| *k == "instrument").count() == 1, "one instrument key");
    assert!(platform.is_some() == from.is_some(), "a platform rests on evidence");
    Value::Object(fields)
}

/// The no-match hint as one sentence, as JSON and the library give it: only a no-match carries one,
/// as in the text report.
pub fn hint_sentence(report: &Report) -> Option<String> {
    if !matches!(report.verdict, Verdict::NoMatch) {
        return None;
    }
    let lines = crate::output::text::no_match_hint(report);
    let sentence = (!lines.is_empty()).then(|| {
        lines.iter().map(|l| l.trim().trim_start_matches("hint: ")).collect::<Vec<_>>().join(" ")
    });
    assert!(sentence.as_ref().is_none_or(|s| !s.is_empty()), "a hint says something");
    assert!(sentence.as_ref().is_none_or(|s| !s.contains('\n')), "a hint is one line");
    sentence
}

/// The JSON object for a long-read `report`, with the same keys as a short-read one plus the
/// long-read evidence. `file` is the one sample's file; it holds barcode and cDNA alike.
pub fn build_long(
    report: &LongReport,
    file: &Path,
    whitelists: &Path,
) -> Result<Value, JsonRefusal> {
    let (verdict, confidence, layout, candidates, contradictions, hint) = match &report.verdict {
        LongVerdict::Match { kit, confidence, contradictions } => {
            let layout = layout_value(kit.layout, ReadKind::Long, whitelists)?;
            (
                "match",
                Some(confidence_word(*confidence)),
                layout,
                Vec::new(),
                contradictions.clone(),
                None,
            )
        }
        LongVerdict::Ambiguous { kits, .. } => {
            let names = kits.iter().map(|k| text(k.layout.name)).collect();
            ("ambiguous", None, Value::Null, names, Vec::new(), None)
        }
        LongVerdict::NoMatch { reason } => {
            ("no_match", None, Value::Null, Vec::new(), Vec::new(), Some(reason.clone()))
        }
    };
    let files = if verdict == "match" {
        Value::Object(vec![("barcode_read", path(file)?), ("cdna_read", path(file)?)])
    } else {
        Value::Null
    };
    assert!(report.evidence.reads > 0, "a long report examined reads");
    assert!(candidates.is_empty() || verdict == "ambiguous", "only an ambiguity lists candidates");
    Ok(report_object([
        ("verdict", text(verdict)),
        ("confidence", confidence.map_or(Value::Null, text)),
        ("layout", layout),
        ("candidates", Value::Array(candidates)),
        ("files", files),
        ("cdna_strand", Value::Null),
        ("library", Value::Null),
        ("swapped", Value::Bool(false)),
        ("strand_evidence", Value::Null),
        ("contradictions", Value::Array(contradictions.iter().map(|c| text(c)).collect())),
        ("notes", Value::Array(Vec::new())),
        ("hint", hint.map_or(Value::Null, Value::String)),
        // Long reads are tested against lists every install has.
        ("not_checked", Value::Array(Vec::new())),
        ("reads_sampled", Value::Integer(report.profile.reads as u64)),
        ("read_kind", text("long")),
        ("long_read_evidence", long_evidence_value(report)),
        // Long reads that match nothing get a reason, not a fingerprint (not yet).
        ("fingerprint", Value::Null),
    ]))
}

fn long_evidence_value(report: &LongReport) -> Value {
    let e = &report.evidence;
    assert!(e.adapter <= e.reads, "adapter reads are a subset of reads");
    assert!(e.forward <= e.adapter, "forward reads are a subset of adapter reads");
    let lists = TENX_LISTS
        .iter()
        .zip(e.list_hits)
        .zip(e.list_exclusive)
        .map(|((id, hits), alone)| {
            Value::Object(vec![
                ("list", text(id.label())),
                ("share", Value::Share(ratio(hits, e.adapter))),
                // Barcodes on this list and no other: what the lists are told apart on.
                ("exclusive_share", Value::Share(ratio(alone, e.adapter))),
            ])
        })
        .collect();
    Value::Object(vec![
        ("reads_examined", Value::Integer(e.reads as u64)),
        ("usable_reads", Value::Integer(e.usable as u64)),
        ("adapter_share", Value::Share(ratio(e.adapter, e.usable))),
        ("adapter_at_start_share", Value::Share(ratio(e.forward, e.adapter))),
        ("barcode_list_shares", Value::Array(lists)),
        ("poly_t_share", Value::Share(ratio(e.poly_t, e.adapter))),
        ("tso_share", Value::Share(ratio(e.tso, e.adapter))),
        ("umi_length", umi_mode(e).map_or(Value::Null, |u| Value::Integer(u as u64))),
        ("read_length_mode", Value::Integer(report.profile.len_mode as u64)),
    ])
}

/// A confidence as the JSON names it.
fn confidence_word(confidence: Confidence) -> &'static str {
    match confidence {
        Confidence::High => "high",
        Confidence::Medium => "medium",
        Confidence::Low => "low",
    }
}

/// The verdict word, a match's confidence and layout, and an ambiguity's candidate names.
fn verdict_parts(
    report: &Report,
) -> (&'static str, Option<&'static str>, Option<&Score>, Vec<&'static str>) {
    match &report.verdict {
        Verdict::Match { best, confidence } => {
            let Some(score) = report.scores.get(*best) else {
                unreachable!("verdict index {best} past {} scores", report.scores.len())
            };
            ("match", Some(confidence_word(*confidence)), Some(score), Vec::new())
        }
        Verdict::Ambiguous { candidates, .. } => {
            let names: Vec<&'static str> = candidates
                .iter()
                .filter_map(|&i| report.scores.get(i))
                .map(|s| s.layout.name)
                .collect();
            assert_eq!(names.len(), candidates.len(), "every candidate index names a score");
            ("ambiguous", None, None, names)
        }
        Verdict::NoMatch => ("no_match", None, None, Vec::new()),
    }
}

fn layout_value(layout: &Layout, kind: ReadKind, whitelists: &Path) -> Result<Value, JsonRefusal> {
    // A geometry describes both reads, read 1 first, as simpleaf expects.
    assert!(layout.geometry.starts_with("1{"), "{}: geometry opens on read 1", layout.name);
    assert!(layout.geometry.contains("}2{"), "{}: geometry covers read 2", layout.name);
    let lists = layout.lists().iter().map(|id| path(&whitelists.join(id.file_name())));
    Ok(Value::Object(vec![
        ("name", text(layout.name)),
        // simpleaf / piscem / alevin-fry geometry string.
        ("geometry", text(layout.geometry)),
        ("barcode_lists", Value::Array(lists.collect::<Result<_, _>>()?)),
        ("note", text(layout.note)),
        // Oxford Nanopore wf-single-cell's --kit, for 10x kits.
        ("long_read_kit", layout.long_read_kit.map_or(Value::Null, text)),
        ("validation", validation_value(validation::of(layout.name, kind))),
        ("family_validation", family_validation_value(layout.name, kind)),
    ]))
}

/// The record of the layout's kit family (`validation::FAMILIES`): how often a HIGH call naming
/// any layout of the family is that family. Null for a layout in no family.
fn family_validation_value(layout: &str, kind: ReadKind) -> Value {
    let Some((family, record)) = validation::of_family(layout, kind) else { return Value::Null };
    let Value::Object(mut fields) = validation_value(record) else {
        unreachable!("validation_value builds an object")
    };
    // The family record counts no family-only calls: every right call is right for the family.
    let before = fields.len();
    fields.retain(|(key, _)| *key != "family_only_calls");
    assert_eq!(fields.len() + 1, before, "the family object drops exactly family_only_calls");
    fields.insert(0, ("family", text(family)));
    Value::Object(fields)
}

/// How far the layout's HIGH calls have been checked on held-out runs. It
/// describes the layout, not this run, and vouches for HIGH calls only (`applies_to`).
fn validation_value(v: Validation) -> Value {
    assert!(v.studies <= v.held_out_high_calls, "more studies than runs: {v:?}");
    let level = v.level();
    let bound = v.error_bound();
    assert!(
        bound.is_some()
            == matches!(level, Level::Checked | Level::Validated | Level::ErrorUnder3Pct),
        "a bound exactly when studies were checked and nothing was wrong"
    );
    Value::Object(vec![
        ("applies_to", text("high")),
        ("level", text(level.label())),
        ("held_out_high_calls", Value::Integer(u64::from(v.held_out_high_calls))),
        ("studies", Value::Integer(u64::from(v.studies))),
        ("family_only_calls", Value::Integer(u64::from(v.family_only_calls))),
        ("disputed", Value::Integer(u64::from(v.disputed))),
        ("wrong", Value::Integer(u64::from(v.wrong))),
        ("error_bound_95", bound.map_or(Value::Null, Value::Share)),
    ])
}

/// Which file holds the barcode read and which the cDNA, in the orientation detection chose.
fn files_value(
    report: &Report,
    score: &Score,
    r1: &Path,
    r2: Option<&Path>,
) -> Result<Value, JsonRefusal> {
    assert!(!report.swapped || r2.is_some(), "a swap needs two files");
    let read_of = |mate: Mate| -> Option<&Path> {
        let labelled_r1 = matches!(mate, Mate::R1) != report.swapped;
        if labelled_r1 { Some(r1) } else { r2 }
    };
    let other = if matches!(score.layout.read, Mate::R1) { Mate::R2 } else { Mate::R1 };
    let Some(barcode_read) = read_of(score.layout.read) else { return Ok(Value::Null) };
    let cdna_read = match read_of(other) {
        Some(cdna) => path(cdna)?,
        None => Value::Null,
    };
    if r2.is_none() {
        assert!(cdna_read == Value::Null, "one file is the barcode read, so no cDNA file is named");
    }
    Ok(Value::Object(vec![("barcode_read", path(barcode_read)?), ("cdna_read", cdna_read)]))
}

/// Strand for the mapper: the measured cDNA strand when it is clear, else the layout's tail.
pub fn mapper_strand(report: &Report, chosen: Option<&Score>) -> Option<&'static str> {
    assert!(chosen.is_some(), "only a match gives the mapper a strand");
    let measured = report.cdna_strand.as_ref().and_then(StrandEvidence::call);
    let strand = measured
        .map(|(strand, _)| match strand {
            Strand::Sense => "forward",
            Strand::Antisense => "reverse",
        })
        .or_else(|| match chosen.and_then(|s| s.layout.tail()) {
            Some(Tail::PolyT) => Some("forward"),
            Some(Tail::Tso) => Some("reverse"),
            None => None,
        });
    assert!(
        strand.is_none_or(|s| s == "forward" || s == "reverse"),
        "a strand word is forward or reverse"
    );
    strand
}

fn library_value(library: &LibraryEvidence) -> Value {
    let shares = [
        ("constant_region_share", library.constant_region),
        ("housekeeping_share", library.housekeeping),
        ("scaffold_share", library.scaffold),
        ("tag_start_cover", library.tag_cover),
    ];
    // Each is a count of reads over the reads examined; the writer checks the same again.
    for (key, share) in shares {
        assert!(share.is_finite(), "{key} is a finite share, not {share}");
        assert!((0.0..=1.0).contains(&share), "{key} lies in 0..=1, not {share}");
    }
    let kind = library.call();
    let mut fields = vec![
        ("kind", text(kind.label())),
        ("gene_expression", Value::Bool(kind == LibraryType::GeneExpression)),
    ];
    fields.extend(shares.map(|(key, share)| (key, Value::Share(share))));
    Value::Object(fields)
}

fn strand_value(strand: &StrandEvidence) -> Value {
    let shares = [
        ("poly_a", strand.poly_a),
        ("poly_t", strand.poly_t),
        ("chrm_sense", strand.mito_forward),
        ("chrm_antisense", strand.mito_reverse),
        ("nuclear_sense", strand.nuclear_sense),
        ("nuclear_antisense", strand.nuclear_antisense),
    ];
    // Each is a count of reads over the reads sampled; the writer checks the same again.
    for (key, share) in shares {
        assert!(share.is_finite(), "{key} is a finite share, not {share}");
        assert!((0.0..=1.0).contains(&share), "{key} lies in 0..=1, not {share}");
    }
    Value::Object(Vec::from(shares.map(|(key, share)| (key, Value::Share(share)))))
}

/// `value` on one line: the pretty form with its line breaks and indents taken out. Strings never
/// hold a raw line break (the writer escapes it), so every break removed is layout.
pub fn compact(value: &Value) -> String {
    let pretty = pretty(value);
    let one_line: String = pretty.lines().map(str::trim_start).collect();
    one_line
}

/// `value` pretty-printed, without a trailing newline.
pub fn pretty(value: &Value) -> String {
    let mut out = String::with_capacity(OUTPUT_CAPACITY_BYTES);
    write_value(&mut out, value, 0);
    assert!(!out.is_empty(), "every value prints something");
    assert!(!out.ends_with('\n'), "the caller adds the final newline");
    out
}

/// Appends `value`'s display form; a `String` sink cannot fail.
fn push_display(out: &mut String, value: impl Display) {
    let before = out.len();
    let Ok(()) = write!(out, "{value}") else { unreachable!("writing to a String cannot fail") };
    // Every caller writes a number or an escape, never nothing.
    assert!(out.len() > before, "a displayed value adds text");
}

fn write_value(out: &mut String, value: &Value, depth: usize) {
    assert!(depth <= DEPTH_MAX, "JSON nested {depth} deep, past {DEPTH_MAX}");
    let before = out.len();
    match value {
        Value::Null => out.push_str("null"),
        Value::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
        Value::Integer(n) => push_display(out, n),
        Value::Share(x) => write_share(out, *x),
        Value::String(s) => write_string(out, s),
        Value::Array(items) => {
            write_items(out, depth, ('[', ']'), items.iter().map(|item| (None, item)));
        }
        Value::Object(fields) => {
            write_items(
                out,
                depth,
                ('{', '}'),
                fields.iter().map(|(key, item)| (Some(*key), item)),
            );
        }
    }
    assert!(out.len() > before, "every JSON value prints at least one character");
}

/// One array or object: `[]`/`{}` when empty, else one item a line, indented two spaces a level.
fn write_items<'a>(
    out: &mut String,
    depth: usize,
    (open, close): (char, char),
    items: impl ExactSizeIterator<Item = (Option<&'static str>, &'a Value)>,
) {
    // Items sit one level deeper, which must still fit.
    assert!(depth < DEPTH_MAX, "a container at depth {depth} leaves no room for its items");
    let object = open == '{';
    assert!(
        (open, close) == ('{', '}') || (open, close) == ('[', ']'),
        "brackets pair: {open}{close}"
    );
    let count = items.len();
    out.push(open);
    for (index, (key, item)) in items.enumerate() {
        assert_eq!(key.is_some(), object, "object items have keys, array items none");
        out.push_str(if index == 0 { "\n" } else { ",\n" });
        push_indent(out, depth + 1);
        if let Some(key) = key {
            write_string(out, key);
            out.push_str(": ");
        }
        write_value(out, item, depth + 1);
    }
    if count > 0 {
        out.push('\n');
        push_indent(out, depth);
    }
    out.push(close);
}

fn push_indent(out: &mut String, depth: usize) {
    assert!(depth <= DEPTH_MAX, "indent {depth} past {DEPTH_MAX}");
    let before = out.len();
    for _ in 0..depth {
        out.push_str("  ");
    }
    assert_eq!(out.len(), before + 2 * depth, "two spaces a level");
}

/// A string with `"`, `\` and control characters escaped; other characters as they are.
fn write_string(out: &mut String, s: &str) {
    let before = out.len();
    out.push('"');
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            c if u32::from(c) < 0x20 => push_display(out, format_args!("\\u{:04x}", u32::from(c))),
            c => out.push(c),
        }
    }
    out.push('"');
    assert!(out.len() >= before + s.len() + 2, "an escaped string never shrinks");
    let Some(written) = out.get(before..) else { unreachable!("{before} is a char boundary") };
    assert!(written.chars().all(|c| u32::from(c) >= 0x20), "control characters are escaped");
}

/// A share in its shortest round-trip digits, laid out like `serde_json`: plain decimals from
/// 1e-5 up, exponent form below (`2.5e-7`), whole numbers keep `.0`.
fn write_share(out: &mut String, x: f64) {
    assert!(x.is_finite(), "a share is a finite fraction, not {x}");
    assert!((0.0..=1.0).contains(&x), "a share lies in 0..=1, not {x}");
    if x == 0.0 {
        out.push_str("0.0");
        return;
    }
    // `{:e}` gives the shortest digits that round-trip, as `d.ddde<exp>`.
    let scientific = format!("{x:e}");
    let Some((mantissa, exponent)) = scientific.split_once('e') else {
        unreachable!("{{:e}} always writes an exponent: {scientific}")
    };
    let Ok(exponent) = exponent.parse::<i32>() else {
        unreachable!("{{:e}} exponent: {scientific}")
    };
    let digits: String = mantissa.chars().filter(char::is_ascii_digit).collect();
    assert!(!digits.is_empty(), "shortest digits of {x} are empty");
    assert!(!digits.starts_with('0'), "shortest digits of nonzero {x} start with 0: {digits}");
    write_digits(out, &digits, exponent);
}

/// Lays out `digits` (no leading zero) with the decimal point after `exponent + 1` of them.
fn write_digits(out: &mut String, digits: &str, exponent: i32) {
    let Ok(length) = i32::try_from(digits.len()) else { unreachable!("17 digits at most") };
    assert!(
        (1..=F64_DIGITS_MAX).contains(&length),
        "an f64 has 1 to 17 significant digits, not {length}"
    );
    let point = exponent + 1; // Digits before the decimal point.
    let zeros = |out: &mut String, count: i32| {
        for _ in 0..count {
            out.push('0');
        }
    };
    if (length..=PLAIN_POINT_MAX).contains(&point) {
        out.push_str(digits);
        zeros(out, point - length);
        out.push_str(".0");
    } else if (1..=PLAIN_POINT_MAX).contains(&point) {
        let Ok(split) = usize::try_from(point) else { unreachable!("point {point} is positive") };
        let (head, tail) = digits.split_at(split);
        assert!(!tail.is_empty(), "a point inside the digits leaves digits after it");
        push_display(out, format_args!("{head}.{tail}"));
    } else if (-PLAIN_LEADING_ZEROS_MAX..=0).contains(&point) {
        out.push_str("0.");
        zeros(out, -point);
        out.push_str(digits);
    } else {
        let (first, rest) = digits.split_at(1);
        out.push_str(first);
        if !rest.is_empty() {
            out.push('.');
            out.push_str(rest);
        }
        // A positive exponent carries its sign, as serde_json writes it: 1e+16, 2.5e-7.
        let exponent = point - 1;
        push_display(out, format_args!("e{}{exponent}", if exponent > 0 { "+" } else { "" }));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn share(x: f64) -> String {
        let mut out = String::new();
        write_share(&mut out, x);
        out
    }

    fn digits(digits: &str, exponent: i32) -> String {
        let mut out = String::new();
        write_digits(&mut out, digits, exponent);
        out
    }

    #[test]
    fn shares_print_like_serde_json() {
        let cases = [
            (0.0, "0.0"),
            (1.0, "1.0"),
            (0.5, "0.5"),
            (0.1234, "0.1234"),
            (1e-5, "0.00001"),
            (3e-5, "0.00003"),
            (1.5e-5, "0.000015"),
            (1e-6, "1e-6"),
            (2.5e-7, "2.5e-7"),
            (1.0 / 3.0, "0.3333333333333333"),
            (0.02, "0.02"),
        ];
        for (x, expected) in cases {
            assert_eq!(share(x), expected, "{x:e}");
        }
    }

    #[test]
    fn digits_lay_out_on_every_side_of_the_point() {
        assert_eq!(digits("123", 2), "123.0");
        assert_eq!(digits("1", 15), "1000000000000000.0");
        assert_eq!(digits("1", 16), "1e+16");
        assert_eq!(digits("15", 16), "1.5e+16");
        assert_eq!(digits("123456789", 5), "123456.789");
        assert_eq!(digits("25", -1), "0.25");
        assert_eq!(digits("1", -5), "0.00001");
        assert_eq!(digits("1", -6), "1e-6");
    }

    #[test]
    fn strings_escape_and_trees_indent() {
        let value = Value::Object(vec![
            ("a", text("q\"b\\n\n\u{1}é")),
            ("e", Value::Array(vec![])),
            ("o", Value::Object(vec![])),
            ("l", Value::Array(vec![Value::Integer(7), Value::Null, Value::Bool(true)])),
        ]);
        let expected = "{\n  \"a\": \"q\\\"b\\\\n\\n\\u0001é\",\n  \"e\": [],\n  \"o\": {},\n  \
                        \"l\": [\n    7,\n    null,\n    true\n  ]\n}";
        assert_eq!(pretty(&value), expected);
    }

    /// The published schema for this major version exists and admits this exact version, so a
    /// major bump cannot ship without its schema file.
    #[test]
    fn the_schema_file_matches_the_schema_version() {
        let (major, minor) = SCHEMA_VERSION.split_once('.').unwrap_or_default();
        assert!(major.parse::<u32>().is_ok(), "major version '{major}' is a number");
        assert!(minor.parse::<u32>().is_ok(), "minor version '{minor}' is a number");
        let file = format!("{}/schema/report.v{major}.json", env!("CARGO_MANIFEST_DIR"));
        let schema = std::fs::read_to_string(&file).unwrap_or_default();
        assert!(!schema.is_empty(), "{file} is missing or empty");
        let pattern = format!("\"pattern\": \"^{major}\\\\.[0-9]+$\"");
        assert!(schema.contains(&pattern), "{file} pins schema_version to {major}.x");
    }
}
