// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! How far each layout's HIGH calls have been checked on runs that played no part in designing
//! it. The counts are `data/validation.tsv`, which `bench/coverage` writes from one run of the
//! binary over every held-out run the benchmarks have; the binary embeds that file rather than
//! keep a second copy of the numbers.
//!
//! What is counted is the precision of a HIGH call: of the held-out runs this layout was called
//! at HIGH, how many the run's own record confirms. MEDIUM and LOW calls are not vouched for.
//! The level comes from studies, not runs: runs from one study share a lab, a protocol and
//! often one library, so they are not independent trials (the rule of three assumes they are).

use crate::ReadKind;

/// The table, as the coverage tally wrote it: comment lines, a header, one row per layout and
/// read kind with any HIGH call checked. A layout without a row has none.
const TABLE: &str = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/validation.tsv"));
const COLUMNS: [&str; 8] = [
    "layout",
    "read_kind",
    "exact_runs",
    "exact_studies",
    "family_runs",
    "read_label_runs",
    "disputed",
    "wrong",
];
/// At most one row per layout and read kind; a longer table is not this table.
const ROWS_MAX: usize = 2 * super::LAYOUTS.len();
/// The family record, in the same columns (`layout` holds the family's name):
/// `data/validation_family.tsv`, written by the same tally.
const FAMILY_TABLE: &str =
    include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/data/validation_family.tsv"));
/// Kit families whose labels rarely name the version or bead type, though the reads separate
/// the layouts: their family record asks, of the HIGH calls naming any layout of the family,
/// how many runs are that family. Same names as `bench/coverage/scripts/tally.py` `FAMILIES`.
pub const FAMILIES: &[(&str, &[&str])] = &[
    (
        "BD Rhapsody",
        &[
            "BD Rhapsody (original V1 beads)",
            "BD Rhapsody (Enhanced Beads)",
            "BD Rhapsody (Enhanced V2/V3 beads)",
        ],
    ),
    ("Parse Evercode WT", &["Parse Evercode WT", "Parse Evercode WT v3"]),
    ("SeekGene SeekOne DD", &["SeekGene SeekOne DD 3'", "SeekGene SeekOne DD 5'"]),
    (
        "Singleron GEXSCOPE",
        &["Singleron GEXSCOPE V1", "Singleron GEXSCOPE V2", "Singleron GEXSCOPE V3"],
    ),
];
const FAMILY_ROWS_MAX: usize = 2 * FAMILIES.len();
/// The rule of three: no error in n independent trials bounds the error rate below 3/n at 95%
/// confidence. The levels are where that bound falls under 3% and 10%.
const RULE_OF_THREE: f64 = 3.0;
const STUDIES_ERROR_UNDER_3PCT: u32 = 100;
const STUDIES_VALIDATED: u32 = 30;

/// One layout's held-out record for one kind of reads.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
#[non_exhaustive]
pub struct Validation {
    /// HIGH calls on held-out runs whose own record names exactly this layout.
    pub held_out_high_calls: u32,
    /// The studies those runs come from.
    pub studies: u32,
    /// HIGH calls on held-out runs whose record names only the kit family (a BD run that does
    /// not say which beads): right at the family level, not counted towards the level.
    pub family_only_calls: u32,
    /// HIGH calls right only against a label read off the reads themselves (barcode-list
    /// matching, linker positions): the same kind of signal the tool uses, so not counted towards
    /// the level (bench/VALIDATION.md, section 4).
    pub read_label_calls: u32,
    /// HIGH calls whose label names another kit and only the reads' own barcode evidence
    /// supports the call: neither counted nor treated as errors, and reported.
    pub disputed: u32,
    /// HIGH calls on held-out runs that the run's own record contradicts.
    pub wrong: u32,
}

/// What a layout's record supports, from the number of independent studies.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
#[non_exhaustive]
pub enum Level {
    /// No held-out HIGH call checked at the exact layout yet.
    Unvalidated,
    /// Checked on fewer than 30 studies: the counts are stated, no rate is claimed.
    Checked,
    /// 30 studies or more without a confident error: the error rate is under 10%.
    Validated,
    /// 100 studies or more without a confident error: the error rate is under 3%.
    ErrorUnder3Pct,
    /// At least one held-out HIGH call was wrong.
    Contradicted,
}

impl Level {
    /// The level as the JSON report and the text report name it.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Level::Unvalidated => "unvalidated",
            Level::Checked => "checked",
            Level::Validated => "validated",
            Level::ErrorUnder3Pct => "error_under_3pct",
            Level::Contradicted => "contradicted",
        }
    }
}

impl Validation {
    /// The level the record earns; any confident error outranks every count.
    #[must_use]
    pub fn level(self) -> Level {
        assert!(self.studies <= self.held_out_high_calls, "more studies than runs: {self:?}");
        let level = if self.wrong > 0 {
            Level::Contradicted
        } else if self.studies >= STUDIES_ERROR_UNDER_3PCT {
            Level::ErrorUnder3Pct
        } else if self.studies >= STUDIES_VALIDATED {
            Level::Validated
        } else if self.studies >= 1 {
            Level::Checked
        } else {
            Level::Unvalidated
        };
        assert!(
            (level == Level::Unvalidated) == (self.wrong == 0 && self.studies == 0),
            "unvalidated exactly when nothing exact was checked: {self:?}"
        );
        level
    }

    /// The 95% upper bound on the error rate of a HIGH call (3 / studies), when no call was
    /// wrong and at least one study was checked.
    #[must_use]
    pub fn error_bound(self) -> Option<f64> {
        assert!(self.studies <= self.held_out_high_calls, "more studies than runs: {self:?}");
        let bound = (self.wrong == 0 && self.studies > 0)
            .then(|| (RULE_OF_THREE / f64::from(self.studies)).min(1.0));
        assert!(bound.is_none_or(|b| b > 0.0 && b <= 1.0), "a bound is a rate: {bound:?}");
        bound
    }
}

/// The table's rows, parsed once (a few dozen; the table is built into the binary).
fn rows() -> &'static [(&'static str, ReadKind, Validation)] {
    static ROWS: std::sync::OnceLock<Vec<(&'static str, ReadKind, Validation)>> =
        std::sync::OnceLock::new();
    let rows = ROWS.get_or_init(|| parse_table(TABLE, ROWS_MAX));
    assert!(!rows.is_empty(), "the embedded table has rows");
    assert!(rows.len() <= ROWS_MAX, "at most one row per layout and read kind");
    rows
}

/// The family table's rows, parsed once.
fn family_rows() -> &'static [(&'static str, ReadKind, Validation)] {
    static ROWS: std::sync::OnceLock<Vec<(&'static str, ReadKind, Validation)>> =
        std::sync::OnceLock::new();
    let rows = ROWS.get_or_init(|| parse_table(FAMILY_TABLE, FAMILY_ROWS_MAX));
    assert!(rows.len() <= FAMILY_ROWS_MAX, "at most one family row per family and read kind");
    assert!(
        rows.iter().all(|(_, _, v)| v.family_only_calls == 0),
        "a family record has no family-only calls: every right call is right for the family"
    );
    rows
}

/// Comment lines, the header, then at most `rows_max` rows.
fn parse_table(table: &'static str, rows_max: usize) -> Vec<(&'static str, ReadKind, Validation)> {
    let mut lines = table.lines().filter(|l| !l.starts_with('#'));
    let header = lines.next().unwrap_or_default();
    assert!(header.split('\t').eq(COLUMNS), "validation table header is {header:?}");
    let rows: Vec<_> = lines.take(rows_max + 1).map(parse_row).collect();
    assert!(rows.len() <= rows_max, "a validation table has more than {rows_max} rows");
    rows
}

/// The family `layout` belongs to, if any (`FAMILIES`).
#[must_use]
pub fn family_of(layout: &str) -> Option<&'static str> {
    assert!(!layout.is_empty(), "a layout has a name");
    let found = FAMILIES.iter().find(|(_, members)| members.contains(&layout)).map(|(f, _)| *f);
    assert!(
        FAMILIES.iter().filter(|(_, members)| members.contains(&layout)).count() <= 1,
        "{layout} belongs to one family at most"
    );
    found
}

/// The family record for `layout`'s family on `kind` reads (all zeros when none of its HIGH
/// calls has been checked); None for a layout in no family.
#[must_use]
pub fn of_family(layout: &str, kind: ReadKind) -> Option<(&'static str, Validation)> {
    let family = family_of(layout)?;
    let found = family_rows().iter().find(|(name, k, _)| *name == family && *k == kind);
    let record = found.map_or_else(Validation::default, |(_, _, v)| *v);
    let rows_for =
        family_rows().iter().filter(|(name, k, _)| *name == family && *k == kind).count();
    assert!(rows_for <= 1, "data/validation_family.tsv has {rows_for} rows for {family}");
    assert_eq!(record.family_only_calls, 0, "{family}: a family record has no family-only calls");
    Some((family, record))
}

/// The record for `layout` (its name in `LAYOUTS`) on `kind` reads; all zeros when none of its
/// HIGH calls has been checked on held-out runs.
#[must_use]
pub fn of(layout: &str, kind: ReadKind) -> Validation {
    assert!(!layout.is_empty(), "a layout has a name");
    let found = rows().iter().find(|(name, k, _)| *name == layout && *k == kind);
    let record = found.map_or_else(Validation::default, |(_, _, v)| *v);
    let rows_for = rows().iter().filter(|(name, k, _)| *name == layout && *k == kind).count();
    assert!(rows_for <= 1, "data/validation.tsv has {rows_for} rows for {layout}");
    record
}

/// One row of the embedded table. The table is built into the binary, so a malformed row is a
/// broken build, not untrusted input: it halts.
fn parse_row(line: &'static str) -> (&'static str, ReadKind, Validation) {
    let fields: Vec<&str> = line.split('\t').collect();
    assert_eq!(fields.len(), COLUMNS.len(), "data/validation.tsv row {line:?}");
    let number = |index: usize| -> u32 {
        let field = fields.get(index).copied().unwrap_or_default();
        field.parse().unwrap_or_else(|_| unreachable!("{line:?}: column {index} is {field:?}"))
    };
    let kind = match fields.get(1).copied() {
        Some("short") => ReadKind::Short,
        Some("long") => ReadKind::Long,
        other => unreachable!("{line:?}: read kind {other:?}"),
    };
    let validation = Validation {
        held_out_high_calls: number(2),
        studies: number(3),
        family_only_calls: number(4),
        read_label_calls: number(5),
        disputed: number(6),
        wrong: number(7),
    };
    assert!(validation.studies <= validation.held_out_high_calls, "{line:?}: studies > runs");
    (fields.first().copied().unwrap_or_default(), kind, validation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::layouts::LAYOUTS;

    /// A level word is a JSON enum value: lower case, digits and `_` only, so a pipeline can
    /// match it exactly, and no two levels share one.
    #[test]
    fn level_words_are_json_words() {
        let levels = [
            Level::Unvalidated,
            Level::Checked,
            Level::Validated,
            Level::ErrorUnder3Pct,
            Level::Contradicted,
        ];
        let words: Vec<&str> = levels.iter().map(|l| l.label()).collect();
        for word in &words {
            let plain =
                word.bytes().all(|b| b.is_ascii_lowercase() || b == b'_' || b.is_ascii_digit());
            assert!(!word.is_empty() && plain, "level word {word:?}");
        }
        let mut distinct = words.clone();
        distinct.sort_unstable();
        distinct.dedup();
        assert_eq!(distinct.len(), words.len(), "every level has its own word");
    }

    use std::fmt::Write as _;

    const README: &str = concat!(env!("CARGO_MANIFEST_DIR"), "/README.md");
    const BEGIN: &str = "<!-- coverage:begin -->\n";
    const END: &str = "<!-- coverage:end -->";
    const FAMILY_BEGIN: &str = "<!-- family:begin -->\n";
    const FAMILY_END: &str = "<!-- family:end -->";

    /// The README's coverage table, rendered from the embedded table with the binary's own levels
    /// and bounds, so the README cannot state a level the binary would not print.
    fn readme_table() -> String {
        let mut out = String::from(
            "| Layout | Reads | Level | Held-out HIGH calls (studies) | Error bound (95%) | \
             Kit family only | Read-derived label only | Disputed | Wrong |\n\
             |---|---|---|---|---|---|---|---|---|\n",
        );
        for (name, kind, v) in rows() {
            let kind = if *kind == ReadKind::Short { "short" } else { "long" };
            let bound = v.error_bound().map_or("-".to_owned(), |b| format!("{:.1}%", b * 100.0));
            let _ = writeln!(
                out,
                "| {name} | {kind} | {} | {} ({}) | {bound} | {} | {} | {} | {} |",
                v.level().label(),
                v.held_out_high_calls,
                v.studies,
                v.family_only_calls,
                v.read_label_calls,
                v.disputed,
                v.wrong
            );
        }
        out
    }

    /// The README's family table, rendered from the embedded family table.
    fn readme_family_table() -> String {
        let mut out = String::from(
            "| Kit family | Layouts | Reads | Level | Held-out HIGH calls (studies) | \
             Error bound (95%) | Read-derived label only | Disputed | Wrong |\n\
             |---|---|---|---|---|---|---|---|---|\n",
        );
        for (name, kind, v) in family_rows() {
            let members = FAMILIES.iter().find(|(f, _)| f == name).map_or(0, |(_, m)| m.len());
            let kind = if *kind == ReadKind::Short { "short" } else { "long" };
            let bound = v.error_bound().map_or("-".to_owned(), |b| format!("{:.1}%", b * 100.0));
            let _ = writeln!(
                out,
                "| {name} | {members} | {kind} | {} | {} ({}) | {bound} | {} | {} | {} |",
                v.level().label(),
                v.held_out_high_calls,
                v.studies,
                v.read_label_calls,
                v.disputed,
                v.wrong
            );
        }
        out
    }

    /// Replaces the text between `begin` and `end` in `text` with `table`.
    fn with_table(text: &str, begin: &str, end: &str, table: &str) -> String {
        let (before, rest) = text.split_once(begin).unwrap_or_else(|| panic!("no {begin:?}"));
        let (_, after) = rest.split_once(end).unwrap_or_else(|| panic!("no {end:?}"));
        format!("{before}{begin}{table}{end}{after}")
    }

    /// `tchem --help` gives the number of layouts: it must be the number that tchem knows.
    #[test]
    fn the_help_counts_the_layouts() {
        let usage = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/src/cli/usage.rs"));
        let said = format!("It knows {} layouts", super::super::LAYOUTS.len());
        assert!(usage.contains(&said), "src/cli/usage.rs does not say '{said}'");
        assert!(super::super::LAYOUTS.len() > 1, "a registry of layouts");
    }

    /// The README's tables are the embedded tables'; regenerate with `TCHEM_BLESS=1`.
    #[test]
    fn the_readme_states_what_the_binary_prints() {
        let readme = std::fs::read_to_string(README).unwrap_or_default();
        let want = with_table(&readme, BEGIN, END, &readme_table());
        let want = with_table(&want, FAMILY_BEGIN, FAMILY_END, &readme_family_table());
        if std::env::var_os("TCHEM_BLESS").is_some() {
            std::fs::write(README, &want).unwrap_or_else(|e| panic!("{README}: {e}"));
        }
        assert!(readme == want, "README.md's coverage tables are stale; rerun with TCHEM_BLESS=1");
    }

    /// Every family names real layouts, each in one family only, and every family row names a
    /// family once per read kind.
    #[test]
    fn families_name_layouts_and_rows_name_families() {
        for (family, members) in FAMILIES {
            assert!(members.len() >= 2, "{family}: a family has at least two layouts");
            for m in *members {
                assert!(LAYOUTS.iter().any(|l| l.name == *m), "{family}: unknown layout {m:?}");
                assert_eq!(family_of(m), Some(*family), "{m} belongs to {family}");
            }
        }
        let mut seen = Vec::new();
        for (name, kind, _) in family_rows() {
            assert!(FAMILIES.iter().any(|(f, _)| f == name), "unknown family {name:?}");
            assert!(!seen.contains(&(*name, *kind)), "{name} {kind:?} appears twice");
            seen.push((*name, *kind));
        }
        assert_eq!(family_of("10x 3' v3/v3.1"), None, "a layout outside every family");
        assert_eq!(of_family("10x 3' v3/v3.1", ReadKind::Short), None);
    }

    /// Every row names a layout, once per read kind, and long-read rows only for layouts the
    /// long-read detector offers.
    #[test]
    fn every_row_names_a_layout_once() {
        let mut seen = Vec::new();
        for line in TABLE.lines().filter(|l| !l.starts_with('#')).skip(1) {
            let (name, kind, _) = parse_row(line);
            let layout = LAYOUTS.iter().find(|l| l.name == name);
            assert!(layout.is_some(), "data/validation.tsv names unknown layout {name:?}");
            if kind == ReadKind::Long {
                assert!(layout.is_some_and(|l| l.long_read_kit.is_some()), "{name}: long row");
            }
            assert!(!seen.contains(&(name, kind)), "{name} {kind:?} appears twice");
            seen.push((name, kind));
            assert_eq!(of(name, kind), parse_row(line).2, "lookup finds the row");
        }
        assert!(!seen.is_empty(), "the table has rows");
    }

    #[test]
    fn levels_follow_studies_and_errors_outrank_counts() {
        let record = |runs, studies, wrong| Validation {
            held_out_high_calls: runs,
            studies,
            wrong,
            ..Validation::default()
        };
        assert_eq!(record(0, 0, 0).level(), Level::Unvalidated);
        assert_eq!(record(0, 0, 0).error_bound(), None);
        assert_eq!(record(49, 25, 0).level(), Level::Checked, "runs alone do not validate");
        assert_eq!(record(49, 30, 0).level(), Level::Validated);
        assert_eq!(record(400, 100, 0).level(), Level::ErrorUnder3Pct);
        assert_eq!(record(400, 100, 1).level(), Level::Contradicted);
        assert_eq!(record(400, 100, 1).error_bound(), None);
        assert_eq!(record(1, 1, 0).error_bound(), Some(1.0), "a bound is capped at 1");
        assert_eq!(record(60, 30, 0).error_bound(), Some(0.1));
        assert_eq!(of("no such layout", ReadKind::Short), Validation::default());
    }
}
