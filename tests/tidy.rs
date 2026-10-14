// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! Tidy: the TigerStyle rules a compiler does not check, checked over this crate's own source on
//! every `cargo test` (after TigerBeetle's `tidy.zig`):
//!
//! - every line is at most 100 columns,
//! - every production function is at most 70 lines,
//! - every file opens with the copyright and Apache-2.0 licence lines,
//! - every production function holds at least two assertions (`assert!`, `assert_eq!`,
//!   `assert_ne!`, `debug_assert*!`, and `unreachable!` for a path an invariant rules out),
//!   except the few in `EXEMPT`, each with the reason an assertion there would be noise.
//!
//! Production means `src/` outside `#[cfg(test)]` modules; the line and header rules also cover
//! `tests/`. Functions are found by a small scanner that
//! skips strings, characters and comments when matching braces; it understands this crate's
//! code, not all of Rust, and fails loudly on anything it cannot follow.

// Test code: a panic is how a test reports failure (clippy.toml allows the same inside
// #[cfg(test)], which an integration test's helper functions are not marked with).
#![allow(clippy::unwrap_used, clippy::indexing_slicing, clippy::panic)]

use std::path::Path;

const COLUMNS_MAX: usize = 100;
const FUNCTION_LINES_MAX: usize = 70;
const FUNCTION_ASSERTIONS_MIN: usize = 2;
/// Source files read per tree; the crate has about 30.
const FILES_MAX: usize = 64;
/// Directory levels walked under `src/` (`src/signals/kmers.rs` is one level down).
const DEPTH_MAX: usize = 4;

/// Functions allowed fewer than two assertions: `path::function` (the file's path under `src/`
/// without `.rs`) and why. Each entry must name a
/// function that exists and still has fewer than two, so the list cannot go stale.
const EXEMPT: &[(&str, &str)] = &[
    ("cli/refusal::fmt", "Display plumbing: one write! per refusal variant"),
    ("emit/mod::fmt", "Display plumbing: one write! per refusal variant"),
    ("emit/mod::from_name", "lookup: a name to its pipeline, or none"),
    ("emit/mod::name", "lookup: a pipeline to its fixed name (tested: names_are_lower_case)"),
    ("emit/mod::unsupported", "constructor: the refusal from its parts"),
    (
        "emit/mod::strand",
        "one precondition (short reads); the refusal handles an unmeasured strand",
    ),
    ("remote::url", "guard: the refusal above is the check (tested: urls_stay_on_the_host)"),
    (
        "layouts/validation::label",
        "lookup: a level to its JSON word (tested: level_words_are_json_words)",
    ),
    ("remote::fmt", "Display plumbing: one write! per refusal variant"),
    ("remote::source", "Error plumbing: the reads refusal under a stream's, if any"),
    ("output/json::fmt", "Display plumbing: one write! per refusal variant"),
    ("main::fmt", "Display plumbing: writes the message as given"),
    ("lib::fmt", "Display plumbing: one write! per refusal variant"),
    ("long_read/detect::fmt", "Debug plumbing: writes the kit's name"),
    ("long_read/detect::eq", "identity: two kits are equal when they are the same layout"),
    ("long_read/search::base_code", "lookup table: a base's 2-bit code or none"),
    ("signals/structure::acgt", "predicate: whether every byte is one of four bases"),
    ("input/fastq::format_of", "wrapper: opens a file for decoded, which asserts the format"),
    ("signals/instrument::label", "lookup table: one fixed name per platform"),
    ("detect/decide::label", "lookup table: one fixed word per confidence level"),
    ("output/json::compact", "wrapper: pretty asserts its output; removing layout adds no break"),
    ("emit/mod::printable", "validator: refuses a control character; the refusal is the check"),
    (
        "emit/mod::distinct",
        "validator: refuses one file given as both reads; the refusal is the check",
    ),
    ("emit/mod::unsplit", "validator: refuses a value holding a comma; the refusal is the check"),
    ("emit/mod::stored", "classifier: the reader's format to Stored or a refusal"),
    ("emit/starsolo::read_command", "lookup: a pair of storages to its read command, or a refusal"),
    ("input/fastq::decoded", "classifier: three magic checks; the magics differ at compile time"),
    ("input/bam::record_name", "field parser: the name bytes handed to kept_name, which asserts"),
    ("signals/instrument::from_read_group", "lookup: a PL tag to its platform, or none"),
    ("signals/instrument::original", "parser: the archive's accession token dropped, if any"),
    (
        "signals/instrument::dnbseq_flowcell",
        "predicate: whether an id has a DNBSEQ flowcell's shape",
    ),
    ("signals/instrument::mgi", "shape parser: a match, then named(), which asserts the id"),
    ("signals/instrument::ultima", "shape parser: a match, then named(), which asserts the id"),
    ("signals/instrument::pacbio", "shape parser: a match, then named(), which asserts the id"),
    ("signals/instrument::nanopore", "shape parser: a match, then named(), which asserts the id"),
    ("output/json::confidence_word", "lookup table: one word per confidence level"),
    ("input/fastq::fmt", "Display plumbing: one write! per refusal variant"),
    ("input/whitelist::fmt", "Display plumbing: one write! per refusal variant"),
    ("input/fastq::source", "Error plumbing: the io error under a refusal, if any"),
    ("input/whitelist::source", "Error plumbing: the io error under a refusal, if any"),
    ("output/json::text", "one-line constructor: a &str becomes a JSON string value"),
    ("hash::add", "the Fx mixing step itself: any word and any state are valid"),
    ("hash::write_u32", "Hasher forwarding to add, on the hottest lookup path"),
    ("hash::write_u64", "Hasher forwarding to add, on the hottest lookup path"),
    ("hash::write_usize", "Hasher forwarding to add, on the hottest lookup path"),
    ("hash::finish", "Hasher forwarding: returns the state"),
    ("signals/kmers::complement", "lookup table: every byte has a complement (N for non-ACGT)"),
    ("detect/profile::base_index", "lookup table: every byte has an index (4 for non-ACGT)"),
    ("detect/score::two_bit_code", "lookup table: every byte has a 2-bit code (A for non-ACGT)"),
    ("input/whitelist::code", "lookup table on the hottest path: a base's 2-bit code or none"),
    ("lists/mod::file_name", "forwards to entry, which asserts the registry's order"),
    ("lists/mod::label", "forwards to entry, which asserts the registry's order"),
    ("lists/registry::url", "lookup table: one fixed address per host"),
    ("lists/mod::optional", "forwards to entry, which asserts; checks no public fallback"),
    ("layouts/mod::label", "lookup table: one fixed label per tail"),
    ("lists/mod::barcode_length", "forwards to entry, which asserts the registry's order"),
    ("signals/library::label", "lookup table: one fixed string per library type"),
    ("signals/strand::poly_call", "forwards to lean, which asserts its inputs and its result"),
    ("signals/strand::mito_call", "forwards to lean, which asserts its inputs and its result"),
    ("signals/strand::nuclear_call", "forwards to lean, which asserts its inputs and its result"),
];

/// Functions on the per-read or per-barcode path whose checks are `debug_assert`s: they run in the
/// debug build, which the gate runs on every parity oracle input, and not in a release build.
/// Each entry must name a function that has fewer than two release checks and at least two
/// checks in the debug build, so the list cannot go stale.
const DEBUG_ONLY: &[(&str, &str)] = &[
    ("detect/score::slot", "per-read list lookup"),
    ("detect/score::get", "per-read list lookup"),
    ("detect/score::mismatches", "per-read base comparison"),
    ("detect/score::fixed_mismatches", "per-read linker comparison"),
    ("detect/score::fixed_allowed", "per-read linker tolerance"),
    ("detect/score::tail_allowed", "per-read tail tolerance"),
    ("detect/score::misfit", "per-read offset search"),
    ("detect/score::tally_tail", "per-read tail test"),
    ("detect/score::pack", "per-read window packing"),
    ("input/whitelist::pack", "per-barcode packing"),
    ("input/whitelist::hit", "per-barcode list lookup"),
    ("input/whitelist::hit_by_search", "per-barcode substitution search"),
    ("input/whitelist::bucket_of", "per-barcode bucket index"),
    ("layouts/mod::expected", "per-read tail bases"),
    ("layouts/mod::len", "per-read segment length"),
    ("lists/mod::index", "per-read list index; both properties also proved at compile time"),
];

struct Function {
    file: String,
    name: String,
    line_first: usize,
    line_count: usize,
    assertions: usize,
    debug_assertions: usize,
}

/// Every source file opens with the copyright and licence lines, so a file copied out of the repo
/// still says whose it is and on what terms.
const HEADER: &str = "// Copyright 2026 Tussl Ltd\n// SPDX-License-Identifier: Apache-2.0\n";

#[test]
fn files_carry_the_licence_header() {
    let missing: Vec<String> = all_sources()
        .into_iter()
        .filter(|(_, source)| !source.starts_with(HEADER))
        .map(|(file, _)| file)
        .collect();
    assert!(missing.is_empty(), "files without the licence header:\n{}", missing.join("\n"));
}

#[test]
fn lines_fit_100_columns() {
    let mut wide = Vec::new();
    for (file, source) in all_sources() {
        for (index, line) in source.lines().enumerate() {
            if line.chars().count() > COLUMNS_MAX {
                wide.push(format!("{file}:{}: {} columns", index + 1, line.chars().count()));
            }
        }
    }
    assert!(wide.is_empty(), "lines over {COLUMNS_MAX} columns:\n{}", wide.join("\n"));
}

#[test]
fn functions_fit_70_lines() {
    let long: Vec<String> = functions()
        .iter()
        .filter(|f| f.line_count > FUNCTION_LINES_MAX)
        .map(|f| format!("{}:{}: {} is {} lines", f.file, f.line_first, f.name, f.line_count))
        .collect();
    assert!(long.is_empty(), "functions over {FUNCTION_LINES_MAX} lines:\n{}", long.join("\n"));
}

#[test]
fn functions_hold_two_assertions() {
    let functions = functions();
    let key = |f: &Function| format!("{}::{}", f.file.trim_end_matches(".rs"), f.name);
    let short: Vec<String> = functions
        .iter()
        .filter(|f| f.assertions < FUNCTION_ASSERTIONS_MIN)
        .filter(|f| !EXEMPT.iter().chain(DEBUG_ONLY).any(|(name, _)| *name == key(f)))
        .map(|f| format!("{}:{}: {} has {}", f.file, f.line_first, f.name, f.assertions))
        .collect();
    assert!(short.is_empty(), "functions under two assertions:\n{}", short.join("\n"));
    for (name, why) in EXEMPT {
        assert!(!why.is_empty(), "exemption {name} needs a reason");
        let matching: Vec<&Function> = functions.iter().filter(|f| key(f) == *name).collect();
        assert!(!matching.is_empty(), "exemption {name} names no function");
        let earned = matching.iter().all(|f| f.assertions >= FUNCTION_ASSERTIONS_MIN);
        assert!(!earned, "exemption {name} is stale: it now holds two assertions");
    }
    for (name, why) in DEBUG_ONLY {
        assert!(!why.is_empty(), "debug-only entry {name} needs a reason");
        assert!(!EXEMPT.iter().any(|(n, _)| n == name), "{name} is in both lists");
        let matching: Vec<&Function> = functions.iter().filter(|f| key(f) == *name).collect();
        assert!(!matching.is_empty(), "debug-only entry {name} names no function");
        let earned = matching.iter().all(|f| f.assertions >= FUNCTION_ASSERTIONS_MIN);
        assert!(!earned, "debug-only entry {name} is stale: it holds two release checks");
        let checked = matching.iter().all(|f| f.debug_assertions >= FUNCTION_ASSERTIONS_MIN);
        assert!(checked, "debug-only entry {name} has fewer than two checks in the debug build");
    }
}

/// Every layout name the bench runner, the gold scorer and the validate table expect is one the
/// code defines, so renaming a layout cannot leave them scoring against a name nothing produces.
#[test]
fn bench_names_real_layouts() {
    let defined: Vec<String> = sources("src/layouts")
        .iter()
        .flat_map(|(_, source)| quoted_after(source, "name: "))
        .collect();
    assert!(defined.len() >= 20, "found only {} layout names", defined.len());
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let read = |path: &str| std::fs::read_to_string(root.join(path)).unwrap();
    let runner = read("bench/heldout/scripts/run.py");
    let start = runner.find("PARSE = (").unwrap();
    let end = start + runner[start..].find("\n}\n").unwrap();
    let mut wanted: Vec<(String, &str)> = quoted(&runner[start..end])
        .into_iter()
        .filter(|(_, next)| *next != ':')
        .map(|(name, _)| (name, "bench/heldout/scripts/run.py"))
        .collect();
    let scorer = read("bench/vendor_gold/scripts/score.py");
    wanted.extend(quoted_after(&scorer, "return ").into_iter().map(|n| (n, "score.py")));
    let table = read("scripts/validate_kits.tsv");
    for row in table.lines().filter(|l| !l.starts_with('#') && !l.starts_with("name\t")) {
        let want = row.split('\t').nth(1).unwrap_or_default();
        let name = want.split("  (confidence").next().unwrap_or_default();
        wanted.push((name.to_string(), "scripts/validate_kits.tsv"));
    }
    assert!(wanted.len() >= 20, "found only {} expected names", wanted.len());
    let unknown: Vec<String> = wanted
        .iter()
        .filter(|(name, _)| !defined.contains(name))
        .map(|(name, file)| format!("{file}: {name}"))
        .collect();
    assert!(unknown.is_empty(), "names no layout defines:\n{}", unknown.join("\n"));
}

/// Every layout has a base sample in the chaos suite, so a layout cannot ship without its damaged
/// copies being checked for a wrong HIGH call.
#[test]
fn chaos_covers_every_layout() {
    let defined: Vec<String> = sources("src/layouts")
        .iter()
        .flat_map(|(_, source)| quoted_after(source, "name: "))
        .collect();
    assert!(defined.len() >= 20, "found only {} layout names", defined.len());
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let chaos = std::fs::read_to_string(root.join("bench/chaos/scripts/chaos.py")).unwrap();
    let start = chaos.find("BASES = [").unwrap();
    let bases = &chaos[start..start + chaos[start..].find("\n]\n").unwrap()];
    let missing: Vec<&String> =
        defined.iter().filter(|name| !bases.contains(&format!("\"{name}\", "))).collect();
    assert!(missing.is_empty(), "layouts with no chaos base: {missing:?}");
}

/// The coverage tally and the binary agree on the kit families, so the family record the tally
/// writes is the one the binary prints.
#[test]
fn tally_and_binary_name_the_same_families() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let read = |path: &str| std::fs::read_to_string(root.join(path)).unwrap();
    let block = |text: &str, start: &str, end: &str| {
        let from = text.find(start).unwrap();
        let mut names: Vec<String> = quoted(&text[from..from + text[from..].find(end).unwrap()])
            .into_iter()
            .map(|q| q.0)
            .collect();
        names.sort();
        names
    };
    let tally = block(&read("bench/coverage/scripts/tally.py"), "FAMILIES = {", "\n}\n");
    let binary = block(&read("src/layouts/validation.rs"), "pub const FAMILIES", "\n];");
    assert!(tally.len() >= 8, "found only {} family names in tally.py", tally.len());
    assert_eq!(tally, binary, "tally.py FAMILIES and validation.rs FAMILIES differ");
}

/// Double-quoted strings in `text`, each with the character after its closing quote.
fn quoted(text: &str) -> Vec<(String, char)> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(open) = rest.find('"') {
        let after = &rest[open + 1..];
        let Some(close) = after.find('"') else { break };
        let next = after[close + 1..].trim_start().chars().next().unwrap_or(' ');
        out.push((after[..close].to_string(), next));
        rest = &after[close + 1..];
    }
    out
}

/// The double-quoted string right after each occurrence of `marker` (a layout's `name: "..."`, a
/// scorer's `return "..."`); occurrences followed by anything else are skipped.
fn quoted_after(text: &str, marker: &str) -> Vec<String> {
    text.match_indices(marker)
        .filter_map(|(at, _)| {
            let rest = text[at + marker.len()..].strip_prefix('"')?;
            rest.find('"').map(|close| rest[..close].to_string())
        })
        .collect()
}

/// Every `.rs` file under `<crate>/<tree>`, as (path relative to the tree, text), in path order.
fn sources(tree: &str) -> Vec<(String, String)> {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join(tree);
    let mut out = Vec::new();
    let mut pending = vec![(root.clone(), 0usize)];
    while let Some((dir, depth)) = pending.pop() {
        assert!(depth <= DEPTH_MAX, "{} is nested past {DEPTH_MAX} levels", dir.display());
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                pending.push((path, depth + 1));
            } else if path.extension().is_some_and(|e| e == "rs") {
                let name = path.strip_prefix(&root).unwrap().to_string_lossy().into_owned();
                out.push((name, std::fs::read_to_string(&path).unwrap()));
            }
        }
    }
    out.sort();
    assert!(!out.is_empty(), "no sources under {tree}");
    assert!(out.len() <= FILES_MAX, "{} sources under {tree}", out.len());
    out
}

/// Every file the line and header rules cover: the program and its tests.
fn all_sources() -> Vec<(String, String)> {
    let mut out = sources("src");
    out.extend(sources("tests").into_iter().map(|(name, text)| (format!("tests/{name}"), text)));
    out
}

/// Every production function in the crate.
fn functions() -> Vec<Function> {
    let mut out = Vec::new();
    for (file, source) in sources("src") {
        let code = blank_comments_and_literals(&source);
        let lines: Vec<&str> = code.lines().collect();
        // The test module starts at `#[cfg(test)]` over an inline `mod name {`. A lone
        // `#[cfg(test)]` marks one test-only function (skipped below) or a module declared in its
        // own file (`mod tidy;`), which ends nothing here.
        let tests_first = (0..lines.len()).find(|&i| {
            let opens_module = lines.get(i + 1).is_some_and(|l| {
                let l = l.trim();
                l.starts_with("mod ") && l.ends_with('{')
            });
            lines[i].trim() == "#[cfg(test)]" && opens_module
        });
        let production = &lines[..tests_first.unwrap_or(lines.len())];
        let mut index = 0;
        while index < production.len() {
            let Some(name) = function_name(production[index]) else {
                index += 1;
                continue;
            };
            let Some(line_last) = body_end(production, index) else {
                index += 1; // A declaration without a body: a trait method signature.
                continue;
            };
            let test_only = index > 0 && production[index - 1].trim() == "#[cfg(test)]";
            if test_only {
                index = line_last + 1;
                continue;
            }
            let body = production[index..=line_last].join("\n");
            out.push(Function {
                file: file.clone(),
                name,
                line_first: index + 1,
                line_count: line_last - index + 1,
                assertions: count_assertions(&body),
                debug_assertions: count_debug_build_assertions(&body),
            });
            index = line_last + 1;
        }
    }
    assert!(out.len() > 50, "found only {} functions: the scanner lost its way", out.len());
    out
}

/// The function a line declares, if it declares one.
fn function_name(line: &str) -> Option<String> {
    let mut rest = line.trim_start();
    for prefix in ["pub(crate) ", "pub(super) ", "pub ", "const ", "unsafe "] {
        rest = rest.strip_prefix(prefix).unwrap_or(rest);
    }
    let rest = rest.strip_prefix("fn ")?;
    let name: String = rest.chars().take_while(|c| c.is_alphanumeric() || *c == '_').collect();
    (!name.is_empty()).then_some(name)
}

/// The line holding the closing brace of the function declared at `first`, or None when the
/// declaration ends in `;` before any body opens.
fn body_end(lines: &[&str], first: usize) -> Option<usize> {
    let mut depth = 0i64;
    let mut opened = false;
    for (index, line) in lines.iter().enumerate().skip(first) {
        for c in line.chars() {
            match c {
                '{' => {
                    depth += 1;
                    opened = true;
                }
                '}' => depth -= 1,
                ';' if !opened && depth == 0 => return None,
                _ => {}
            }
            assert!(depth >= 0, "unbalanced braces at line {}", index + 1);
        }
        if opened && depth == 0 {
            return Some(index);
        }
    }
    panic!("function at line {} never closes", first + 1)
}

fn count_assertions(body: &str) -> usize {
    let body = without_debug_only(body);
    let macros = ["assert!(", "assert_eq!(", "assert_ne!(", "unreachable!("];
    // Only checks that run in a release build count: `debug_assert!(` contains `assert!(`, so a
    // match preceded by `debug_` is not counted.
    let mut count = 0;
    for m in macros {
        for (at, _) in body.match_indices(m) {
            let debug = body.get(..at).is_some_and(|before| before.ends_with("debug_"));
            count += usize::from(!debug);
        }
    }
    count
}

/// Checks that run in the debug build: every assertion macro, `debug_` forms included.
fn count_debug_build_assertions(body: &str) -> usize {
    let macros = ["assert!(", "assert_eq!(", "assert_ne!(", "unreachable!("];
    // `debug_assert!(` contains `assert!(`, so each call, debug or not, is counted once.
    let count: usize = macros.iter().map(|m| body.matches(m).count()).sum();
    assert!(count >= count_assertions(body), "the debug build runs every release check");
    count
}

/// The body with each `#[cfg(debug_assertions)]` statement or block blanked: what it holds is
/// compiled out of a release build.
fn without_debug_only(body: &str) -> String {
    const MARK: &str = "#[cfg(debug_assertions)]";
    let mut out = String::with_capacity(body.len());
    let mut rest = body;
    while let Some(at) = rest.find(MARK) {
        let (kept, gated) = rest.split_at(at);
        out.push_str(kept);
        let end = statement_end(&gated[MARK.len()..]).map_or(gated.len(), |e| MARK.len() + e);
        out.push_str(&" ".repeat(end));
        rest = &gated[end..];
    }
    out.push_str(rest);
    assert_eq!(out.len(), body.len(), "blanking keeps the length");
    assert!(!out.contains(MARK), "every debug-only statement is blanked");
    out
}

/// The byte just past the statement that opens `text`: the `;` that ends it at depth 0, or the
/// `}` that closes its outermost block when no `;` follows at depth 0 first.
fn statement_end(text: &str) -> Option<usize> {
    let mut depth = 0usize;
    for (at, c) in text.char_indices() {
        match c {
            '{' | '(' | '[' => depth += 1,
            '}' | ')' | ']' => {
                depth = depth.checked_sub(1)?;
                if depth == 0 && c == '}' && text[..at].trim_start().starts_with('{') {
                    return Some(at + 1);
                }
            }
            ';' if depth == 0 => return Some(at + 1),
            _ => {}
        }
    }
    None
}

/// The source with comments, string and character literals replaced by spaces, keeping lines and
/// columns, so braces and macro names inside them are not mistaken for code.
fn blank_comments_and_literals(source: &str) -> String {
    let chars: Vec<char> = source.chars().collect();
    let mut out = String::with_capacity(source.len());
    let mut index = 0;
    let blank = |c: char| if c == '\n' { '\n' } else { ' ' };
    while index < chars.len() {
        let c = chars[index];
        let next = chars.get(index + 1).copied();
        if c == '/' && next == Some('/') {
            while index < chars.len() && chars[index] != '\n' {
                out.push(' ');
                index += 1;
            }
        } else if c == '"' || (c == 'r' && next == Some('"')) {
            let raw = c == 'r';
            out.push(' ');
            index += if raw { 2 } else { 1 };
            while index < chars.len() && chars[index] != '"' {
                if !raw && chars[index] == '\\' {
                    out.push(blank(chars[index]));
                    index += 1;
                }
                out.push(blank(chars[index]));
                index += 1;
            }
            out.push(' ');
            index += 1;
        } else if c == '\'' && is_char_literal(&chars[index..]) {
            let length = chars[index + 1..].iter().position(|&c| c == '\'').unwrap() + 2;
            let length =
                if chars[index + 1] == '\\' && chars[index + 2] == '\'' { 4 } else { length };
            out.extend(std::iter::repeat_n(' ', length));
            index += length;
        } else {
            out.push(c);
            index += 1;
        }
    }
    out
}

/// A `'` opening a character literal (`'a'`, `'\n'`, `'\''`, `'\u{08}'`), not a lifetime.
fn is_char_literal(chars: &[char]) -> bool {
    matches!(chars, ['\'', '\\', ..] | ['\'', _, '\'', ..])
}
