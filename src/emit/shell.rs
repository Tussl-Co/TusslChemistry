// Copyright 2026 Tussl Ltd
// SPDX-License-Identifier: Apache-2.0

//! An emitted command as a POSIX shell script: a header naming the call and its validation, then
//! the command, one flag per line. Literal arguments are quoted here and nowhere else; a
//! placeholder becomes `"${NAME:?set NAME to <what>}"`, so the command stops with that sentence
//! when the variable is unset instead of running on an empty path.

use super::{Arg, Emitted, Placeholder};
use std::fmt::Write as _;

/// Bytes besides letters and digits an argument may hold unquoted: none means anything to the
/// shell.
const PLAIN: &[u8] = b"_-./:=,+@%";

/// `text` as one shell word: bare when every byte is plain, else in single quotes, with each
/// single quote written as `'\''` (close, escaped quote, reopen). Nothing inside single quotes
/// is special to the shell, so this round-trips any text without a NUL.
#[must_use]
pub fn shell_quote(text: &str) -> String {
    assert!(!text.contains('\0'), "a shell word cannot hold a NUL");
    let quoted = if !text.is_empty()
        && text.bytes().all(|b| b.is_ascii_alphanumeric() || PLAIN.contains(&b))
    {
        text.to_owned()
    } else {
        format!("'{}'", text.replace('\'', r"'\''"))
    };
    assert!(!quoted.is_empty(), "an empty argument is still a word ('')");
    quoted
}

/// `text` made safe for a `#` comment line: every control character (a newline above all, which
/// would end the comment and start a command) is written as its escape. Notes carry file names,
/// and a directory name may hold a newline.
fn comment(text: &str) -> String {
    let safe: String = text
        .chars()
        .flat_map(|c| -> Box<dyn Iterator<Item = char>> {
            if c.is_control() { Box::new(c.escape_default()) } else { Box::new(std::iter::once(c)) }
        })
        .collect();
    assert!(!safe.chars().any(char::is_control), "a comment holds no control character");
    assert!(safe.len() >= text.len(), "escaping only adds characters");
    safe
}

/// A placeholder as a shell word: the variable, or the sentence saying what to set.
fn placeholder_word(p: Placeholder) -> String {
    // The sentence sits inside "${...:?...}": it must not end the expansion or start another.
    assert!(!p.meaning.bytes().any(|b| b"}\"$`\\".contains(&b)), "{}: meaning", p.name);
    // A shell variable: a capital or `_` first, then capitals, digits or `_`.
    let first_ok = p.name.bytes().next().is_some_and(|b| b.is_ascii_uppercase() || b == b'_');
    let rest_ok = p.name.bytes().all(|b| b.is_ascii_uppercase() || b.is_ascii_digit() || b == b'_');
    assert!(first_ok && rest_ok, "{}: a shell variable name", p.name);
    format!("\"${{{}:?set {} to {}}}\"", p.name, p.name, p.meaning)
}

/// The script: comment lines for the header and notes, then the command.
#[must_use]
pub fn render_shell(e: &Emitted) -> String {
    let mut out = String::new();
    let _ = writeln!(
        out,
        "# tchem {}: {} (confidence {})",
        env!("CARGO_PKG_VERSION"),
        comment(e.layout),
        e.confidence
    );
    let _ = writeln!(out, "# validation: {}", comment(&e.validation));
    for p in &e.placeholders {
        let _ = writeln!(out, "# set {}: {}", p.name, comment(p.meaning));
    }
    for note in &e.notes {
        let _ = writeln!(out, "# note: {}", comment(note));
    }
    let header = out.lines().count();
    assert!(out.lines().all(|l| l.starts_with("# ")), "the header is comment lines only");
    let words: Vec<(bool, String)> = e
        .argv
        .iter()
        .map(|arg| match arg {
            Arg::Text(t) => (t.starts_with('-'), shell_quote(t)),
            Arg::Var(p) => (false, placeholder_word(*p)),
        })
        .collect();
    for (index, (flag, word)) in words.iter().enumerate() {
        // A flag starts a new line; its value follows on the same line.
        let separator = match index {
            0 => "",
            _ if *flag => " \\\n  ",
            _ => " ",
        };
        out.push_str(separator);
        out.push_str(word);
    }
    out.push('\n');
    // Only the command follows the header: one line per flag, each continuing the last.
    assert!(
        out.lines().skip(header + 1).all(|l| l.starts_with("  ")),
        "after the header, only the command's continuation lines"
    );
    assert!(out.lines().filter(|l| !l.starts_with('#')).count() >= 1, "a script has a command");
    assert!(!e.argv.is_empty(), "a command has a program");
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::Command;

    /// What `sh` makes of a quoted word: it must be the original text, whatever it holds.
    #[test]
    fn quoted_words_round_trip_through_the_shell() {
        let cases = [
            "plain/path_1.fastq.gz",
            "with space",
            "it's",
            "$HOME `ls` \"q\" \\ ; & | > < * ? ~ # !",
            "line\nbreak",
            "",
            "--soloCBwhitelist",
            "10x 3' v3/v3.1",
        ];
        for text in cases {
            let script = format!("printf %s {}", shell_quote(text));
            let out = Command::new("sh").arg("-c").arg(&script).output().unwrap();
            assert_eq!(String::from_utf8(out.stdout).unwrap(), text, "{script}");
        }
        assert_eq!(shell_quote("a/b.fq"), "a/b.fq", "plain words stay bare");
    }

    /// A note naming a directory with a newline in it stays a comment: run through sh, the text
    /// after the newline is not a command.
    #[test]
    fn notes_cannot_break_out_of_their_comment() {
        let dir = std::env::temp_dir().join(format!("tchem-shell-{}", std::process::id()));
        std::fs::create_dir_all(&dir).unwrap();
        let marker = dir.join("PWNED");
        let e = Emitted {
            pipeline: super::super::Pipeline::Kb,
            argv: vec![Arg::Text("true".into()), Arg::Text("--x".into())],
            placeholders: Vec::new(),
            notes: vec![format!("lane in /data/a\ntouch {}\n", marker.display())],
            layout: "10x 3' v3/v3.1",
            confidence: "HIGH",
            validation: "checked\rtouch x".into(),
        };
        let script = render_shell(&e);
        let out = Command::new("sh").arg("-c").arg(&script).current_dir(&dir).output().unwrap();
        let ran = marker.exists();
        std::fs::remove_dir_all(&dir).unwrap();
        assert!(out.status.success(), "{script}");
        assert!(!ran, "the note ran as a command:\n{script}");
        assert!(script.contains("\\n"), "the newline is written as an escape: {script}");
    }

    /// An unset placeholder stops the command with its sentence; a set one is its value.
    #[test]
    fn placeholders_stop_when_unset() {
        let p = Placeholder { name: "STAR_INDEX", meaning: "a STAR genome index directory" };
        let word = placeholder_word(p);
        let unset = Command::new("sh")
            .arg("-c")
            .arg(format!("printf %s {word}"))
            .env_remove("STAR_INDEX")
            .output()
            .unwrap();
        assert!(!unset.status.success());
        assert!(
            String::from_utf8_lossy(&unset.stderr)
                .contains("set STAR_INDEX to a STAR genome index directory")
        );
        let set = Command::new("sh")
            .arg("-c")
            .arg(format!("printf %s {word}"))
            .env("STAR_INDEX", "/idx dir")
            .output()
            .unwrap();
        assert_eq!(String::from_utf8(set.stdout).unwrap(), "/idx dir");
    }
}
