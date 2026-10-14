# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""The documentation's Simplified Technical English rules (docs/glossary.md), checked:

  1. a sentence that describes has at most 25 words; a step of a numbered list (a procedure) at
     most 20;
  2. a paragraph or a list item has at most 6 sentences;
  3. no word that the glossary replaces occurs in prose (REPLACED);
  4. the glossary defines each term once.

Prose is what remains after code blocks, inline code, tables, headings, HTML comments and link
targets are removed. Inline code counts as one word. A sentence ends at '.', '!' or '?' followed by
a space or the end of the text, except after an abbreviation (ABBREVIATIONS) or inside a number.

It does not test grammar, voice or meaning: a reviewer does.

Usage: python ste_check.py [file...]   (default: FILES)
"""

import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bench"))
from invariant import invariant  # noqa: E402

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
FILES = ["README.md", "CONTRIBUTING.md", "ROADMAP.md", "CHANGELOG.md", "SECURITY.md",
         "bench/VALIDATION.md", "docs/glossary.md", "docs/json.md", "docs/emit.md",
         "docs/validation.md"]
GLOSSARY = "docs/glossary.md"

WORDS_DESCRIPTION_MAX = 25
WORDS_PROCEDURE_MAX = 20
SENTENCES_PARAGRAPH_MAX = 6
FILE_BYTES_MAX = 1 << 20

# Word (lower case, whole word) -> what to write instead. Prose only: inside `code` they are names.
REPLACED = {
    "whitelist": "barcode list",
    "whitelists": "barcode lists",
    "held out": "held-out (an adjective: a held-out run)",
    "chemistry": "kit or layout (Cell Ranger's option is `--chemistry`, in code format)",
    "chemistries": "kits or layouts",
    "stage": "level",
    "stages": "levels",
    "vouch": "the barcode list agrees with",
    "vouches": "the barcode list agrees with",
    "checked": "tested, verified or compared (`checked` in code format is a level)",
    "validated": "tested or verified (`validated` in code format is a level)",
    "unvalidated": "`unvalidated`, in code format",
    "contradicted": "`contradicted` in code format for the level; otherwise 'the reads contradict'",
    "rule 2": "the design-study rule",
    "rule 7": "the design-data rule",
    "a4": "the independent-label rule",
    "e.g.": "for example",
    "i.e.": "that is",
    "etc.": "a complete list",
    "obviously": "(delete it)",
    "simply": "(delete it)",
}
# The glossary names the old rule names once, to map them to the new ones.
REPLACED_ALLOWED = {GLOSSARY: {"rule 2", "rule 7", "a4", "checked", "validated", "whitelist"}}

ABBREVIATIONS = {"e.g.", "i.e.", "etc.", "vs.", "cf.", "approx.", "no.", "fig."}

# (file, the sentence's first words, why it may be longer). Each must still match a long sentence.
EXEMPT = []


def prose_units(text):
    """(line, text, is_step) for each paragraph or list item of prose."""
    lines = text.split("\n")
    invariant(len(lines) >= 1, "a file has a line")
    units, block, start, step, fenced, comment = [], [], 0, False, False, False
    for number, line in enumerate(lines, 1):
        if line.lstrip().startswith("```"):
            fenced = not fenced
        if comment and "-->" in line:
            comment = False
            line = ""
        elif "<!--" in line and "-->" not in line:
            comment = True
        skip = fenced or comment or line.lstrip().startswith(("```", "|", "#", "<!--")) \
            or not line.strip()
        item = re.match(r"^\s*([-*]|(\d+)\.)\s+", line)
        if (skip or item) and block:
            units.append((start, " ".join(block), step))
            block = []
        if skip:
            continue
        if item:
            line, step = line[item.end():], item.group(2) is not None
        elif not block:
            step = False
        if not block:
            start = number
        block.append(line.strip())
    if block:
        units.append((start, " ".join(block), step))
    invariant(not fenced, "every code block is closed")
    return units


def plain(unit):
    """The unit with inline code as one word, link targets and emphasis removed."""
    text = re.sub(r"`[^`]*`", "CODE", unit)
    text = re.sub(r"\]\([^)]*\)", "]", text)
    text = re.sub(r"[*_\[\]]", "", text)
    invariant("`" not in text, f"inline code is closed: {unit[:60]!r}")
    return text


def sentences(text):
    out, words = [], []
    for word in text.split():
        words.append(word)
        end = re.search(r"[.!?][\"')]*$", word) is not None
        if end and word.lower().strip("\"'(") not in ABBREVIATIONS:
            out.append(words)
            words = []
    if words:
        out.append(words)
    invariant(sum(len(s) for s in out) == len(text.split()), "every word is in one sentence")
    return out


def replaced_words(text, allowed):
    found = []
    lower = text.lower()
    for word, instead in REPLACED.items():
        if word in allowed:
            continue
        pattern = r"(?<![\w-])" + re.escape(word) + (r"(?![\w-])" if word[-1].isalnum() else "")
        if re.search(pattern, lower):
            found.append(f"'{word}': write {instead}")
    return found


def glossary_terms(text):
    terms = []
    for line in text.split("\n"):
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if line.startswith("|") and len(cells) == 2 and cells[0] not in ("Term", "---"):
            terms.append(cells[0].strip("`").lower())
    invariant(len(terms) > 0, "the glossary has terms")
    return terms


def check(path, text, exempt_used):
    problems = []
    allowed = REPLACED_ALLOWED.get(path, set())
    for line, unit, step in prose_units(text):
        text_plain = plain(unit)
        found = sentences(text_plain)
        if len(found) > SENTENCES_PARAGRAPH_MAX:
            problems.append(f"{path}:{line}: {len(found)} sentences in one paragraph "
                            f"(at most {SENTENCES_PARAGRAPH_MAX})")
        limit = WORDS_PROCEDURE_MAX if step else WORDS_DESCRIPTION_MAX
        for words in found:
            if len(words) <= limit:
                continue
            sentence = " ".join(words)
            exempt = [e for e in EXEMPT if e[0] == path and sentence.startswith(e[1])]
            if exempt:
                exempt_used.add(exempt[0][:2])
                continue
            problems.append(f"{path}:{line}: {len(words)} words (at most {limit}): {sentence[:90]}")
        for problem in replaced_words(text_plain, allowed):
            problems.append(f"{path}:{line}: {problem}")
    return problems


def main(paths):
    invariant(len(paths) > 0, "files to check")
    problems, exempt_used = [], set()
    for path in paths:
        with open(os.path.join(ROOT, path)) as f:
            text = f.read(FILE_BYTES_MAX + 1)
        invariant(len(text) <= FILE_BYTES_MAX, f"{path}: under {FILE_BYTES_MAX} bytes")
        problems += check(path, text, exempt_used)
        if path == GLOSSARY:
            terms = glossary_terms(text)
            twice = sorted({t for t in terms if terms.count(t) > 1})
            problems += [f"{path}: '{t}' is defined more than once" for t in twice]
    if paths == FILES:
        stale = [e for e in EXEMPT if e[:2] not in exempt_used]
        problems += [f"{e[0]}: exemption no longer needed: {e[1]!r}" for e in stale]
    for problem in problems:
        print(problem)
    print(f"{len(paths)} files, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or FILES))
