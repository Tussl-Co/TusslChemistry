# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Map one free-text sentence to a 10x Genomics layout, by the fixed specification (see README.md).

Standard library only. Public functions:
  is_candidate(sentence) -> bool        the 10x-wording gate (the builder applies it; label() does not)
  guide_numbers(sentence) -> list[str]  10x user-guide numbers (CG0xxxxx), recorded verbatim, never mapped
  analyse(sentence) -> (mapped | None, why)
  label(sentence) -> str | None         rules 1-9 only

Return values of label():
  exact:   "10x 3' v2", "10x 3' v3/v3.1", "10x 3' v4 (GEM-X)", "unsupported (3' v1)",
           "10x 5' v1/v2", "10x 5' v3 (GEM-X)", "10x Multiome GEX", "10x Flex v1", "unclaimed (Flex v2)"
  partial: "3' unversioned", "5' unversioned"
  "ambiguous" (two or more different layouts in one sentence)
  None (no statement)
"""
import re

I = re.I
CONTEXT = re.compile(r"\b10x\b|chromium|single[- ]cell\s*(3|5)|gem-?x|v\(d\)j|immune profiling|multiome|fixed rna|\bflex\b", I)
GUIDE = re.compile(r"\bCG0\d{5}\b", I)

# Rule 2: other vendors' kits. "Parse" and "Fluent" are case-sensitive (they are also English words).
VENDOR_CS = re.compile(r"\bParse\b|\bFluent\b")
VENDOR_CI = re.compile(r"evercode|\bBD\s+Rhapsody\b|\brhapsody\b|\bscale\s?bio\b|dnbelab|pipseq|singleron|gexscope|ddseq|icell8", I)

# Rule 1: masks. Only the version token is replaced by <SW>; tool and sequencer names stay.
TOOLS = (r"cell\s?ranger(?:-arc|-atac|-dna)?|space\s?ranger|long\s?ranger|starsolo|star|seurat|bcl2fastq|mkfastq|kallisto|"
         r"bustools|alevin(?:-fry)?|salmon|scanpy|scvelo|velocyto|cellbender|souporcell|harmony|scvi-tools|scrublet|"
         r"doubletfinder|monocle|signac|archr|macs2|python|r|bioconductor|loupe(?:\s+browser)?")
KITWORD = r"(?:chromium|gem|kit|kits|reagent|reagents|chemistry|immune|profiling|multiome|flex|genomics)\b"
MASK_TOOL = re.compile(
    rf"(?P<pre>\b(?:{TOOLS})\b(?:\.gz)?(?:(?:\s+|\s*\()(?!{KITWORD})[A-Za-z][A-Za-z\-]*,?){{0,4}}?[\s(,:]*)"
    r"(?P<ver>v\.?\s?\d[\w.\-]*|\d+\.\d[\w.\-]*|version\s*\d[\w.\-]*)", I)
PRIMES = "'’′ʹ‘´"      # ' ’ ′ ʹ ‘ ´
MASK_SEQ = re.compile(
    rf"(?P<pre>\b(?:novaseq|nextseq|hiseq|miseq|iseq)\b(?:\s+(?!\S*(?:10x|chromium|gem|immune|multiome|flex|single|prime|[{PRIMES}]))\S+){{0,6}}?\s+)"
    r"(?P<ver>v\d[\w.\-]*)", I)
MASK_3PART = re.compile(r"(?<![\w.])v?\.?\d+\.\d+\.\d+[\w.\-]*", I)


def mask(s):
    s = MASK_TOOL.sub(lambda m: m.group("pre") + "<SW>", s)
    s = MASK_SEQ.sub(lambda m: m.group("pre") + "<SW>", s)
    return MASK_3PART.sub("<SW>", s)


# Rule 3: the end.
FIVE = re.compile(rf"(?<![\d.,])5\s?[{PRIMES}]|\b5[\s-]?prime\b|\bfive[\s-]prime\b|\b5p\b|immune\s+profiling|"
                  r"single[- ]cell\s+v\(d\)j|v\(d\)j\s+(?:reagent|enrichment)\s+kit", I)
THREE = re.compile(rf"(?<![\d.,])3\s?[{PRIMES}]|\b3[\s-]?prime\b|\bthree[\s-]prime\b|\b3p\b", I)

# Rule 4: the version. v1.1 and v3.1 are tried before v1 and v3; a trailing ".0" is allowed; other v<d>.<d> are rejected.
VER = re.compile(r"(?<![\w.])(?:v\.?\s?|version\s?)(1\.1|3\.1|[1-4])(?:\.0)?(?![\w]|\.\d)", I)
GEMX = re.compile(r"\bgem-?x\b", I)

# Rules 5 and 6.
MULTI = re.compile(r"multiome|\barc\b", I)
ATAC = re.compile(r"atac", I)
GEX_WORD = re.compile(r"gene\s+expression|\bgex\b|rna|transcript", I)
FLEX2 = re.compile(r"gem-?x\s+flex\b|\bflex\s+v2\b|\bflex\s+gene\s+expression\s+v2\b", I)
FLEX1 = re.compile(r"\bflex\b|fixed\s+rna\s+profiling", I)

TABLE = {
    ("3", "2"): "10x 3' v2", ("3", "3"): "10x 3' v3/v3.1", ("3", "3.1"): "10x 3' v3/v3.1",
    ("3", "4"): "10x 3' v4 (GEM-X)", ("3", "gemx"): "10x 3' v4 (GEM-X)", ("3", "1"): "unsupported (3' v1)",
    ("5", "1"): "10x 5' v1/v2", ("5", "1.1"): "10x 5' v1/v2", ("5", "2"): "10x 5' v1/v2",
    ("5", "3"): "10x 5' v3 (GEM-X)", ("5", "gemx"): "10x 5' v3 (GEM-X)",
}
EXACT = set(TABLE.values()) | {"10x Multiome GEX", "10x Flex v1", "unclaimed (Flex v2)"}
PARTIAL = {"3' unversioned", "5' unversioned"}


def is_candidate(sentence):
    return bool(CONTEXT.search(sentence))


def guide_numbers(sentence):
    return [g.upper() for g in GUIDE.findall(sentence)]


def analyse(sentence):
    """Return (mapped, why). mapped is None when the sentence makes no statement."""
    v = VENDOR_CS.search(sentence) or VENDOR_CI.search(sentence)
    if v:
        return None, f"skipped: other vendor ({v.group(0)})"
    s = mask(sentence)
    layouts = set()
    flex2 = FLEX2.search(s)
    if flex2:
        layouts.add("unclaimed (Flex v2)")
        s_ver = FLEX2.sub(" ", s)          # "GEM-X Flex" / "Flex v2" is not an end version
    else:
        s_ver = s
        if FLEX1.search(s):
            layouts.add("10x Flex v1")
    if MULTI.search(s) and not (ATAC.search(s) and not GEX_WORD.search(s)):
        layouts.add("10x Multiome GEX")
    ends = [e for e, rx in (("3", THREE), ("5", FIVE)) if rx.search(s)]
    versions = {m.group(1) for m in VER.finditer(s_ver)}
    if GEMX.search(s_ver):
        versions.add("gemx")
    for e in ends:
        got = {TABLE[(e, x)] for x in versions if (e, x) in TABLE}
        layouts |= got or {f"{e}' unversioned"}
    if not layouts:
        return None, "no end" + (f"; versions {sorted(versions)}" if versions else "")
    if len(layouts) > 1:
        return "ambiguous", "; ".join(sorted(layouts))
    return layouts.pop(), ""


def label(sentence):
    return analyse(sentence)[0]


TESTS = [
    ("Chromium Single Cell 5′ v2", "10x 5' v1/v2"),
    ("5′ v2", "10x 5' v1/v2"),
    ("5-prime v2", "10x 5' v1/v2"),
    ("five prime v1.1", "10x 5' v1/v2"),
    ("Chromium Single Cell Immune Profiling Next GEM v2", "10x 5' v1/v2"),
    ("Immune Profiling v2", "10x 5' v1/v2"),
    ("GEM-X 3' v4", "10x 3' v4 (GEM-X)"),
    ("GEM-X Universal 3’ Gene Expression", "10x 3' v4 (GEM-X)"),
    ("GEM-X 5′ v3", "10x 5' v3 (GEM-X)"),
    ("3' v3.1 HT", "10x 3' v3/v3.1"),
    ("Chromium Next GEM Single Cell 3' LT v3.1", "10x 3' v3/v3.1"),
    ("Cell Ranger v3.0.2 and 3' v2", "10x 3' v2"),
    ("Data were processed with Cell Ranger v3 from Single Cell 3' Reagent Kit v2 libraries", "10x 3' v2"),
    ("Cell Ranger using Chromium Single Cell 3' v2", "10x 3' v2"),
    ("NovaSeq 6000 with Chromium 3' v3", "10x 3' v3/v3.1"),
    ("Sequenced on NovaSeq 6000 S4 Reagent Kit v1.5 after 10x 3' library prep", "3' unversioned"),
    ("v2 or v3", None),
    ("Single Cell 3' v2 or v3", "ambiguous"),
    ("Parse Evercode WT v2", None),
    ("Parse Evercode WT v2, Chromium 3' v2 not used", None),
    ("BD Rhapsody 3' WTA with 10x comparison v3", None),
    ("Chromium Single Cell 3' Reagent Kits v1", "unsupported (3' v1)"),
    ("Chromium Single Cell 3' Library & Gel Bead Kit", "3' unversioned"),
    ("10x Genomics 5' kit", "5' unversioned"),
    ("Chromium Single Cell V(D)J Reagent Kit", "5' unversioned"),
    ("V(D)J sequences were assembled with 10x software", None),
    ("Chromium Next GEM Single Cell Multiome ATAC + Gene Expression", "10x Multiome GEX"),
    ("Chromium Single Cell ATAC v1.1 with multiome-style nuclei", None),
    ("Chromium Fixed RNA Profiling Reagent Kit", "10x Flex v1"),
    ("10x Flex Gene Expression", "10x Flex v1"),
    ("GEM-X Flex Gene Expression", "unclaimed (Flex v2)"),
    ("10x Flex v2", "unclaimed (Flex v2)"),
    ("10x Genomics v3 chemistry", None),
    ("10x 3 prime version 3", "10x 3' v3/v3.1"),
    ("10x 3p v3.0 libraries", "10x 3' v3/v3.1"),
    ("10x 3' v2.5", "3' unversioned"),
    ("10x 3' and 5' libraries", "ambiguous"),
    ("Single Cell 3ʹ v3 and 5´ v2", "ambiguous"),
    ("libraries made per CG000204", None),
    ("Seurat v4 analysis of 10x 3‘ v3 data", "10x 3' v3/v3.1"),
]

if __name__ == "__main__":
    bad = 0
    for s, want in TESTS:
        got = label(s)
        ok = got == want
        bad += not ok
        print(("PASS" if ok else "FAIL"), repr(s), "->", repr(got), "" if ok else f"(want {want!r})")
    # gate checks: the context gate is separate from label()
    gate = [("5′ v2", False), ("Chromium Single Cell 5′ v2", True), ("10x Genomics", True),
            ("10xGenomics", False), ("Single Cell 3' v3", True), ("Immune Profiling", True), ("CG000204", False)]
    for s, want in gate:
        ok = is_candidate(s) == want
        bad += not ok
        print(("PASS" if ok else "FAIL"), "is_candidate", repr(s), "->", is_candidate(s))
    ok = analyse("Parse Evercode WT v2, Chromium 3' v2 not used")[1].startswith("skipped: other vendor")
    bad += not ok
    print(("PASS" if ok else "FAIL"), "vendor skip reason recorded")
    ok = guide_numbers("per user guide CG000204 Rev D") == ["CG000204"]
    bad += not ok
    print(("PASS" if ok else "FAIL"), "guide number recorded")
    n = len(TESTS) + len(gate) + 2
    print(f"{n - bad}/{n} passed")
    raise SystemExit(1 if bad else 0)
