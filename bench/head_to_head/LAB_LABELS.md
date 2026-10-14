# Lab-stated 10x labels (independent of the reads)

What each held-out 10x run's own records (ENA fields, the linked GEO sample) say its chemistry is,
mapped by a rule fixed before any outcome was seen (`scripts/lab_labels/label10x.py`, self-test:
`python3 -I label10x.py`). Built 9 Oct 2026 from public metadata only: no read, no tool call
and no Arc label was looked at. Coverage uses an exact label as a run's truth, with referee
"the lab's own record" (`bench/coverage/scripts/sources.py`, `head_to_head_cases`); a run
without one keeps Arc's list, which counts as read-derived (bench/VALIDATION.md, the
independent-label rule).

Results: `results/lab_labels.tsv` (one row per run), `results/lab_label_statements.tsv` (every
candidate sentence, quoted, with its URL), `results/lab_label_reasons.tsv` (hand-written reasons
for each conflict), `results/lab_label_runs.tsv` (the runs asked about).

To rebuild (raw records kept under `$TUSSL_CHEMISTRY_DATA/lab_labels/raw`):

```sh
W=$(mktemp -d); cp results/lab_label_runs.tsv $W/input_runs.tsv
python3 -I scripts/lab_labels/fetch_ena.py $W && python3 -I scripts/lab_labels/fetch_geo.py $W
#   (or: ln -s $TUSSL_CHEMISTRY_DATA/lab_labels/raw $W/raw, to reuse the records fetched on 9 Oct)
mkdir $W/out && cp results/lab_label_reasons.tsv $W/out/reasons.tsv
python3 -I scripts/lab_labels/build.py $W $W/out   # labels.tsv and statements.tsv, byte-identical
```

The working notes of the build follow.

# Independent 10x chemistry labels from public metadata

Inputs: the raw ENA and GEO records already fetched by the audit in `../../audit/tenx/` (`raw/ena/`, `raw/geo/`,
`input_runs.tsv`; only the srx, srr and study columns are read). Nothing was refetched. No sequence data, no
tool binary, no calls or results file was opened. `score.py` is not imported; `label10x.py` is standard library only.

## Implementation decisions (declared pre-run, before the first build)

The specification in the brief is applied as written. Where it left a choice, the choice is:

- **Gate vs mapping.** `is_candidate()` is the 10x-wording gate (the brief's regex, case-insensitive). `label()`
  applies rules 1-9 only. The builder calls both.
- **Sentences.** Each field is split into lines, then at `[.;]` + whitespace. Tabs collapse to one space in the quote.
- **Rule 1 masks** replace only the version token with `<SW>`; tool and sequencer names stay. Tools: Cell Ranger
  (incl. -arc/-atac/-dna), Space Ranger, Long Ranger, STARsolo, STAR, Seurat, bcl2fastq, mkfastq, kallisto, bustools,
  alevin(-fry), salmon, scanpy, scvelo, velocyto, CellBender, souporcell, Harmony, scvi-tools, Scrublet,
  DoubletFinder, Monocle, Signac, ArchR, MACS2, Python, R, Bioconductor, Loupe (Browser). A tool may be followed by
  up to 4 alphabetic words and `(`/`,`/`:` before its version; those words may not be 10x kit words (chromium, gem,
  kit, reagent, chemistry, immune, profiling, multiome, flex, genomics). Sequencers: NovaSeq, NextSeq, HiSeq, MiSeq,
  iSeq followed by up to 6 tokens and `v<n>`; those tokens may not contain 10x, chromium, gem, immune, multiome,
  flex, single, prime or a prime mark. Any `v?d+.d+.d+` is masked.
- **Rule 2 vendors:** `Parse` and `Fluent` case-sensitive whole words; case-insensitive Evercode, BD Rhapsody,
  Rhapsody, ScaleBio/Scale Bio, DNBelab, PIPseq, Singleron, GEXSCOPE, ddSEQ, ICELL8. The skip is written in
  `mapped` as `skipped: other vendor (...)`.
- **Rule 3 ends.** A prime mark may follow the digit after one optional space; the digit may not follow a digit, `.`
  or `,`. "5 prime"/"5-prime"/"5prime", "five prime"/"five-prime", whole-word `5p`; the same for 3. V(D)J kit name =
  `single[- ]cell V(D)J` or `V(D)J (Reagent|Enrichment) Kit`. Bare "V(D)J" is not an end.
- **Rule 4 versions.** `v`, `v.`, `v ` or `version` + one of 1.1, 3.1, 1, 2, 3, 4 (the two-part ones tried first),
  optionally followed by `.0`; any other `v<d>.<d>` (v2.5, v2.1) is not a version. "GEM-X"/"GEMX" is a version.
  A version that has no entry for the stated end (e.g. 5' v4) is ignored, so that end is "unversioned".
- **Several ends or versions in one sentence**: every (end, version) pair is mapped; two or more different results
  (including "3' unversioned" with "5' unversioned", or a Flex/Multiome layout with an end) is `ambiguous`.
  A bare "v2 or v3" with no end is no statement (rule 7), not ambiguous.
- **Rule 5.** `multiome` or whole-word `ARC`. ATAC-only wording = the sentence mentions ATAC and none of
  "gene expression", GEX, RNA, transcript; then rule 5 does not fire.
- **Rule 6.** Flex v2 wording ("GEM-X Flex", "Flex v2", "Flex Gene Expression v2") is removed before versions are
  read, so its GEM-X / v2 do not act as an end version.
- **Rule 9.** CG0xxxxx numbers make a sentence a recorded row even without 10x wording; they never map.
  The numbers stay visible in the quote, and `mapped` reads `none (guide CG0...)` when they are the only thing.
- **Exact layouts** are the 6 table values plus `10x Multiome GEX`, `10x Flex v1`, `unsupported (3' v1)` and
  `unclaimed (Flex v2)`. Partial = `3' unversioned`, `5' unversioned`. `ambiguous` and no-statement never count.
- **statements.tsv** has one row per candidate sentence (10x wording or a CG number), stating or not, including
  study_title rows (which never count). `mapped` is one of: a layout, a partial, `ambiguous (<layouts>)`,
  `none`, `none (<why>)`, `skipped: other vendor (<word>)`. `mapped` is computed on the full sentence; a sentence
  over 300 characters is quoted as a 300-character window that starts 40 characters before a decisive match.
(After run 3: the window is the first one that keeps the label and every CG number.) The builder counts
  rows where `label(quote) != label(sentence)` and reports them below.
- **labels.tsv.** `n_statements` = sentences (not fields) in the run's own fields (study_title excluded) whose
  mapped value is exact or partial. Label precedence: one exact value -> that value; two or more exact -> `conflict`;
  no exact, partials all one end -> `end_only (3')` / `end_only (5')`; else `not stated`. Partials are ignored once
  an exact value exists; runs where a partial disagrees with the exact label, and runs where partials disagree with
  no exact value, are counted below.
- **Linkage** as in the audit: ENA fields experiment_title, sample_title, library_name,
  library_construction_protocol, sample_alias, study_title, experiment_alias of the input SRR's read_run row; GSM from
  sample_alias / experiment_alias / experiment_title; GEO keys title, description, extract_protocol, characteristics,
  data_processing (each `!Sample_...` line is its own field occurrence). Raw files are read as UTF-8.

## Commands (Python 3, always `-I`, run from /tmp)

```
E=<scratch>/exec/tenx
A=<scratch>/audit/tenx
python3 -I $E/label10x.py            # self-test: 49/49 passed (40 label() cases, 7 gate cases, 2 checks)
python3 -I $E/build.py $A $E         # statements.tsv labels.tsv; reads reasons.tsv for conflict_reason; output in build_output.txt
```

The build was run four times. Run 1 was before the conflict review. Run 2 came after `reasons.tsv` was written.
Run 3 added a check that every CG number in a long sentence is kept in its quote; one was lost. Run 4 changed the
quote window only: the builder now takes the first window, from the earliest decisive match onward, that keeps the
label and every CG number. `label10x.py` did not change at any step. All four builds give the same mapped values
and labels.

## Conflict review (hand, quoted sentences only)

7 conflict runs in 4 studies. Every conflict is a real disagreement in the text, and none comes from a misreading
by the rule. In each case, a shared, study-wide protocol describes two or more assays. For PRJNA1208088 (4 runs),
the counted 3' side comes from the GEO description "10X 3prime v3.1". It is read through the pre-run "3prime"
allowance (no separator), which is wider than the spec's literal "3 prime"/"3-prime". The sentence that names both
assays maps to ambiguous and does not count. The one-line reasons are in
`reasons.tsv` and in labels.tsv `conflict_reason`.

## Rule changes (post conflict review)

None. No conflict came from a misreading, so the rule was not changed and nothing was re-run for a rule change.

## Counts (labels.tsv, 465 runs)

| label | runs | distinct studies |
|---|---|---|
| not stated | 215 | 63 |
| 10x 3' v3/v3.1 | 78 | 55 |
| 10x 5' v1/v2 | 57 | 41 |
| end_only (5') | 50 | 16 |
| 10x 3' v2 | 35 | 26 |
| end_only (3') | 12 | 11 |
| conflict | 7 | 4 |
| 10x Multiome GEX | 7 | 2 |
| 10x 5' v3 (GEM-X) | 2 | 2 |
| 10x 3' v4 (GEM-X) | 2 | 1 |

statements.tsv: 1429 candidate-sentence rows. mapped: none 739, 10x 3' v3/v3.1 192, 5' unversioned 149,
10x 5' v1/v2 138, 10x 3' v2 72, ambiguous 46, 3' unversioned 41, 10x Multiome GEX 37, 10x 5' v3 (GEM-X) 8,
10x 3' v4 (GEM-X) 4, skipped: other vendor 2, 10x Flex v1 1.

Checks: no record missing (nothing refetched); 50 sentences over 300 characters, 0 quote windows whose label
differs from the full sentence, 0 with a CG number lost from the quote;
1 sentence with mojibake (`2‚®â150 cycle`, a read-length "x", no prime mark involved).

## Recorded, not changed (observed outside the conflict review)

- 6 runs (SRX27087122, -123, -128, -129, -132, -133) are labelled `10x Multiome GEX` because their experiment and
  sample titles start with "Multiome", while their protocol says "Chromium Single Cell 5' ver.3 GEM ... kit" and a
  V(D)J enrichment kit. "ver.3" is not a version token in the specification, so those sentences are
  `5' unversioned` and do not count against the exact label. The specification was applied as written.
- 2 vendor-skipped rows (SRX28550141) name Parse, Fluent and 10x Fixed RNA Profiling in one fixation sentence.

## Readings beyond the spec's literal tokens (all declared pre-run, listed here for review)

"3prime"/"5prime" with no separator; an optional space before a prime mark; a `.0` suffix on versions; `v.`/`v `
forms; a version with no table entry for that end (e.g. 5' v4) read as unversioned; "V(D)J Enrichment Kit" read
as a 5' kit name; masks remove the version token only, not the tool name; CG-only sentences recorded as rows; bare
"v2 or v3" = no statement, "3' v2 or v3" = ambiguous.

Sample of 5 statements.tsv rows: `sample_rows.txt`, `random.Random(20261009).sample(all 1429 rows, 5)`, ENA URLs
shortened for display only.
