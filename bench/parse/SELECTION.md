# Parse catalogue: how it was extended on 9 Oct 2026

`results/catalogue.tsv` is append-only. The 147 rows drawn on 8 Oct, with their 9 design runs, are kept
exactly as drawn, because their runs had already been called. The 9 Oct draw (`scripts/select.py`,
whose docstring records each rule revision and its reason) adds 128 held-out rows from 67 studies
not already in the catalogue. Each of those studies names the Parse WT kit. The rules were revised
before any added study was run.

| kit | runs | studies |
|---|---|---|
| parse_wt (WT named, no single version; new kit, no design run) | 101 | 53 |
| parse_wt_v1 | 8 | 4 |
| parse_wt_v2 | 9 | 5 |
| parse_wt_v3 | 10 | 5 |

Held-out Parse runs are scored at the family level (`bench/heldout/scripts/run.py`, `PARSE`), so
a run whose protocol states no version tests the family record, not an exact layout.

The working notes of the 9 Oct draw and of the audit behind it follow.

## Parse Evercode WT: revised draw, 9 Oct 2026 (appendix to the catalogue)

`select.py` is a revised copy of `bench/parse/scripts/select.py`. Its docstring records each change,
the reason for it and the study-level evidence. `appended.tsv` holds only the rows of studies that
are not in `bench/parse/results/catalogue.tsv`. It uses that file's column order (kit, study, run,
r1_url, r2_url, split, evidence), and every row is `held_out`. The 147 existing rows and their
9 design runs do not change: the catalogue is append-only.

No added study has been run through tchem. The revision is made before any run, as VALIDATION.md
section 5, step 1 allows.

## Run

From `bench/parse/scripts/` once copied there (the import path is relative, as in
`bench/mgi/scripts/select.py`):

    python -I -B select.py ../results/catalogue.tsv <appended.tsv>

To run this copy from the scratch directory, the repository's `bench/heldout/scripts` was put on
`sys.path` before the module was loaded. Read lengths were fetched only for studies not already in
the catalogue (`catalogue.read_length`, first 64 KB), and only for the 2 runs of each that were
kept.

## Counts

| kit | runs | studies |
|---|---|---|
| parse_wt (no version stated; new kit, no design run) | 101 | 53 |
| parse_wt_v1 | 8 | 4 |
| parse_wt_v2 | 9 | 5 |
| parse_wt_v3 | 10 | 5 |
| **total** | **128** | **67** |

Check: the 128 (study, kit, run) triples are exactly those in `../proposed_added.tsv`, with none
missing and none extra. No study in `appended.tsv` is in the existing catalogue.

To use these rows:
- add `"parse_wt": PARSE` to `bench/heldout/scripts/run.py` `EXPECTED`, and to any coverage or
  tally code keyed on kit names;
- give `parse_wt` no design run (both Parse layouts already have theirs). Because of this, the
  appendix is not passed through `catalogue.drawn_sound`'s rule of one design run per kit.

## What the revision does not take

These were decided before the draw:
- the 11 studies that name Parse without a product;
- the 3 studies that only mention 10x (PRJNA1094867, PRJNA1169626, PRJNA993024), which
  `names_other_vendor` refuses;
- PRJNA1390207, refused by name: its first paired runs are bulk QIAseq.

Of the 67 studies, 66 are certainly WT. PRJNA1120008 is WT Mini v2 modified for scMPRA, so its
library is uncertain. In PRJEB65606 the barcode read looks filed as R1 (R1 86, R2 74).

## Side finding: 7 non-WT studies already in the catalogue (recorded decision)

ENA run metadata shows that 14 held-out rows in 7 catalogue studies are not Parse WT
gene-expression runs:

| study | rows | runs are | tchem (checked by the coordinator) |
|---|---|---|---|
| PRJNA1208204 | 2 | Parse CRISPR Detect barcode sublibrary | Parse at HIGH |
| PRJNA1304277 | 2 | "Mini_Sublibrary1/2_CRISPR" | Parse at HIGH |
| PRJNA1456597 | 2 | gRNA ClonMapper barcode capture sublibraries | Parse at HIGH |
| PRJNA1463191 | 2 | "Parse CRISPR guide RNA sublibrary" | Parse at HIGH |
| PRJNA1465771 | 2 | KDM1A CUT&RUN (GENOMIC) | no call |
| PRJNA1449947 | 2 | strategy OTHER, source GENOMIC | no call |
| PRJNA975472 | 2 | 10x Visium spatial | no call |

**Decision: these rows stay as they are.**
- The CRISPR guide sublibraries carry Parse's barcode read on R2. A HIGH call of Parse is right
  for the layout tchem reads.
- The Visium, CUT&RUN and genomic runs got no call. A run with no call is not a HIGH call, so it
  adds nothing to the precision record and no error.

For new studies, rule 4 in `select.py` (non-WT products refused) keeps guide libraries out, and
rule 5 (other vendors refused) keeps out studies like PRJNA975472.

## Audit: Why the Parse catalogue draws far fewer studies than ENA holds

Inputs: the audit lists 189 studies by text (157 by strict regex), from
`research/atlases/data/ena/parse_evercode.studies.txt`. Today's catalogue has 78 studies
(147 rows), and all 78 are audit studies. **111 audit studies are missing.** 84 of the 111 are
in the strict set.

Method: `analyse.py` mirrors `catalogue.draw` with `select.py`'s QUERY and `classify`. It uses
ENA as of 9 Oct 2026 (`cache/query_rows.json`). `select.py` passes no `shape` and does not call
`names_other_vendor`, so today neither one refuses any Parse study. Both are reported as notes and
are used only in the proposal. Read lengths come from `catalogue.read_length` (first 64 KB,
lengths only), for at most the first 2 paired runs of each study. Some extra run metadata (library
source/strategy, titles) was read from the ENA portal. No sequence content was read, `tchem` was
not run, and no results file was opened.

## Reasons (first failing step of today's draw)

| reason | studies |
|---|---|
| classify returned None | 58 |
| ... no kit version stated anywhere | 38 |
| ... a version is stated, but not beside "WT"/"whole transcriptome" (e.g. "Evercode V2 kit", "Version 2 WT Kit", "Kit (v1.3.1)" written so the regex misses it) | 19 |
| ... more than one WT version (PRJNA1255131: "WT + TCR v1" and "v2 WT") | 1 |
| ... version only beside "fixation" | 0 |
| not returned by QUERY | 28 |
| no paired run (2 FASTQ files) | 17 |
| passes today's rules, but absent from the 8 Oct draw | 8 |
| the protocol names another vendor (not applied to Parse today) | 0 as first reason; noted for 14 |
| shape refused (no shape for Parse today) | 0 |

Notes:
- **QUERY**: ENA does not return 28 studies for `"*Parse*Bio*"` or `"*Evercode*"`, even though
  most of them say "Parse Biosciences". (The audit found the same with `"*parsebio*"`.)
  `"*parse*" OR "*evercode*"` returns all 189.
- **8 that pass today**: their runs were first public between 6 Aug and 28 Sep 2026, and last
  updated before 8 Oct. ENA metadata does not show why the 8 Oct draw missed them. A re-draw under today's rules would add all 8. The proposal adds 7, because
  PRJNA1284247 names 10x and also used it.
- **Other vendor** (from the protocol text, judged by hand): 14 missing studies name 10x
  (PRJNA1255131 also names ScaleBio and Fluent; PRJNA1309117 also names ScaleBio). 3 only mention 10x (PRJNA1094867, PRJNA1169626 and PRJNA993024 cite
  10x dissociation or nuclei protocols). The others used that vendor's kit for libraries too.
  PRJNA1060422 credits a 10x Chromium kit to "Parse Biosciences", but it is a 10x study.
- **Non-WT Parse products** among the missing studies: PRJNA1455443 (TCR only), PRJNA1255131
  (Evercode TCR "WT + TCR v1" in a multi-vendor benchmark), PRJNA1244860 (WT + CRISPR guide
  library) and PRJNA1427207 (WT v3 + CRISPR Detect). No BCR or Gene Select appears anywhere in
  the 230 QUERY2 studies.

## Proposal (`proposed_rule.md`, `rules2.py`; not applied)

The revision scores the family, so the version is not needed. It accepts a protocol that names
the Parse WT kit with no version, as new kit `parse_wt`. It widens QUERY to `*parse*`. It refuses
other vendors, non-WT products and companion libraries, and adds a family shape rule
(R2 ≥ 58). It applies only to studies not already in the catalogue, and every new row is held out.

It adds **67 studies** (128 runs), listed in `proposed_added.tsv`. That is 67 held out if the
new `parse_wt` kit takes no design row (recommended), or 66 if its earliest study, PRJNA713904,
is made design to keep `drawn_sound`'s one-design-per-kit rule. All 67 name the WT kit and
none is TCR, BCR, CRISPR or Gene Select. 66 are certain WT. 1 (PRJNA1120008, scMPRA on WT Mini v2)
is uncertain. PRJNA1390207 is refused by name because its first paired runs are bulk QIAseq.

## Side finding: today's catalogue holds runs that are not Parse WT

The catalogue has no shape or metadata guard, and some held-out rows are other library types. ENA
run metadata shows:

| study | rows | what the runs are |
|---|---|---|
| PRJNA1208204 | 2 | "Parse CRISPR Detect single-cell barcode detection" (source OTHER) |
| PRJNA1304277 | 2 | "Mini_Sublibrary1/2_CRISPR" |
| PRJNA1456597 | 2 | "hU6-targeted gRNA ClonMapper barcode capture sublibrary" |
| PRJNA1463191 | 2 | "Illumina Parse CRISPR guide RNA sublibrary" |
| PRJNA1465771 | 2 | "KDM1A CUTRUN", GENOMIC (CUT&RUN, not Parse) |
| PRJNA1449947 | 2 | strategy OTHER, source GENOMIC ("neuron, E17.5_NFIB_WT_R1/2"): not RNA |
| PRJNA975472 | 2 | "Visium5xFAD_…" (10x Visium spatial, not Parse) |

That is 14 held-out rows in 7 studies. The 4 CRISPR guide libraries still carry Parse barcodes
on R2, so they are Parse but not WT gene expression. The other 3 studies (6 rows) are not Parse
libraries. Per VALIDATION.md section 5, the catalogue has already been run (a results file
exists; it was not opened). So this needs a recorded decision, not an edit by this task.

## Files

- `missing.tsv`: one row for each of the 111 missing audit studies. Columns: study, reason,
  detail (the sub-reason, vendor and product notes, and the proposal's outcome), protocol_window,
  n_runs, files_per_run, r1_len, r2_len (comma-separated, for up to 2 runs).
- `proposed_rule.md`: the revision, its counts and options.
- `proposed_added.tsv`: the 67 studies the revision adds, with kit, runs and evidence.
- `rules2.py`: the proposed QUERY2, classifier, shape and named refusal.
- `analyse.py`, `common.py`: the analysis. Run with `python3 -I -B analyse.py`.
- `cache/`: ENA rows (QUERY, QUERY2, the 28 studies outside QUERY), read lengths and run metadata.
  `stage1.txt` is the first-pass protocol dump.
