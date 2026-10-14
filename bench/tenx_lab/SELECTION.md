# 10x lab-labelled catalogue: the draw of 9 Oct 2026

Drawn by a metadata-only agent. The runs come from the 10x pool, and each label comes from the
lab's own record (ENA experiment fields and the linked GEO sample), mapped by the committed
`bench/head_to_head/scripts/lab_labels/label10x.py`. Every row is held out, and its referee is "the
lab's own record". `results/statements_drawn.tsv.gz` holds every mapped sentence for the drawn runs,
with its URL. The scripts are kept as run: they read the agent's scratch downloads, which are not
committed.

## 10x held-out catalogue, labelled from the labs' own records (built 9 Oct 2026)

Sampling frame: the `srx` column of `arc_bench/pool.tsv` and `pool_recent.tsv` (8,000 unique SRX). No other
pool column was read. Labels come from each run's ENA fields and linked GEO sample, mapped by the committed
`bench/head_to_head/scripts/lab_labels/label10x.py` (imported unchanged; self-test 49/49). No sequence content
was read: only modal read lengths, through `catalogue.read_length` (first 64 KB). `tchem` was not run.

## Files

| File | What |
|---|---|
| `catalogue.tsv` | 624 rows in `catalogue.OUT_FIELDS`, all `held_out`, kit = label key, evidence = first sentence stating the run's exact label (at most 120 characters) |
| `labels.tsv` | srx, run, study, label, n_statements: one row per run (6,268 runs) |
| `statements.tsv` | every candidate sentence, quoted, with its URL (columns of `lab_label_statements.tsv`; no run column, so an SRX with several runs repeats its rows) |
| `excluded.tsv` | srx, reason: every pool SRX not drawn (7,531); 7,531 + 469 drawn SRX = 8,000 |
| `examined_wild.tsv` | run, study, what: 30 runs examined in the wild sweep or `validate.sh` before the draw; `prep.py` excludes their studies |
| `fetch_failures.tsv` | empty (header only): no ENA chunk, GEO sample or read-length fetch failed after one retry |
| `scripts/` | prep, fetch_ena, make_input, fetch_geo, build, draw, finalize |
| `work/` | build.diff (build.py copy vs the committed one), input_runs.tsv, run_labels.json, shape_log.tsv (every length check), summary.json |
| `raw/` | ENA chunks (188) and GEO brief texts (2,415) |

## Procedure

1. `prep.py`: pool SRX; drop SRX in `lab_label_runs.tsv`, `studies.tsv`, `arc/results/selection.tsv` (550 pool SRX).
   Excluded studies: those three files' studies (arc studies through ENA; `NRX0000030` has no ENA record), every
   `bench/*/results/catalogue.tsv` study with a `design` or `examined` row (files with a `split` column), and
   `exec2/examined_wild.tsv`. All study IDs are PRJ accessions, like ENA's `study_accession`.
2. `fetch_ena.py`: ENA portal `result=read_run`, 40 SRX per query, fields of the committed script plus `fastq_ftp`,
   `first_public`. `make_input.py`: every run of a remaining SRX whose study is not excluded (6,268 runs).
3. `fetch_geo.py`: GEO brief text for every GSM in those runs' sample_alias / experiment_alias / experiment_title,
   0.34 s throttle, one retry, `^SAMPLE = GSM` check (as committed).
4. `build.py`: the committed build.py, changed only as in `work/build.diff` (input keyed per run, run column,
   no conflict reasons, no 465 assert, records each run's first exact-stating sentence).
5. `draw.py`: keep the six labels; exactly two FASTQ files; per study, runs in accession order, at most
   `SHAPE_CHECKS_PER_STUDY` (6) checked, shape = one file's modal length >= 26 and the other's >= 50, first
   `RUNS_PER_STUDY` (2) that pass. `catalogue.drawn_sound(rows, design=False)` and `catalogue.write`. If the
   300-character quote was over 120 characters, the evidence is the first 120-character window (start, then 40
   before each decisive match) whose `label10x.label` is still the run's label (234 moved, 0 lost).

## Counts

Before the draw (labels.tsv, runs / studies): 3' v3/v3.1 2,594 / 427; not stated 1,363 / 242; 5' v1/v2 715 / 135;
3' v2 579 / 91; end_only (3') 449 / 73; end_only (5') 294 / 49; Multiome GEX 159 / 35; conflict 82 / 19;
5' v3 (GEM-X) 12 / 4; 3' v4 (GEM-X) 12 / 1; Flex v1 8 / 3; unsupported (3' v1) 1 / 1; ambiguous 0 (the build
never gives it at run level).

After the draw (catalogue.tsv, runs / studies): tenx_3p_v3 360 / 219; tenx_5p_v1v2 135 / 79; tenx_3p_v2 94 / 55;
tenx_multiome 29 / 18; tenx_5p_v3 4 / 2; tenx_3p_v4 2 / 1. Total 624 runs, 374 studies, no study in two kits.

4,071 runs had a kept label; 1,945 of them had exactly two files; 636 runs were length-checked, 12 failed.

## Caveats

- The label is the experiment's record, not the library's. A 5' study that also made V(D)J libraries, or a
  Multiome study whose ATAC run was exported as two 50+ bp files, can give a drawn run that is not the gene
  expression library; the length rule does not catch either. These would show as disputes, not as rule errors.
- The 2 tenx_3p_v4 rows are from one study; tenx_5p_v3 has 2 studies.
- The GEO fetch hit a disk-quota error on `/tmp` once; 25 empty or block-truncated GSM files were deleted and
  refetched (final fetch: 2,415 GSM, 0 failed). Integrity pass: all 2,426 GSM files in raw/geo (11 extra are
  from the first, wider run and are not used) end with the closing `!Sample_data_row_count` line; all were
  refetched (`scripts/check_geo.py`, into raw/geo_check/) and compared with `cmp`: 0 differ. (check_geo.py's own
  in-Python comparison reports all as differing because it reads the file in text mode, which folds GEO's
  `\r\n`; the byte comparison is the valid one.) Every ENA chunk and output TSV parses with a constant field
  count and ends in a newline.
- Exclusion reasons are per SRX, first reason only. Study precedence is head_to_head, then arc, then catalogue
  design/examined, then examined_wild, so a study in both head_to_head and a catalogue counts under head_to_head.
- Not excludable by study: `NRX0000030` (arc selection, no ENA record) and 3 rows of `studies.tsv` with no
  study; their SRX are excluded.
- 64 pool SRX have no ENA read_run record (the queries succeeded); they are in excluded.tsv.
- No pool SRX fell in a design/examined catalogue study or `examined_wild.tsv` study after the SRX and
  used-study exclusions (those reasons have 0 rows).

## Re-running (10 Oct 2026)

The scripts first ran with absolute paths on the machine that drew the catalogue. They now resolve
the repository from their own location and read the Arc pool from `$TUSSL_CHEMISTRY_DATA/arc_bench/`.
`examined_wild.tsv`, an input `prep.py` read from outside the repository, is now committed in
`results/`. Run today, `prep.py` gives the same pool, used SRX and query set as at the draw. Its
exclusion set has 3 more studies: SeekGene's design studies (PRJNA913187, PRJNA947149,
PRJNA1064036), whose catalogue was committed with this one. None is in the catalogue or
`excluded.tsv`.
