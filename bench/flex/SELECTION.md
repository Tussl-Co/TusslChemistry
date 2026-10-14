# 10x Flex held-out catalogue: draft 3 (9 Oct 2026)

This is a draft for `bench/flex/scripts/select.py` and `bench/flex/results/catalogue.tsv`.
Nothing in the repository was changed. Draft 3 keeps draft 2's query, its stricter Flex statement and
its merge of the old catalogue. The one change is the shape rule. It is now
`min(r1, r2) >= 28 and max(r1, r2) >= 45`, which replaces `r1 == 28 and r2 >= 45`. A Flex
barcode read sequenced long (150 or 151, with barcode and UMI at its start) is a genuine
Flex run, and either file may carry the barcode read. Earlier draws are kept as
`catalogue_v1.tsv`/`refused_v1.tsv` (first rules) and `catalogue_v2.tsv`/`refused_v2.tsv`
(revised rules, old shape).

## Files

- `select.py`: the draft. The relative `sys.path.insert` line works after you copy the file to
  `bench/flex/scripts/`. The absolute `sys.path.append` line is for this scratch copy only;
  remove it when you copy the file.
  `main(out, old=../results/catalogue.tsv)` reads the old catalogue before it writes. If the
  old file already has an `origin` column, the script reuses that file's `old_catalogue` rows,
  so a redraw over its own output keeps them.
- `catalogue.tsv`: catalogue.py's `OUT_FIELDS` plus `origin` (`draw` or `old_catalogue`).
  Every row is `held_out`.
- `refused.tsv`: one row for each study that the query returned and the draw did not keep. The
  columns are study, run, reason, detail and window. For a shape refusal, detail gives
  `run:R1xR2` for each run that was checked. Old-catalogue studies that the draw refused are
  listed here too, and they are also in the catalogue with origin `old_catalogue`.
- `audit.py`: fetches ENA once (one retry if the fetch fails; no retry was needed), runs
  `catalogue.draw(QUERY, classify, shape, design=False)` and `merge`, and writes both TSVs.
- `old_studies.txt`, `notreturned.txt`, `query_probe.txt`: notes from draft 1.

## Commands

    python -I select.py catalogue.tsv ../results/catalogue.tsv
    python -I audit.py .     # produced the files here

## Counts

ENA returned 327,349 runs from 3,546 studies, the same as for draft 2.

| kit | origin | runs | studies |
|---|---|---|---|
| 10x_flex_v1 | draw | 242 | 128 |
| 10x_flex_v1 | old_catalogue | 5 | 3 |
| 10x_gemx_flex_v2 | draw | 37 | 20 |
| 10x_gemx_flex_v2 | old_catalogue | 0 | 0 |
| total | | 284 | 151 |

Refused studies (3,398):

| reason | studies |
|---|---|
| not paired (dropped by catalogue.draw's pairing rule, before the classifier) | 1,516 |
| no 10x context | 1,302 |
| no flex statement | 404 |
| visium | 156 |
| other vendor | 18 |
| shape | 2 |

No read length of 0 was recorded, so no file was unreadable.

## Compared with draft 2

- 46 of draft 2's 48 shape-refused studies are now drawn. The 2 still refused are
  PRJNA1375788 (27x89, so one read is shorter than 28) and PRJNA743932 (38x38 and 36x40, so
  neither read reaches 45).
- The other refusal counts did not change.
- The 9 QuantStudio/Qubit/E8-medium studies are still refused as "no flex statement".

## The old catalogue (26 studies)

- 23 are in the draw, each with the same kit as before. The 9 studies that were merged in
  draft 2 because of shape are now drawn: PRJNA1182275, PRJNA1208002, PRJNA1225771,
  PRJNA1247492, PRJNA1285404, PRJNA1285691, PRJNA1329031, PRJNA1365166 and PRJNA1494591.
- 3 are merged as `old_catalogue` rows (5 runs):
  - PRJNA1439199: "no flex statement". It says "Fixed RNA Sample Preparation Kit", not
    "Fixed RNA Profiling";
  - PRJNA1122589: not in the query result, because its ENA protocol is empty;
  - PRJNA1330221: not in the query result.

The coordinator corrected an earlier statement about the merged rows. They did not all have
28-base R1s. They were checked against GEO records, and the tool called them in earlier
coverage runs, at the lengths that the old catalogue records. For the 3 studies merged now,
those lengths are: PRJNA1439199 R1 28, R2 90; PRJNA1122589 R1 28, R2 90; PRJNA1330221 R1 150,
R2 150. This draft did not check the coverage runs (`results*.tsv` was not opened).
