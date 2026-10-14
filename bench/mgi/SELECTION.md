# MGI catalogue: how it was extended on 9 Oct 2026

`results/catalogue.tsv` is append-only. The 54 rows drawn on 8 Oct are kept exactly as drawn,
because their runs had already been called. Rows from the 9 Oct draw are added only for studies
not already in it (24 studies, 48 held-out rows). A fresh draw with today's ENA would pick other
runs in two old studies (PRJNA1149752, PRJNA1502579), so `scripts/select.py` alone does not
reproduce the old rows.

One drawn study was left out by hand, for a stated reason: PRJNA1050786's protocol says
"Self-designed using MGI C4 droplet generation device". That is a home-made library on MGI's
device, not the kit, so its read structure is not the layout's.

The working notes of the 9 Oct draw follow.

# Revised MGI draw, 9 Oct 2026 (scratch only)

This is a draft. It does not change `bench/mgi/`. None of the drawn runs has been run through
`tchem`. This draw read only ENA metadata and the modal read lengths from `catalogue.read_length`
(first 64 KB, lengths only).

## Files

- `select.py`: a revised copy of `bench/mgi/scripts/select.py`. Its docstring records clauses H,
  B, G and T and the other-vendor check (applied to new clauses only), each with its reason. It
  also records that clause A was not adopted and that `mgi_c_series` has no design run.
  - It imports the repository's `catalogue.py` and `invariant.py` from an absolute path, because
    the copy lives outside the repository. If the copy goes into the repository, change this back
    to the relative path.
- `run_draw.py`: runs `select.draw()` once on one ENA snapshot, and records every read length it
  takes. It writes `catalogue.tsv`, `refused.tsv` and `snapshot.json`, so the refusals explain the
  same draw.
- `catalogue.tsv`: 106 rows, with columns `catalogue.OUT_FIELDS`.
- `refused.tsv`: 39 studies, each with its reason and detail. Shape refusals list the R1/R2 length
  of every run checked, up to 6.

## Design runs

`select.split_without_c_series_design` changes any design row whose kit is not `mgi_c4` or
`mgi_c_ht` to held out. It then checks that `mgi_c4` and `mgi_c_ht` have exactly one design run
each, and that `mgi_c_series` has none.

The design rows match today's: `mgi_c4` SRR25532325 (PRJNA1002278) and `mgi_c_ht` SRR25317687
(PRJNA995566).

## Counts

| Kit | Design | Held out | Studies |
|---|---|---|---|
| mgi_c4 | 1 | 23 | 12 |
| mgi_c_ht | 1 | 51 | 26 |
| mgi_c_series | 0 | 30 | 15 |
| Total | 2 | 104 | 53 |

The draw has 53 studies, against 28 in today's catalogue (+25):

- **4 studies pass the 8 Oct rules but are absent from the 8 Oct catalogue:** PRJNA1292106,
  PRJNA1501277, PRJNA1515196 and PRJNA1515281.
- **21 studies come from the revision:**
  - The 20 studies forecast in `../proposed_rule.md`.
  - PRJNA1468071 (`mgi_c_series`). Its first 4 runs are R1 20 / R2 100. Its runs 5 and 6 are 30 /
    100, and the shape check reaches them because the draw checks up to 6 runs; the forecast
    checked only 2. Its protocol names "DNBelab C Series kit (MGI, #940-001818-00)". Other
    protocols give that catalogue number to the HT RNA kit.
- **Other studies in `mgi_c_series` to note:**
  - PRJNA1050786 names no commercial set ("Self-designed using MGI C4 droplet generation
    device").
  - PRJNA1299265 and PRJNA1327849 name "C TaiM4 V3.0".

## Check against today's catalogue

All 28 studies are still drawn, each with the same kit. 52 of the 54 runs have the same kit and
split. There are 2 exceptions, both in studies that today's catalogue holds with 1 run:

- **PRJNA1149752** (`mgi_c4`):
  - Today's catalogue holds SRR30296180.
  - The redraw takes SRR30296166 and SRR30296168.
- **PRJNA1502579** (`mgi_c_ht`):
  - Today's catalogue holds SRR39847436.
  - The redraw takes SRR39847434 and SRR39847435.

The revision does not cause these. Both studies classify under the 8 Oct rules, and the shared
draw rules are unchanged. Today, the earlier runs have R1 30 / R2 100, so the sorted shape check
takes them first. SRR30296180 is not among the first 6 sorted runs, so the current code cannot
draw it from today's ENA. ENA metadata does not show why the 8 Oct draw took different runs.

## Refusals (39 studies)

| Reason | Studies |
|---|---|
| shape refused (up to 6 runs checked) | 15 |
| no paired run | 12 |
| classify None: ATAC (unchanged rule) | 8 |
| classify None: names another vendor's kit (10x): PRJNA1258882, PRJNA902352 | 2 |
| classify None: names both, and C4 names a kit: PRJNA1266728, PRJNA1453268 | 2 |
