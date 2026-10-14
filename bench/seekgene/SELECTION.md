# SeekGene catalogue: the draw of 9 Oct 2026

Drawn by a metadata-only agent by the rules in `scripts/select.py`, revised once after the first
draft and before any run (the docstring records each revision). Its working notes follow.

## SeekGene SeekOne held-out catalogue (draft, 9 Oct 2026)

This catalogue was drawn by `select.py` with the rules in its docstring. They were fixed before the draw.

- **Inputs.** The draw used ENA metadata and read lengths only. `catalogue.read_length` takes the modal line length from the first 64 KB of each FASTQ; the bases are not used.
- **Not done.** No FASTQ content was inspected, and tchem was not run.
- **How to run it.** Outside the repository:
  `TUSSL_BENCH=<repo>/bench python3 -I -B select.py catalogue.tsv refused.tsv`
  Inside the repository it would sit at `bench/seekgene/scripts/` and find `heldout/scripts/catalogue.py` by relative path. The run took about 10 minutes, almost all of it fetching read lengths.

## Files

| File | Content |
|---|---|
| `catalogue.tsv` | `catalogue.OUT_FIELDS`. 110 rows from 57 studies. |
| `examined.tsv` | The three design runs read before the draw. |
| `refused.tsv` | One row per query study left out, with the reason. 16 studies. |

## Counts (after the revisions)

| Kit | Claimed by tool | Design | Examined | Held out (runs / studies) | Total (runs / studies) |
|---|---|---|---|---|---|
| seekgene_dd3 | yes (DDV2) | 1 (SRR27533148) | 1 (SRR27533149) | 57 / 29 | 59 / 30 |
| seekgene_dd5 | yes (DD5V1) | 1 (SRR23931570) | 1 (SRR23931571) | 6 / 3 | 8 / 4 |
| seekgene_mm | no | 1 (SRR22784277) | 1 (SRR22784278) | 26 / 13 | 28 / 14 |
| seekgene_ddq | no | 0 | 0 | 7 / 4 | 7 / 4 |
| seekgene_dd (end not stated) | no | 0 | 0 | 8 / 5 | 8 / 5 |

The query (`*seekgene*` OR `*seekone*`) returned 73 studies: 57 drawn and 16 refused.

## Refusals

| Reason | Studies | Accessions |
|---|---|---|
| Names another vendor's kit | 8 | PRJNA1523965 (Singleron), PRJNA1063271, PRJNA1222787, PRJNA1036114, PRJNA1268354, PRJNA1279856, PRJNA1253968, PRJNA1167743 (10x) |
| Classified, but no paired run | 5 | PRJNA1511341, PRJNA1412541, PRJNA1529080 (dd3); PRJEB106495 (dd5); PRJNA1164785 (ddq) |
| Mixed: names both an MM kit and a DD kit (revision 1) | 1 | PRJNA1224797 |
| SeekSpace, Multiome or ATAC | 1 | PRJNA1458007 |
| SeekGene for service only | 1 | PRJNA1322363 |

## Revisions (9 Oct 2026, approved, made before any catalogue run was read)

These revisions are recorded in the `select.py` docstring. Draft 1 is kept as `catalogue_draft1.tsv`, `refused_draft1.tsv` and `select_draft1.py`.

1. **Mixed kits.** A protocol that names both an MM kit and a DD kit is refused. This refuses PRJNA1224797, which uses the MM kit for scRNA and the DD 3' kit for snRNA.
2. **Shape.** A run is also accepted when R1 >= 29 and R2 >= 40. This admits PRJEB93881 (R1 29, R2 43).
3. **Evidence.** The window is cut from the protocol itself, around the deciding match, instead of from the joined mention windows.

## Flag kept: PRJDB36774

PRJDB36774 stays `seekgene_dd` (end not stated). Its text says it is multiplexed with HEATag labels "on SeekOne DD", so it may be a SimpliPlex or DDVS-type library, with a sample barcode on R1. Check it before any score for `seekgene_dd` is read. The group is not claimed by the tool.

## Checks done

- Within each query study, every run's protocol classifies to the same label, so the draw's first-run-only classification loses nothing.
- Each design run was drawn as the first run of its study. No design row had to be added.
- Every query study appears in either `catalogue.tsv` or `refused.tsv`.
