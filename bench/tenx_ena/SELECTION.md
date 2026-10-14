# 10x catalogue from ENA protocols: the draw of 9 Oct 2026

Drawn by a metadata-only agent by the rules in `scripts/select.py`, from every ENA run whose library
protocol names 10x, Chromium, GEM-X or Multiome. Labels come from the protocol text, mapped by the
committed `bench/head_to_head/scripts/lab_labels/label10x.py`. 1,500 studies were sampled (seed
20261009) from a frame of 6,251. A second batch of 2,100 studies (10 Oct 2026) continued the same
order; `results/not_sampled.tsv` lists the 2,618 studies still left. The batch 1 working notes follow, then
batch 2's.

## 10x single-cell RNA catalogue from ENA protocol text: the draw of 9 Oct 2026

Drawn by a metadata-only agent from all of ENA. Each label comes from the run's own `library_construction_protocol`, mapped by the committed `bench/head_to_head/scripts/lab_labels/label10x.py` (imported unchanged at commit c157518 of this repository; self-test 49/49). No sequence content was read: only modal read lengths, through `catalogue.read_length` (first 64 KB). `tchem` was not run. Rules: `select.py` docstring.

## Counts (catalogue.tsv, every row `held_out`)

| kit | runs | studies |
|---|---:|---:|
| tenx_3p_v2 | 486 | 252 |
| tenx_3p_v3 | 1820 | 931 |
| tenx_3p_v4 | 68 | 36 |
| tenx_5p_v1v2 | 310 | 158 |
| tenx_5p_v3 | 21 | 12 |
| tenx_multiome | 221 | 111 |
| total | 2926 | 1500 |

## Sampling

The frame (studies passing every rule but the shape rule) held 6,251 studies, over the cap of 1,500: 10x 3' v3/v3.1 3,876, 10x 3' v2 1,015, 10x 5' v1/v2 713, 10x Multiome GEX 435, 10x 3' v4 (GEM-X) 153, 10x 5' v3 (GEM-X) 59. The frame's accessions were sorted, shuffled with `random.Random(20261009)`, and walked in that order, each study shape-checked, until 1,500 had a kept run. 1,518 studies were walked (18 failed the shape rule); the 4,733 after the stop were not checked and are listed in `work/not_sampled.tsv` (not refused). 3,045 runs were length-checked, 119 failed. So per-kit study counts follow the frame's mix, not a per-kit quota.

## Refused studies (refused.tsv, one row per study)

| reason | studies |
|---|---:|
| label: not stated | 4,635 |
| no paired run (exactly two FASTQ files) with the study's label | 4,417 |
| used (already in a bench source; excluded, not a rule refusal) | 1,968 |
| label: end_only (3') | 1,389 |
| another assay: ATAC without Multiome | 891 |
| label: conflict | 625 |
| label: end_only (5') | 613 |
| another assay: Visium or spatial | 435 |
| another assay: V(D)J only, no gene expression | 192 |
| another assay: Flex or fixed RNA profiling | 124 |
| runs carry more than one label | 88 |
| another vendor named | 80 |
| label: unsupported (3' v1) | 53 |
| shape: no paired run passed | 18 |
| total | 15,528 |

## Queries (ENA portal, result=read_run, 9 Oct 2026)

| query | runs |
|---|---:|
| `library_construction_protocol="*10x*"` | 737,473 |
| `library_construction_protocol="*chromium*"` | 548,822 |
| `library_construction_protocol="*gemx*"` | 120 |
| `library_construction_protocol="*multiome*"` | 32,287 |
| `library_construction_protocol="*single?cell?3*"` | 37 |
| `library_construction_protocol="*single?cell?5*"` | 390 |
| `library_construction_protocol="*singlecell3*"` | 1,799 |
| `library_construction_protocol="*singlecell5*"` | 126 |
| `union of the eight terms` | 786,534 |
| `library_construction_protocol="*gem-x*"` (count probe, not used) | 0 |
| `library_construction_protocol="*GEM-X*"` (count probe, not used) | 0 |
| `library_construction_protocol="*gem x*"` (count probe, not used) | 0 |
| `library_construction_protocol="*gem?x*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single cell 3'*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single cell 3’*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single cell 5'*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single cell 5’*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single-cell 3*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single-cell 5*"` (count probe, not used) | 0 |
| `library_construction_protocol="*single cell 3*"` (count probe, not used) | 0 |
| `library_construction_protocol="*3 prime*"` (count probe, not used) | 0 |
| `library_construction_protocol="*immune profiling*"` (count probe, not used) | 0 |
| `library_construction_protocol="*v(d)j*"` (count probe, not used) | 0 |
| `library_construction_protocol="*multi?ome*"` (count probe, not used) | 2 |

ENA wildcards match nothing across a space, hyphen or apostrophe, so the phrase queries return nothing; their runs are reached through `*10x*` and `*chromium*`. The union was fetched as 12 queries by `first_public` year (`work/queries.json`); each partition's row count equals its count, and the partitions add up to the union.

## Five random rows (`random.Random(9).sample(rows, 5)`)

| kit | study | run | evidence |
|---|---|---|---|
| tenx_3p_v3 | PRJNA1205713 | SRR31874125 | acturer's protocol CG000315 Rev C for the Chromium Next GEM Single Cell 3' Kit v3.1 (10X Genomics, Pleasanton, CA, USA). |
| tenx_5p_v1v2 | PRJNA1113333 | SRR29079234 | wed by encapsulation in a lipid droplet (Single Cell 5′kit V2, 10X Genomics) followed by cDNA and library generation acc |
| tenx_3p_v3 | PRJNA853194 | SRR19859694 |  cDNAs with Chromium Next GEM Single Cell 3' Reagent Kits v3.1. Briefly, blood cell suspensions were loaded on a 10X Gen |
| tenx_3p_v3 | PRJNA966751 | SRR24413215 | nual of the Chromium Next GEM Single-cell 3' Reagent Kits v3.1. The cell suspensions were loaded on the Single-cell G Ch |
| tenx_3p_v3 | PRJNA757967 | SRR15618885 | using the Chromium Single Cell 3' Library and Gel Bead Kit v3 on the Chromium Controller (10X Genomics) following manufa |

## Checks

- Recount from the files (`work/recount.py`): 2,926 rows, all `held_out`, 1,500 studies, at most 2 runs per study, no study
  in two kits, no run twice; refused.tsv has 15,528 studies, none in the catalogue or in `work/not_sampled.tsv`;
  15,510 refused before the frame + 6,251 in the frame = 21,761 studies returned by the queries.
- Exclusion: no catalogue study, run, experiment or secondary study accession is in the used set (3,156 PRJ studies mapped from
  the sources in the docstring by `work/tokens.py` and `work/map_used.py`). 45 accessions could not be mapped through the ENA portal: 17 SRX (from tenx_lab's excluded.tsv; one is also in
  head_to_head) and 28 GSE series (long_read and flex rows, each beside its PRJ accession in the same row, so the study is
  excluded anyway). After the draw the 17 SRX were looked up again through NCBI E-utilities runinfo and the ENA browser XML
  (`work/unmapped_srx.json`): 11 map to PRJNA1103675 and PRJNA1192233, neither in the catalogue; 6 (SRX14779154,
  SRX14779158, SRX19902442, SRX28844507, SRX28844510, SRX28844511) name no study in either source, so their studies are
  unknown and could not be excluded. The 3 ENA 404s are in fetch_failures.tsv.
- Evidence: every window is a verbatim slice of its run's protocol (all 2,926 also without the tab/CR/LF replacement), and
  re-fetched from ENA for 10 random rows (seeds 20261009 and 9, `work/verify.py`): 10/10 verbatim, study matches.
  `label10x.label(window)` gives the row's label for 2,924 of 2,926 rows. The 2 others (SRR32034434, SRR32034435, ENCODE
  umbrella PRJNA66167) are labelled Multiome from a protocol URL, `.../Chromium_Multiome_ATAC_GEX.pdf`; the full
  sentence maps to Multiome GEX, a 120-character window does not. Their modal lengths are 28 and 90, a gene-expression shape.
- fetch_failures.tsv: 3 post-draw SRX lookups (above), and 1 read-length fetch (SRR31815262, R1) failed on both attempts; that run (PRJNA1202085) counted as failing the shape rule, and two other runs of its study were drawn. No ENA
  portal query or accession lookup made for the draw failed.

## Batch 2: 10 Oct 2026 (append-only)

Drawn by the same kind of agent and under the same rules, metadata only. The rules were fixed before its walk and are in `select.py` (docstring section "Batch 2"). `select.py <dir>` still reproduces batch 1 exactly: the same frame, the same 15,510 refusals, the same order and the same 2,926 rows.

- Order: batch 1's own (the batch-1 frame sorted, then shuffled with `random.Random(20261009)`), not reshuffled.
  - Batch 1 walked positions 0 to 1,517.
  - Batch 2 walked positions 1,518 to 3,632 until 2,100 more studies had a kept run: 2,100 were drawn and 15 refused for shape (`results/refused_batch2.tsv`, with position).
  - `results/not_sampled.tsv` now holds the remaining positions 3,633 to 6,250 (2,618 studies), in order.
- Size: 2,100 studies. This was raised from 1,100 before any run, once disk was freed. A first pass was interrupted at about 900 drawn studies and kept nothing, so the walk was re-run whole from position 1,518.
- Frame re-derived from today's ENA (`results/queries_batch2.json`): the same 786,534 runs with identical fields, only in another row order. No walked study changed and none is new to the frame.
- Exclusion: today's used set adds batch 1 and every other bench catalogue at 3759832, plus PRJNA1103675 and PRJNA1192233 (batch 1's post-draw SRX check), for 4,658 studies in all. No drawn study or run is in it, and none is in batch 1.

| kit | runs | studies |
|---|---:|---:|
| tenx_3p_v2 | 683 | 352 |
| tenx_3p_v3 | 2485 | 1278 |
| tenx_3p_v4 | 117 | 60 |
| tenx_5p_v1v2 | 488 | 248 |
| tenx_5p_v3 | 27 | 14 |
| tenx_multiome | 293 | 148 |
| total | 4093 | 2100 |

Checks, recounted from the files outside the agent:
- 4,093 rows, all `held_out`, at most 2 runs per study, no run twice, and no study in two kits.
- No study or run overlaps any bench file.
- The 2,115 walked studies are exactly the first 2,115 of batch 1's not_sampled list, and are the drawn plus refused studies.
- Evidence windows: 5 random rows (seed 4242) re-fetched from ENA, 5/5 verbatim with the study matching; the agent checked 10 more, also 10/10.
- `label10x.label(window)` gives the row's label for 4,091 of 4,093 rows. The 2 others (SRR32481539 and SRR32481540, PRJNA1227983) state "v3.0 chemistry" more than 120 characters from their "3'".
- One read-length fetch failed on both attempts (SRR37573406, R2; `results/fetch_failures.tsv`). That run counted as failing the shape rule.
- Shape checks are appended to `results/shape_log.tsv`.

## Caveats

- The label is the experiment's protocol, not the library's. A 5' study with V(D)J libraries, or a Multiome study whose
  ATAC runs were exported as two long files, can give a drawn run that is not the gene-expression library when the same
  protocol text covers both; the shape rule does not catch either. PRJNA66167 (ENCODE, many assays) is one such risk.
  For scale: 167 of 221 tenx_multiome rows (90 studies) have both modal lengths >= 40, so neither file is a 24-28 base
  barcode read; the same holds for 846 of 1,820 tenx_3p_v3 rows (long R1 sequencing is common), so this is not by itself
  a sign of ATAC.
- Stricter than build.py: a sentence that `label10x` calls ambiguous makes the run "conflict".
- "not stated" runs beside a kept label do not refuse a study (rule 4); 23 frame studies have such runs. They are never drawn.
- Refused for "another assay" is over the whole protocol text: a 3' study that also mentions spatial or ATAC data is
  refused even if its RNA runs are clean.

## Files

| File | What |
|---|---|
| `select.py` | the rules (docstring, 9 Oct 2026) and the draw; imports catalogue.py and label10x.py through `TUSSL_BENCH` |
| `catalogue.tsv` | the draw, `catalogue.OUT_FIELDS` |
| `refused.tsv` | study, reason, run labels, first run, first 200 characters of its protocol |
| `fetch_failures.tsv` | every fetch that failed after its retry |
| `work/` | fetch.py, queries.json, tokens.py, used_tokens.json, map_used.py, used_studies.json, shape_log.tsv, not_sampled.tsv, summary.json, recount.py, verify.py, readme.py |
| `raw/` | the ENA query results, one gzipped TSV per year |
