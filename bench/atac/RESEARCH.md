> Research notes written before this kit's catalogue was drawn (8 Oct 2026), kept as the record of
> sources and of which runs' reads were examined (those runs' studies are never held out: see
> results/examined.tsv). Scripts these notes name were scratch work and are not kept; the
> catalogue's rules are the ones in scripts/select.py, which say where they differ from these notes.

# 10x single-cell ATAC (scATAC v1 / v1.1 / v2, Multiome ATAC): held-out catalogue research

Researched 2026-10-08. Nothing in the repo was changed. Working files are in `scratchpad/atac_work/`. Unsourced claims are marked **to verify**.

## 1. Read structure

| Kit | Read 1N | i7 | i5 = barcode read | Read 2N |
|---|---|---|---|---|
| Single Cell ATAC v1, CG000168 Rev A [G1] | 50 | 8 | **16** | 50 |
| Next GEM Single Cell ATAC v2, CG000496 Rev A (PN-1000390/1000406) [G2] | 50 | 8 | **16** | 50 |
| Multiome ATAC, CG000338 Rev F [G3] | 50 | 8 | **24** ("16 bp 10x Barcode and 8 bp Spacer") | 49 |

- **File names.** Cell Ranger ATAC names the files R1 = Read 1N, **R2 = i5 barcode**, R3 = Read 2N, I1 = sample index. Source: `BCL2FASTQ_SEQNAMES` in cellranger-atac 1.2, `lib/python/constants.py` [S1].
- **16-bp Multiome reads.** Multiome on a NextSeq 500/550 needs "8 dark cycles and 16 nt cycles on i5" [G3], which gives a 16-bp read. On ENA, 6 of 92 Multiome candidate runs (4 studies) had a 16-bp barcode read; only 1 of those studies used a NextSeq 500. So a 16-bp read does not imply scATAC.
- **Versions to verify.**
  - v1.1: no v1.1 guide was found. Teichlab's "Next GEM v1.1" example uses cratac [T1].
  - Multiome v2: Cell Ranger 10.1.0 `chemistry_defs.json` describes ARC-v1 as "Multiome ATAC + Gene Expression v1 or v2".

**Orientation.** I read 64 KB of the barcode read from each of 19 runs and scored its last 16 bp against every list in both orientations. No run matched the wrong list.
- scATAC (11 runs) matched cratac:
  - forward in 3 runs, all NovaSeq 6000;
  - reverse complement in 8 runs: NovaSeq 6000, X Plus, NextSeq 500/2000. Match rates were 0.84–0.96.
- Multiome (8 runs: NovaSeq 6000, X, X Plus) matched the ARC ATAC list reverse complement in every run (0.92–0.97). Positions 1–8 were always `CAGACGCG`.
- Teichlab's explanation [T1, T2]:
  - forward-strand i5 machines (MiSeq, HiSeq 2000/2500, NovaSeq 6000 v1.0) read the barcode as listed;
  - reverse-complement machines (NextSeq, HiSeq X/3000/4000, NovaSeq 6000 v1.5) read its reverse complement.
- CG000338 allows pooling Multiome ATAC with other ATAC libraries "only when using forward strand Illumina workflow" [G3].
- Cell Ranger ATAC 1.2 hard-codes `barcode_reverse_complement: False` [S1]. Auto-detection in later versions (first or last 16 bp, either orientation) [W1] is **to verify**, because the 10x site returned 403/429.

## 2. Barcode lists

All entries in this table were verified by set comparison.

| List | Count | Where |
|---|---|---|
| `737K-cratac-v1` (scATAC v1/v1.1/v2) | 737,280 | cellranger-atac repo [S2]; md5 `ab645450…`, identical to Teichlab [T3]. Not in the Cell Ranger 10.1.0 tarball. |
| ARC ATAC `737K-arc-v1` (`lib/python/atac/barcodes` in cellranger-arc) | 736,320 | Teichlab `atac_737K-arc-v1` [T4]. **The local `10x_arc_v1_gex_label_only.txt` is this ATAC list.** It shares 688 entries with cratac and none with the GEX list in either orientation. |
| ARC GEX `737K-arc-v1` | 736,320 | Cell Ranger 10.1.0 `lib/python/cellranger/barcodes/`. Same set as the local `10x_arc_v1_gex_permit.txt`. |

- The two ARC lists pair by line, one gel bead per line [T2]. Not independently verified.
- simpleaf maps `10x-atac-v1` and `10x-atac-v2` to the same file `10x_atac_v1_v11_v2.txt` (cratac), and maps `10x-arc-atac-v1` separately [C1].
- Terms: cellranger-atac is MIT with a clause limiting use to "a 10x Genomics Product … or data generated using" one [S3]. Whether tussl may redistribute the lists is **to verify**.

## 3. How ENA presents these runs

Portal searches use `library_strategy="ATAC-seq"`. Protocol matching is case-insensitive.

| Protocol filter | Runs | Studies | Runs by number of FASTQ files (0 / 1 / 2 / 3 / 4) |
|---|---|---|---|
| `*chromium*` | 19,191 | 983 | 2,181 / 556 / 14,675 / 351 / 1,428 |
| broad (§4) | 42,262 | 1,451 | 3,043 / 995 / 36,096 / 546 / 1,582 |

Multiome runs filed under other strategies (21,478 runs) have no 4-file runs, so they add nothing.

Read lengths, first 64 KB of each file:
- **2-file runs (85% of all runs)** lacked the barcode read in all 12 sampled; they were genomic pairs (34–151 bp). One exception is known: the repo's `bench/arc/results/results_wild.tsv` lists "10x-ATAC" 2-file runs of 89/28, which are GEX-shaped and probably misfiled.
- **3-file runs `ACC.fastq.gz` + `_1` + `_2`** had no barcode in 7 of 7 sampled.
- **1-file runs** were the barcode read alone in 7 of 9 sampled.
- **4-file runs and 3-file runs named by index (`_I1/_R1/_R2` or `_1/_2/_3`)** held the barcode read in **114 of 114 candidate runs**, with a 16-, 24- or 26-bp file in any position.
- Every barcode-bearing run was first public between 2025-08-08 and 2026-09-03.

## 4. Proposal

**Query.** `library_strategy="ATAC-seq"`, with `library_construction_protocol` matching any of: `*chromium*`, `*10x*`, `*multiome*`, `*scATAC*`, `*snATAC*`, `*single cell ATAC*`, `*single-cell ATAC*`, `*single nucle*`, `*cellranger*`, `*cell ranger*`. Add `submitted_ftp,instrument_model` to FIELDS.

**Classifier** (full regexes in `atac_work/propose.py`):
1. Require a 10x word (`10x genomics|chromium|cell ?ranger|\b10x\b`, not "10x buffer/TD/PBS") or a Multiome word.
2. Refuse a protocol that names another vendor's kit: `rhapsody|evercode|parse bio|scale ?bio|surecell|ddseq|hydrop|sci-?atac|dscatac|dnbelab|\bmgi\b|pipseq|snare-seq|share-seq|paired-seq|s3-atac|txci|spear-atac|issaac`.
   - Do not refuse on bare `BD` or `Bio-Rad`. These words name sorters and reagents in genuine 10x protocols; they caused 6 of 9 false refusals in a first pass.
3. If the protocol matches `multiome|multi-ome|atac\s*\+\s*(gene expression|gex)|cell ?ranger[- ]arc|1000283|1000285`, label the run **10x Multiome ATAC**.
4. Else if it matches `single[- ](cell|nucle(i|us)) atac|\bs[cn]-?atac|chromium.{0,60}atac|cell ?ranger[- ]atac`, label it **10x scATAC**.
   - v1, v1.1 and v2 share one list and one read structure, so they get one label.
   - Record the version in `evidence` only (v2: `\bv2\b|1000390|1000406`).

**Shape rule** (3 or 4 files):
- exactly one file has a modal length of 16, 24 or 26 bp (the barcode);
- exactly two files are 30 bp or longer (the genomic reads, in file-name order);
- an 8- or 10-bp file (the sample index) is optional;
- 3-file runs that include a bare `ACC.fastq.gz` are dropped;
- a scATAC label with a 24/26-bp barcode is left out. This drops 2 studies: PRJNA1305282, which cites CG000505, the Multiome nuclei protocol, and PRJNA1285399, which says "v1.1" (**to verify**). A Multiome label with a 16-bp barcode is allowed.

**Changes to `draw()`:**
- add a `files=(2,)` argument; ATAC passes `(3, 4)`;
- `shape` receives file names as well as lengths;
- rows carry `urls` (all files), `barcode_url`, `r1_url` and `r2_url`, and `drawn_sound` checks that all URLs differ;
- add an **`exclude` set of run accessions** (the "Runs examined" list below). Without it, `draw()` takes the first runs by accession and can pick examined runs again.

**Pool**, at most 2 runs per study, with the shape rule applied to every candidate run. One run per kit becomes the design run (Multiome SRR24653523, PRJNA974167; scATAC SRR36097788, PRJNA1366375), so held-out counts are one less.

| Exclusion level | Multiome ATAC | scATAC |
|---|---|---|
| Examined runs excluded (the coordinator's rule) | 50 studies, 92 runs | 11 studies, 20 runs |
| Every examined study excluded (stricter) | 44 studies, 80 runs | 5 studies, 8 runs |

## 5. Risks and open questions

- **Orientation:** both occur for scATAC on the NovaSeq 6000. tchem must score both orientations and both ends of 24/26-bp reads.
- **Missing barcodes:** 85% of runs are 2-file with the barcode lost. The pool holds only 2025-08 onward, mostly NovaSeq, with almost no older NextSeq or HiSeq runs.
- **Thin scATAC pool:** the scATAC pool is small, and my research examined much of it.
- **26-bp reads** (2 Multiome runs): spacer position is **to verify**.
- **DNBSEQ-T7:** one candidate was sequenced on this instrument. Its i5 orientation is **to verify**.
- **Repo mislabel:** `bench/arc/scripts/label_pool.py` maps "737K-arc-v1" to `10x_arc_v1_gex_label_only.txt`, which is the ARC **ATAC** list, not GEX. This may be a mislabel; reported only, nothing changed.

## Sources
- [G1] CG000168 Rev A (Teichlab mirror): https://github.com/Teichlab/scg_lib_structs/blob/master/data/10X-Genomics/CG000168_ChromiumSingleCellATAC_ReagentKits_UserGuide_RevA.pdf
- [G2] CG000496 Rev A: https://github.com/Teichlab/scg_lib_structs/blob/master/data/10X-Genomics/CG000496_Chromium_NextGEM_SingleCell_ATAC_ReagentKits_v2_UserGuide_RevA.pdf
- [G3] CG000338 Rev F: https://github.com/Teichlab/scg_lib_structs/blob/master/data/10X-Genomics/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevF.pdf
- [S1] cellranger-atac 1.2 source: https://github.com/10XGenomics/cellranger-atac (`lib/python/constants.py`, `mro/atac/stages/processing/setup_chunks/__init__.py`, `attach_bcs/__init__.py`)
- [S2] https://raw.githubusercontent.com/10XGenomics/cellranger-atac/main/lib/python/barcodes/737K-cratac-v1.txt
- [S3] https://github.com/10XGenomics/cellranger-atac/blob/main/LICENSE
- [T1] https://github.com/Teichlab/scg_lib_structs/blob/master/docs/source/epi/10x_scATAC-seq.md
- [T2] https://github.com/Teichlab/scg_lib_structs/blob/master/docs/source/multi/10x_multiome.md
- [T3] https://github.com/Teichlab/scg_lib_structs/blob/master/data/10X-Genomics/737K-cratac-v1.txt.gz
- [T4] https://github.com/Teichlab/scg_lib_structs/tree/master/data/10X-Genomics (`atac_737K-arc-v1.txt.gz`, `gex_737K-arc-v1.txt.gz`)
- [C1] https://github.com/COMBINE-lab/simpleaf/blob/main/resources/permit_list_info.json
- [W1] search snippet only, not fetched (403): https://www.10xgenomics.com/support/software/cell-ranger-atac/latest/algorithms-overview/algorithms-overview
- Local: `$TUSSL_CHEMISTRY_DATA/cellranger/cellranger-10.1.0.tar.gz` (`lib/python/cellranger/barcodes/737K-arc-v1.txt.gz`, `lib/python/cellranger/chemistry_defs.json`)
- ENA portal API: https://www.ebi.ac.uk/ena/portal/api/search

## Runs examined during research
Read content beyond lengths was seen for these runs. "Prefix" means the first 30 bp of the first read in each file was printed; "BC" means the barcode read (64 KB) was matched against the lists in both orientations. None of these runs should be held out. The stricter pool estimate excludes every study listed here.

DRR394953 (PRJDB13870) prefix
SRR10315838 (PRJNA578617) prefix
SRR11692088 (PRJNA630477) prefix+BC
SRR12594185 (PRJNA661557) prefix
SRR13702718 (PRJNA701673) prefix
SRR15679287 (PRJNA759228) prefix+BC
SRR20711371 (PRJNA647723) prefix
SRR22324498 (PRJNA902762) prefix
SRR22988542 (PRJNA762175) prefix+BC
SRR23637209 (PRJNA939339) prefix
SRR23803112 (PRJNA943238) prefix
SRR23866829 (PRJNA944664) prefix
SRR24107517 (PRJNA953340) prefix
SRR24686814 (PRJNA975137) prefix
SRR24766995 (PRJNA977704) prefix
SRR25434104 (PRJNA998795) prefix
SRR26963284 (PRJNA1045639) prefix
SRR29818502 (PRJNA1135080) prefix
SRR30171295 (PRJNA1145873) prefix
SRR31634734 (PRJNA1195001) prefix+BC
SRR31999832 (PRJNA876286) prefix
SRR32259282 (PRJNA1220064) prefix+BC
SRR32763893 (PRJNA1237904) prefix+BC
SRR32948987 (PRJNA1245282) prefix
SRR33916982 (PRJNA1274482) prefix
SRR35330318 (PRJNA1322435) prefix+BC
SRR36147346 (PRJNA1368583) prefix
SRR36278221 (PRJNA1372762) prefix+BC
SRR36481951 (PRJNA1389920) prefix+BC
SRR37808349 (PRJNA1443547) prefix
SRR38019679 (PRJNA1451081) prefix+BC
SRR38154313 (PRJNA1454847) prefix+BC
SRR38742674 (PRJNA1468008) prefix+BC
SRR39152339 (PRJNA1478579) prefix+BC
SRR39379634 (PRJNA1483640) prefix+BC
SRR39463065 (PRJNA1490244) prefix+BC
SRR39525705 (PRJNA1492128) prefix+BC
SRR39607136 (PRJNA1494556) prefix+BC
SRR39819763 (PRJNA63443) prefix+BC
SRR39820076 (PRJNA63443) prefix
SRR39820248 (PRJNA63471) prefix+BC
SRR8893798 (PRJNA532774) prefix
