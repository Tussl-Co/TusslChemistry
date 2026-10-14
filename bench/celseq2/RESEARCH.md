> Research notes written before this kit's catalogue was drawn (8 Oct 2026), kept as the record of
> sources and of which runs' reads were examined (those runs' studies are never held out: see
> results/examined.tsv). Scripts these notes name were scratch work and are not kept; the
> catalogue's rules are the ones in scripts/select.py, which say where they differ from these notes.

# CEL-seq2 family (CEL-seq2, SORT-seq, MARS-seq): catalogue research, 8 Oct 2026

## 1. Read structure

| kit | barcode read | other read | source |
|---|---|---|---|
| CEL-seq2 (Yanai) | R1: UMI 6, cell barcode (CB) 6, then polyT (24 T, then V) | R2: cDNA | Supplement MOESM4 primer `...GACGATC NNNNNN [6 bp BC] T24V`; Teichlab page; `bc_demultiplex.py` slices UMI first (`umi_strt=0`, then the BC) |
| CEL-seq2 protocol run | R1 15 bases, index 7, R2 36 | | MOESM4 ("15 bases for read 1 ... 36 bases for read 2") |
| SORT-seq, Muraro 2016 | R1: "24 bp polyT, 4 bp UMI, 8 bp CB", "256 possible UMIs", 2×75 on a NextSeq | R2: cDNA | Muraro full text (PMC5092539) |
| SORT-seq / "CEL-seq2", Utrecht labs (current) | R1: UMI 6 then CB 8 (`--lenumi 6 --umifirst`, 384 8-mers), then polyT | R2: cDNA | anna-alemany/transcriptomics `extractBC.sh`, `bc_celseq2.tsv`. The Muraro 4+8 layout vs the later 6+8 layout is **to verify** on a design run |
| MARS-seq (2014) | R1: 4 bp pool barcode + 3 N, then cDNA | R2: CB 6, UMI 4, then dT | Teichlab MARS-seq page |
| MARS-seq2.0 | R1: 4 bp pool barcode + 5 N, then cDNA | R2: CB 7, UMI 8, then dT | Teichlab page; Keren-Shaul 2019 Supp. Tables 1 and 2 |

All CEL-seq2/SORT-seq barcode reads start straight at the UMI: there is no fixed linker before the polyT. MARS-seq splits its barcodes across both reads. R2 still has CB+UMI at fixed offsets, so the layout fits "fixed segments on one read" only if the tool accepts the barcode read in file `_2`. tchem already reports `swapped()` (src/main.rs:233), so it does handle a barcode read in `_2`.

Sources: https://teichlab.github.io/scg_lib_structs/methods_html/CEL-seq_family.html, https://teichlab.github.io/scg_lib_structs/methods_html/MARS-seq.html, https://github.com/yanailab/CEL-Seq-pipeline/blob/stable/bc_demultiplex.py, https://pmc.ncbi.nlm.nih.gov/articles/PMC5092539/, https://github.com/anna-alemany/transcriptomics/blob/master/mapandgo/bwamap/extractBC.sh

## 2. Barcode lists

- **CEL-seq2, 96 × 6-mers**: the paper's Additional file 4 (MOESM4 docx, https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-0938-8/MediaObjects/13059_2016_938_MOESM4_ESM.docx) lists 96 primers 1s–96s. The text says "48 primers". The 96 barcodes are identical, in the same order, to `example/barcodes_cel-seq_umis96.tab` in https://github.com/yanailab/celseq2 (BSD-3-Clause). The list has a minimum Hamming distance of 2, as the protocol states. yanailab/CEL-Seq-pipeline `barcode_umis.tab` (GPL-3.0) holds the first 48 of these, but its config says `umi_length = 5`, the older 5+6 "CEL-Seq + UMI" design (Supp. Table S2, MOESM5). Licence of the Genome Biology supplement: BMC articles are normally CC BY 4.0; Crossref gave no licence (**to verify**). Do not use `minoda-lab/universc` `bc_celseq2.txt`: it shares only 72 of the 96.
- **SORT-seq / Utrecht "celseq2", 384 × 8-mers**: `bc_celseq2.tsv` in anna-alemany/transcriptomics (GPL-3.0) and anna-alemany/VASAseq (no licence). The two files are identical. It has 384 unique barcodes with a **minimum Hamming distance of 1**, so matching must be exact. I found no primary publication of this list, so its provenance and licence are a risk. It is not a superset of the Yanai 96: 21 of the 96 appear as its first 6 bases. `bc_celseq1.tsv` (96 × 8, CEL-seq1 4+8) shares no barcode with it.
- **MARS-seq2.0**: Keren-Shaul et al., Nat Protoc 2019, Supp. Table 1 gives 384 RT1 primers with 7-mer CBs and 8 N UMIs. Supp. Table 2 gives 32 pool adaptors (4-mers + 5 N). Mirrored at https://github.com/Teichlab/scg_lib_structs/tree/master/data/MARS-seq. Crossref licence: Springer TDM only, so redistribution is **to verify**. The original MARS-seq lists (Jaitin 2014, Tables S7 and S8) were not checked.

## 3. ENA

The portal's wildcards return 0 rows when the term contains a hyphen (`*CEL-seq*` gives 0). Without the hyphen, `*CELseq*` (5,361 runs) still matches "CEL-seq" text. Queries `*CEL*seq*`, `*SORT*seq*` and `*MARS*seq*`, filtered locally by `cel[-_ ]?seq|sort[-_ ]?seq|mars[-_ ]?seq`, give **9,208 runs in 121 studies**, 7,100 of them paired (2 files). Of the 9,208, 2,484 runs name CEL-seq2, 721 name SORT-seq and 3,345 name MARS-seq.

- **Paired runs**: 72 studies have at least one run with 2 files. 49 studies hold no paired run (almost all single-file), typical of per-cell or per-sample demultiplexed submissions with the barcode stripped. The draw already drops these (e.g. PRJNA732658, 638 runs titled `Lane3_RPI2_11s`, which use the supplement's barcode names; most MARS-seq PRJNA1021xxx studies).
- **Per-cell runs that are paired**: PRJNA432538 (374 runs titled `B_Lymph_A9`, median 123k reads per run, 32/70) and PRJNA343284 (2,684 per-well runs; plain CelSeq). Each file holds one barcode, so concentration-based detection sees a single barcode. List matching still works.
- **Read lengths** (modal, first 64 KB, files `_1`/`_2`, 2 runs per paired study; full table in `dl/shapes2.tsv`):

| study | kit text | `_1`/`_2` |
|---|---|---|
| PRJNA486448 | CEL-seq2 | 63/14 |
| PRJNA481966, 676207, 767348, 769669, 1018873 | CEL-seq2 | 32/71–72 |
| PRJNA686092 | CEL-seq2 | 15/100 |
| PRJNA1065811, 1174744 | CEL-seq2 | 26/60 |
| PRJNA476320 | SORT-seq/CELseq2 | 60/14 |
| PRJNA646712, 801598 | SORT-seq | 63/16 |
| PRJNA834641, 939161, 1429303 (scdiscoveries) | SORT-seq | 26/60 |
| PRJNA474920, 506282 | SORT-seq | 76/76 |
| PRJNA762008 / 1088061 | MARS-seq | 69/15 / 75/15 |
| PRJNA1003003, 1187607, 704688 | MARS-seq | 75/8 |

The 26 and 14–16 base barcode reads fit UMI 6 + CB 8 (+ T). The 32/72 group fits several layouts, and the shape rule cannot tell which. MARS-seq's 8-base `_2` file is probably an index read (**to verify**).

## 4. Proposal

- **Query**: `library_construction_protocol="*CEL*seq*" OR "*SORT*seq*" OR "*MARS*seq*"`, run as three calls and merged, AND `library_strategy="RNA-Seq"`. This filter removes the Hi-C runs of PRJNA625734, which pass on protocol text alone (**to verify** that the field filter works on read_run).
- **Classifier** (first match wins):
  1. Refuse `dam ?& ?t|damid|vasa|hiv|degenerate|sort_seq|reporter|hi-?c\b|pbat|methylom|atac` and `names_other_vendor(...)`. 10x appears in PRJNA507788, 603489 and 1195951.
  2. sortseq: `\bsort[- ]?seq\b|single cell discoveries|scdiscoveries`.
  3. celseq2: `\bcel[- ]?seq[- ]?2\b`.
  4. marsseq: `\bmars[- ]?seq`.
  5. Refuse a bare "CEL-seq": that is v1 (8 bp CB, no UMI) or 5+6.
- **Shape** (order-free; tchem handles swaps): let `lo, hi = sorted(lengths)`.
  - celseq2: `12 <= lo <= 40 and hi >= 35`
  - sortseq: `14 <= lo <= 40 and hi >= 35`
  - marsseq: `14 <= lo <= 40 and hi >= 50`

  This drops equal long pairs (76/76, 151/151), because the shape rule cannot tell whether the barcode survived export. It also drops 6/63 and 75/8.
- **Pool** (classifier + shape applied to the measured runs; a lower bound, since only 2 runs per study were checked, not 6):

| kit | runs | studies | held out | design run |
|---|---|---|---|---|
| celseq2 | 27 | 16 | 26 | SRR7704914, PRJNA486448, 2019 |
| sortseq | 25 | 13 | 24 | SRR7347826, PRJNA476320, 2019 |
| marsseq | 8 | 5 | 7 | PRJNA435744, 2018; excludes PRJNA625734 (Hi-C) |

## 5. Risks and open questions

1. "CEL-seq2" in ENA text covers two layouts: Yanai 6+6 with 96 barcodes, and Utrecht/SCD 6+8 with 384 barcodes. The kit name does not fix the barcode list, and the two kits may need to merge.
2. The Muraro-era SORT-seq layout (4+8) may differ from later runs (6+8).
3. The 384 list has no primary publication or clear licence.
4. Paired per-cell runs (PRJNA432538) hold one barcode each.
5. Bulk CEL-seq/MARS-seq studies (PRJNA743007, 1088061) may pool few barcodes.
6. MARS-seq needs a barcode read in `_2` and a 4 bp pool barcode on the cDNA read. The 75/8 runs probably lost R2.
7. The classifier reads only the first run's protocol, so multi-assay studies (PRJNA603489) can mislead it.

## Runs examined during research

None. Only read lengths were measured: modal length of the first 64 KB per file, the same check the shape rule runs. No sequence content, barcodes or polyT positions were inspected.
