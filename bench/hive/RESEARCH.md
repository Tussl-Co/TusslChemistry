# Honeycomb HIVE: research, deferred

## What the design runs showed (9 Oct 2026)

Design runs SRR32339204 (CLX, PRJNA1223471) and SRR30495558 (v1, PRJNA1154687):
- Read 1 is a 12-base barcode at bases 1-12, then 14 random bases (a UMI or random linker) at
  13-26, then polyT from about base 27.
- There is no public barcode list and no fixed linker.

**Deferred.** Without a list or a linker, the most a HIVE layout could reach is MEDIUM. Its
segments (barcode, random bases, polyT) are also what a 26-base 10x read looks like, so the layout
would be a clean rival on 10x runs and could take HIGH away from correct 10x calls. HIVE waits on
the same vouching check as SeekGene MM, or on a published list.

> Research by a metadata-only agent on 9 Oct 2026, kept as written. Its tables are here; its
> scripts, downloads and vendor documents stay out of the repository (they read local downloads,
> and some documents may not be redistributed). Paths it names under `dl/`, `raw/`, `docs/`,
> `src/` or `scripts/` are not committed.


## Agent research: Honeycomb HIVE (v1, CLX, "v2"): research notes

Researched 2026-10-09. Metadata and vendor documents only. No FASTQ, BAM or other sequence data was downloaded or opened, and tchem was not run.

## Answer in one line

No public primary source defines the HIVE R1 segment positions. Honeycomb's documents give only these facts:
- R1 has 25 cycles and is the "Cell Barcode" read.
- The cell barcode is 12 bases (BeeNet guide).
- The bead oligo order is UPS, Cell Barcode, Random Linker Sequence, poly(dT), with no lengths.
- No barcode list is published.

Illumina's demo page adds that polyT starts at R1 base 27, for an unstated HIVE version. The barcode start, the linker or UMI length and any fixed sequence in R1 must therefore be read from the design runs.

## Files

| file | content |
|---|---|
| `structure.tsv` | Every structural statement found, one row per source, quoted verbatim. Columns: version, barcode_read, segments (`pos:len:kind`, `?` = not stated), linkers, list, source_url, quote. |
| `runs.tsv` | 775 HIVE runs in 21 studies. Columns: study, run, version_stated, version_source, found_by, file count and names, first_public, instrument, spot length from metadata (base_count/read_count), library kind (cDNA/hash), whether the barcode read is in the FASTQ, protocol window. |
| `design_proposal.tsv` | One design run per version, with the rule, the label quote, and the studies that were skipped and why. |
| `lists.tsv` | Not written: no barcode list was found (see below). |
| `scripts/runs.py`, `scripts/scan.py` | Build runs.tsv from `raw/fr/*.tsv` (ENA filereports) and scan the ENA protocol query. |
| `raw/` | ENA search and filereport TSVs, the Zendesk article index, GEO SOFT for GSE329418, PMC full-text XML. |
| `docs/` | Vendor PDFs as downloaded, with their pdftotext output and rendered figure pages. Proprietary: kept for reference only, do not redistribute. |

## 1. Read structure: what the sources say

### Vendor documents (Honeycomb Zendesk, public Help Center API)

**HIVE scRNAseq v1 Processing Kit User Protocol, Rev A (v22.09, Sept 2022), Appendix 2**
- Read 1: "HIVE Read 1 Seq Primer | Target: Cell Barcode | Cycles: 25". Index 1 and Index 2: 8 cycles each. Read 2: "insert", 50 cycles.
- Custom primers are required.

**Getting Started with HIVE scRNAseq v1 (Dec 2023), p.6, "Follow the molecule"**
- Bead oligos: "Universal Primer Sequence (UPS); Cell Barcode; Random Linker Sequence; Transcript Capture Sequence - poly(dT)". Where the custom R1 primer binds is not stated.
- "Read 1 for cell barcode, and Read 2 for transcript identity".
- No lengths are given.

**HIVE CLX Transcriptome Recovery & Library Preparation User Protocol v1.0 (May 2023)**
- Appendix 3 has the same table as v1: "CLX Read 1 Seq Primer | Cell Barcode | 25".
- Appendix 8 shows the same bead oligo diagram as v1.

**BeeNetPLUS v1.0.X Software Guide (Terra; article updated 2023-11-30), Figure 1.** Version not stated: the text names the "HIVE scRNAseq Processing Kit", and its example data is CLX.
- Labels: "Custom Read 1 Sequence | Cell Barcode | Random Linker | cDNA ... Read 1 (25 bases) Read 2 (50 bases)".
- Text: "Read 1 yielding the individual cell barcode and Read 2 yielding the mRNA sequence".

**BeeNet v2.X.X Software Guide (©2023 Revvity; article updated 2025-01-21), for "the HIVE line of products"**
- Figure 1 labels: "Custom Read 1 Sequence | Cell Barcode | Poly T | cDNA". No linker is drawn.
- Output descriptions: "Barcode: unique cell barcode of 12 bases"; TCM headers: "unique transcriptome barcode of 12 bases".
- `bcin`: "extracts the cell barcode from the forward read (R1) and encodes it in the name of the reverse read (R2)".
- No v1-era BeeNet (v1.1.x) guide text was found. The 12-base length therefore comes from v2-era documents. These say they cover the whole HIVE line, but they do not name v1.

**Vendor FAQ (https://honeycomb.bio/faqs/)**
- "R1: 25  R2: 50  I1: 8  I2: 8"

### Disagreements between sources (both sides kept in structure.tsv)
- **R1 length.**
  - Vendor: 25 cycles (protocols, FAQ, BeeNetPLUS figure, sequencer chart).
  - Illumina BaseSpace demo page: "72 + 8 + 8 + 50" and "polyT (starting from 27bp in R1 read)". This is a NextSeq 2000 run, version not stated.
  - Under the vendor recipe, polyT (base 27 or later) would not be in R1 at all.
- **What follows the barcode.**
  - BeeNetPLUS figure: "Random Linker".
  - BeeNet v2 figure: "Poly T".
  - The v1 and CLX bead diagrams have both, in the order Cell Barcode, Random Linker, poly(dT).

### Not found anywhere
- A UMI. No HIVE source uses the term, although third-party papers report "UMI counts" from BeeNet output.
- The start position of the barcode.
- The length of the random linker.
- Any fixed linker sequence in R1.
- The UPS or custom primer sequences.

The arithmetic 12 + (random linker) = 26 bases before the polyT at base 27 is **not** sourced. It is an open question for the design runs. Do not encode it.

### Third-party tools and specs, checked and not found
- **scg_lib_structs** (Teichlab): has `SeqWell_S3.html` and no HIVE entry.
- **seqspec**: no HIVE spec at pachterlab/seqspec `0771046233fdf93542253c28ecda75054b730080`. The only "hive" hits are the word "archive".
- **GitHub code search**: no STARsolo, simpleaf or kallisto configuration for HIVE.
- **FGCZ SUSHI wrapper** (uzh/sushi `b3a6465865020ce1d62313da100412b02c1d393d`, `lib/otherApps/BeeNetApp.rb`): only runs BeeNet 1.1.2 with `numBarcodes`. Its `--qual-bc` mention gives no structure.
- **De Simone 2025 code** (danledinh/scRNAseq-repo `4494fee15adcade431c051e335ef7a02c5183ed3`): runs BeeNet 1.1.3 and defines no structure.
- **Docker Hub** has unofficial images of the proprietary BeeNet binary: dolphinnext/beenet, nephantes/beenet, briglo/beenet, briglo/beenet_v2 and myjbio/beenet, all with empty descriptions. They were not pulled. They are redistributions of licensed software, not a published definition.
- **PRJNA1078067** calls HIVE "a commercial Seq-Well platform". This is a submitter's claim. It does not carry over Seq-Well's 12+8 structure.

### Papers
- **De Simone et al. 2025 NAR** (10.1093/nar/gkae1186, PMC11754665):
  - Kit: "HIVE CLX scRNA Seq Kit v1 (Honeycomb Biotechnologies)".
  - Sequencing config: "[25, 50, 8, 8]".
  - Processing: BeeNet 1.1.3.
  - Data: PRJNA1106903. The Honeycomb runs SRR28867559 and SRR28867560 are single-file with a 50-base spot, so the barcode read is probably missing.
- **ABRF multisite study** (PMC13252917; preprint 10.1101/2025.10.24.684427; GSE329418 = PRJNA1458839; GEO contact F. W. Kolling, Dartmouth). This is the "Lee et al. 2025 (GSE329418)" named in the task:
  - "V1 kit Rev A was used for leukocytes and CLX for the extension study with the PBMCs".
  - Its table labels the PBMC library "HC HIVE v2".
  - Sequencing: "25-8-8-50".
  - GEO labels GSM9703649 "chemistry V1" and GSM9703650 "chemistry V2".
  - **"v2" exists only in these two places, and the paper's text equates it with CLX.** No Honeycomb document names a v2.
- **Clinical neutrophil comparison** (PMC12539239): "using HIVE custom primers and 25-8-8-50 sequencing cycles".

### Catalogue number 10347
- PRJNA1107957 cites "cat 10347 v.0.3 revA".
- No Honeycomb document or web page found names 10347.
- The part numbers that do appear on Zendesk COA/COC/SDS files are 10134, 10135, 10141, 10155, 10156, 10266-P, 10271, 10369, 10370, 10418, 10429, 10479 and 10922.
- **Not identified.** PRJNA1107957 is therefore counted as version unstated.

## 2. Barcode list

**None found.** BeeNet calls cells by rank rather than against a list:
- "--num-barcodes ... after sorting by count of reads mapped to genes, as the number of cell barcodes in the sample".
- "--hamming-bc=#: Barcode hamming reduction distance, 0 or 1 ... combines cell barcodes that are mismatched by one base pair".
- No whitelist option is documented.

The same holds in BeeNetPLUS (`NumBC`) and in De Simone ("Honeycomb requires users to input the expected cell number and returns exactly this number of barcodes"). This points to an `Unlisted` 12-base segment judged by concentration, as with Drop-seq.

**Terms**: the BeeNet guide says "This manual is proprietary to Honeycomb Biotechnologies, Inc., ... shall not be used or distributed for any other purpose without the prior written consent". BeeNet itself is downloaded after accepting a licence (https://download.honeycomb.bio/).

## 3. Public runs

### Queries
- ENA portal API, `result=read_run`, `library_construction_protocol="*honeycomb*"` (664 runs) and `"*hive*"` (2,704 runs).
- Also `"*beenet*"` (0), `"*picowell*"` (235) and `"*pico-well*"` (0).
- A local regex `\bHIVE\b|honeycomb|beenet|pico-?well` (case-insensitive) kept 833 runs in 22 studies.

### Removed as false hits
- PRJEB3411, PRJEB1945 and PRJNA479740: bee hives and honeycomb cells.
- PRJNA1225645: "single-cell sorting honeycomb plate", processed with the BD Rhapsody pipeline.

### Added from outside the protocol query
These would not enter a catalogue drawn from protocol text:
- PRJNA1106903 (De Simone): empty protocol field. Library names Honeycomb_1/2.
- PRJNA1085264 (equine BAL, PMC11760581): empty protocol field. Experiment title "HIVE scRNA-seq equine BAL".
- PRJEB67548: study title "Single cell sequencing with HIVE and HIVE CLX of K562 and mESC cells", empty protocol field.
- Rejected: PRJNA1107508, whose "HIVE" is the study's own acronym (a CITE-seq workflow).

### GEO
E-utilities `gds` search for HIVE plus Honeycomb, "HIVE CLX" and BeeNet. Every non-bee hit maps to a study already above: GSE308157, GSE329418, GSE276100, GSE273791, GSE266615 and GSE277953.

### Multi-kit studies
Only HIVE runs are listed:
- PRJNA1107957: samples named `*Hive*`.
- PRJNA1458839: HC-11 to HC-20, HC-V1-S1 and HC-V2-S2.
- PRJNA1106903: Honeycomb_1/2.

### Study counts (21 studies, 775 runs; 20 studies have paired runs)

| version (how it is stated) | studies | with paired runs | studies |
|---|---|---|---|
| v1, literal in archive metadata | 2 | 2 | PRJNA1154687 (ENA "HIVE scRNAseq V1 Processing workflow"); PRJNA1458839 HC-V1-S1 (GEO "chemistry V1") |
| v1, indirect | 1 | 1 | PRJNA1110737 ("October 2021 protocol (v21.10)", the v1 kit's launch revision) |
| v1, paper only | 2 | 2 | PRJNA1085264 (PMC11760581); PRJNA1458839 HC-11..20 (ABRF table "HC HIVE v1") |
| CLX, literal in ENA | 3 | 3 | PRJNA1223471, PRJNA1425141, PRJNA1164800 (DNBSEQ-G400) |
| CLX, by protocol URL | 1 | 1 | PRJNA1257368 |
| CLX, paper only | 2 | 1 | PRJNA1134095 (PMC12234664); PRJNA1106903 (single-file only) |
| v2 | 1 | 1 | PRJNA1458839 HC-V2-S2 (GEO "chemistry V2"; the paper says CLX) |
| conflicting | 1 | 1 | PRJNA1078067 ("Hive CLX" plus "HIVE scRNA-seq Processing Kit revA", the v1 kit title; 15,000 cells, which is v1's recommended load) |
| mixed, no per-run label | 1 | 1 | PRJEB67548 (v1 + CLX) |
| unstated | 9 | 9 | PRJNA973610, PRJNA820906, PRJNA1107957, PRJNA1143064, PRJNA1241229, PRJNA1243391, PRJNA1115046, PRJNA1285791, PRJNA1329045 (weak: "Honeycomb scRNA isolation") |

PRJNA1458839 appears in three rows.

### Barcode read in the FASTQ
Judged from ENA file counts and names only:
- 2-file `_1/_2` runs: barcode read probably present. This covers every study except the two below.
- PRJNA1106903 Honeycomb runs (SRR28867559, SRR28867560): single file, barcode read probably missing.
- PRJNA973610: SRR24614903 has no FASTQ.
- PRJNA1085264: SRR31657929 and SRR31657930 also have an unpaired third file.

### Spot lengths
Spot lengths from metadata (base_count/read_count) vary: 75, 73.2, 76, 100, 116, 122, 300 and 38 (PRJEB67548). 75 matches 25+50. The others suggest longer reads or trimming, but read lengths were not measured here.

PRJNA1110737 has 27 cell-hash libraries (`library_kind=hash`).

## 4. Design proposal (metadata only; reads not opened)

**Rule**:
- Pick the earliest-published study whose archive metadata names the version unambiguously for the run.
- Within that study, take the earliest `first_public`, then the lowest run accession, among 2-file cDNA runs.
- Studies that are skipped are listed in design_proposal.tsv with quotes, so the caller can override the rule.

| version | run | study | why |
|---|---|---|---|
| v1 | **SRR30495558** | PRJNA1154687 (2025-06-09, NextSeq 550) | ENA: "following the HIVE scRNAseq V1 Processing workflow". The lowest accession overall would be SRR30495556 (first_public 2025-06-23). Skipped earlier studies: PRJNA1110737 (v21.10 only, with hashing), PRJNA1085264 (paper only), PRJEB67548 (mixed). |
| CLX | **SRR32339204** | PRJNA1223471 (2026-01-16, NovaSeq X) | ENA: "a designated HIVE CLX collector". Skipped earlier studies: PRJEB67548 (mixed), PRJNA1106903 (paper only, single-file), PRJNA1078067 (conflicting), PRJNA1134095 (paper only). |
| v2 | **SRR38306962** | PRJNA1458839 (2026-06-20, NovaSeq 6000) | GEO GSM9703650 "chemistry V2". The paper says this library is CLX, so v2 may not be a separate layout. If v2 shares a layout entry with v1 or CLX, Rule 2 removes all of PRJNA1458839 from held-out, including its 22 v1 runs. |

PRJEB67548 was checked for per-run version labels on 2026-10-09 (ENA experiment XML and BioSamples). Its samples are titled I23-1183-0x / I23-1314-0x, with organism mouse, human or mixed, and no version label, so it stays skipped.

What each design run can show, judged from metadata spot length:
- The v1 run (73.2) and the v2 run (76) fit 25+50. Their R1 probably ends near base 25, so they can place the barcode but probably not the end of the linker or a polyT at 27.
- The CLX run (300) probably can show all three.
- No run with a literal v1 label has a longer spot length. The v1-by-paper study PRJNA1085264 has spot length 116.

What the design runs must settle:
- Where the 12-base barcode starts.
- The length and character of the random linker (UMI-like or not).
- Whether polyT starts at base 27 (only visible where R1 is longer than 26).
- Whether v1, CLX and v2 differ at all.

## Sources (pins)

Zendesk attachments are pinned by sha256 as downloaded, with the article's `updated_at`. List API: `https://honeycombbio.zendesk.com/api/v2/help_center/en-us/articles.json` (sha256 4e6ea0f4…).

| doc | URL | article updated | sha256 |
|---|---|---|---|
| v1 Processing Kit protocol Rev A | https://honeycombbio.zendesk.com/hc/article_attachments/21351113668891 | 2023-12-18 | 57fc4d6751a53ac53e1f1f1a931b72150156f6cfbf1778869bb7cc799bf3fa49 |
| Getting Started v1 | https://honeycombbio.zendesk.com/hc/article_attachments/21350616821787 | 2023-12-18 | 1c7cf4e7c4ab876fe6432cf5de415ac65fc426731e72ab7ef95a550b14a4e4db |
| v1 Sample Capture protocol Rev A | https://honeycombbio.zendesk.com/hc/article_attachments/21351300869659 | 2023-12-18 | 36b7adbd15a5a9aa81dd165836c49c5e05636267d49e968cb2893128896c6355 |
| CLX Library Prep protocol v1.0 | https://honeycombbio.zendesk.com/hc/article_attachments/20582312953499 | 2023-11-21 | 453f643c4e10ccc12337604d8ff4cb91813a7aec4bf3f3c1bc2088e35ad71d76 |
| CLX Sample Loading protocol | https://honeycombbio.zendesk.com/hc/article_attachments/20581584307483 | 2023-11-21 | 6cc5dabe3039e796d4180d7155555a6fba30facd873d51af74e9f25e0345d09c |
| BeeNet v2.X.X guide | https://honeycombbio.zendesk.com/hc/article_attachments/33024529416091 | 2025-01-21 | b4ef72c8f1c6e6837abc5a4a7a11bda679c61b5d8cf31e4527b3e5071ecf921c |
| BeeNetPLUS v1.0.X guide | https://honeycombbio.zendesk.com/hc/article_attachments/20811926861339 | 2023-11-30 | d641fcddfe7c846b41908e00e17d6a6ecd7e032576c5fb2851adf1a5191b4fac |
| Sequencer compatibility chart (v1) | https://honeycombbio.zendesk.com/hc/article_attachments/10601667205147 | 2023-04-05 | abe0d2a6b66d83d7806787762124551f646a0a280916c0ed02bf5c937279b6fa |
| CLX troubleshooting | https://honeycombbio.zendesk.com/hc/article_attachments/15173423438747 | 2023-05-04 | b2710fb172946f6fc4cf6dce81e0d6352c11c3fd6d67bef0669a935bc703f9d1 |

Other sources:
- https://developer.illumina.com/news-updates/hive-single-cell-rnaseq-libraries-data-from-honeycomb-biotechnologies-are-now-available-on-illumina-basespacetm-sequence-hub-demo-data-page
- https://honeycomb.bio/faqs/

Not fetched, because they contain sequence data: `https://1226200508.rsc.cdn77.org/HIVE_CLX_InputDS_fastq.tar.gz` and `..._BeeNet.tar.gz`, linked from Zendesk "Example Data for BeeNetPLUS". The test-data HTML report and sample sheet were fetched; they contain QC metrics only.
