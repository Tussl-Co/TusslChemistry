# Illumina Single Cell 3' RNA Prep (PIPseq V): research and the design run

## What the design run showed (9 Oct 2026, read once, now design data)

The design run is SRR33161101. Read 1 holds:
- the four barcode parts, split by linkers `ATG`, `GAG` and `TCGAG`;
- a 3-base binning index;
- an `ACG` spacer, at offsets 0-3.

These are the linker bases Illumina does not publish. The layout `Illumina Single Cell 3' RNA Prep / PIPseq V` requires no polyT tail, because only 28% of the design run's reads had a clean one.

The held-out PIPseq catalogue's classifier (bench/pipseq) misses studies that name only "V" or
the Illumina kit. Those are left for a recorded classifier revision before any of their runs is
read.

> Research by a metadata-only agent on 9 Oct 2026, kept as written. Its tables are here; its
> scripts, downloads and vendor documents stay out of the repository (they read local downloads,
> and some documents may not be redistributed). Paths it names under `dl/`, `raw/`, `docs/`,
> `src/` or `scripts/` are not committed.


## Agent research: Illumina Single Cell 3' RNA Prep (PIPseq V, "v5"): sources for a layout

Gathered 9 Oct 2026 from metadata and documents only. No FASTQ or sequence data was opened. tchem was not run. No form was submitted.

Files:
- `structure.tsv`: Read 1 and Read 2 elements. Each row has a source URL and a quote.
- `lists.tsv`: barcode lists. There is no public v5 list; the v3/v4 lists are given for reference.
- `runs.tsv`: every candidate run, with a study category.
- `design_proposal.tsv`: one design run.
- `studies_table.md`: one line per study.
- `scripts/`: the scripts that made these files.
- `src/`: downloaded pages, PDFs and figures.
- `ena/`, `geo/`: raw metadata, gzipped.

## 1. Read 1 structure

| bases (0-based, after phase) | element | source |
|---|---|---|
| 0-3 before position 0 | phase bases | DRAGEN v4.6 |
| 0-7 | tier 1, 8 nt, one of 96 | DRAGEN v4.6 |
| 8-10 | linker, 3 nt, sequence not published | DRAGEN figure ("L") |
| 11-16 | tier 2, 6 nt, one of 96 | DRAGEN |
| 17-19 | linker, 3 nt, not published | DRAGEN figure |
| 20-25 | tier 3, 6 nt, one of 96 | DRAGEN |
| 26-30 | linker, 5 nt, not published | DRAGEN figure |
| 31-38 | tier 4, 8 nt, one of 96 | DRAGEN |
| 39-41 | binning index, 3 nt, random. Not a UMI. | DRAGEN; preprint Fig 7 |
| 42-44 | spacer `ACG` | preprint Fig 7A ("BI + Spacer NNNACG") |
| 45- | poly(dT)30V | Illumina KB 9449 figure; preprint |

Read 2 is cDNA, at least 72 nt. Its first base is a constant T, which DRAGEN trims. TSO is trimmed at the 5' end and polyA at the 3' end (DRAGEN). Index reads are 10 + 10 nt.

Read 1 must be at least 45 nt. All four user guides (T2, T10, T20 and T100) say this, and Illumina's 100-cycle recipe is 45+10+10+72. No source shows a structure that differs by size. DRAGEN has one setting for all sizes. **M1: not found** in the training packet, the documentation page or the product page.

Sources:
- DRAGEN v4.6 Illumina scRNA: https://help.connected.illumina.com/dragen/dragen-v4.6/product-guides/dragen-v4.6/dragen-single-cell-pipeline/dragen-scrna-illumina.md
  - "Each read in R1 includes 0 to 3 bases at the start indicating the "phase" of the read, a barcode sequence with 4 tiers, which identifies the individual PIP, followed by a 3-base binning index (BI) sequence. ... Each barcode tier is composed of 96 possible sequences ... By default, the barcode positions are set to 0_7+11_16+20_25+31_38 ... the UMI positions are set to 39_41."
  - "Note that the binning index position is set using the `--scrna-umi-position` option."
  - "IMIs are grouped in one of 64 bins, based on the 3-base binning index."
  - The figure `src/fig_r1.png` labels positions 8-10, 17-19 and 26-30 "L: Linker". DRAGEN v4.4 gives the same positions.
- Fontanez et al., bioRxiv 10.1101/2024.10.04.616561 (Fluent authors): https://www.biorxiv.org/content/10.1101/2024.10.04.616561v1.full
  - Fig 7A legend: "PIPseq V core particles with oligonucleotide sequences carrying 5' PCR handle, cellular barcodes, 3N base binning index, ACG spacer, and poly A mRNA capture moiety."
  - Methods: "v4.0PLUS and UCP libraries utilize a 54-base Read 1 ... PIPseq V (IMI) libraries differ with a shorter 45-base R1".
- Illumina KB 9449: https://knowledge.illumina.com/library-preparation/illumina-single-cell/library-preparation-illumina-single-cell-reference_material-list/000009449.md
  - Its library figure (`src/kb_lib.png`) shows "Barcodes 42-45 | 30Poly(dT)V".
- Training packet 200065159 v01: https://support.illumina.com/content/dam/illumina-support/courses/illumina-single-cell-3-rna-prep-training-packet/illumina-single-cell-3-rna-prep-training-packet-200065159-01.pdf
  - "Illumina Single Cell 3' RNA Prep (formerly known as PIPseq V)"
  - "limited base diversity of Read 1 starting at cycle 46 (poly A region)"
- KB 9448, archived 15 May 2026 (now 404): http://web.archive.org/web/20260515105303/https://knowledge.illumina.com/library-preparation/illumina-single-cell-rna-prep/library-preparation-illumina-single-cell-rna-prep-reference_material-list/000009448
  - Title: "Differences between V4.0PLUS PIPseq and Illumina Single Cell (formerly PIPseq V)"
  - "Intrinsic Molecular Identifiers rather than UMIs"
  - "10 bp indexes"
- Fluent PIPseq V user guides (T2, T10, T20), Rev 1.5, mirrored by NCI, for example https://crtp.ccr.cancer.gov/wp-content/uploads/fsgc-protocols/PIPseq-V-T2-3-Single-Cell-RNA-Kit-User-Guide-1-4.pdf
  - Read 1 at least 45 and Read 2 at least 72.
  - 10-base UDIs.
  - Appendix II lists i5/i7, WTA, TSO and the ligation adapter. It lists **no bead oligo**, so the linker sequences are not published.

One conflict, left open:
- The ACG spacer is in the preprint but not in DRAGEN's positions or KB 9449's figure.
- The training packet's "cycle 46 (poly A region)" matches the spacer at phase 0 (3 + 39 + 3 = 45 bases before polyT).
- The design run's reads settle it.

## 2. How v5 differs from v3/v4, and whether the same reads fit both

- **Same:** the cell barcode region. There are 0-3 phase bases and four tiers of 8/6/6/8 at the same offsets, and linker gaps of 3/3/5 at the same places as ATG/GAG/TCGAG.
  - The v5 linker sequences are not published anywhere I found, so ATG/GAG/TCGAG are unconfirmed for v5. The design run must settle this.
- **Changed:**
  - The 12-nt UMI is replaced by a 3-nt binning index and the fixed `ACG`. polyT starts at 45 + phase instead of 51 + phase.
  - The kit has no UMI. Molecules are counted by IMI (R2 fragment start) within 64 BI bins.
  - The 54-nt R1 of v4.0PLUS becomes 45 nt.
  - Index reads become 10 nt. v3 used 8; the preprint says v4.0PLUS used "optional 8 or 10 base indexes".
- **Can one read fit both layouts? By position, in one direction only.**
  - Under the v3/v4 layout, a long v5 read places NNN + ACG + TTTTTT in the 12-nt UMI slot (positions 39-50), with polyT after it. By position alone, it fits.
  - Even so, the current tool called every v5 run "no known layout" (`bench/pipseq/results/results.tsv`). In practice the v3/v4 layout does not accept them now.
  - A v3/v4 read under a v5 layout fails: positions 42-44 are random UMI bases, and ACG appears there by chance only 1 time in 64.
  - The deciding bases are positions 42-50 after phase. v5 shows `ACG` then T; v3/v4 shows random bases.
  - These bases are visible only when R1 is longer than 45 + phase. On Illumina's own 45-cycle R1, a phase-3 read ends at the BI and a phase-0 read ends at ACG, with no polyT at all. A v5 layout must not require polyT on 45-46 nt R1s.
  - R1 length itself is a strong hint: 45-46 nt for v5, 54 nt for v4.0PLUS.

## 3. Barcode lists

**No public v5 or Illumina list was found.**
- DRAGEN ships `<PIPSEQ_BARCODE_SEQ_LIST>` inside the software: "packaged with DRAGEN".
- PIPseeker sits behind a form.
- DRAGEN states 4 tiers x 96.
- The public v3 tier lists are 96 each, 8/6/6/8 nt. They are identical in pipspeak (MIT), scg_lib_structs (no licence stated) and MethodsDev (MIT). scg_lib_structs says v3 = v4.
- No source says v5 reuses them. Hashes and pinned URLs are in `lists.tsv`.
- Unless the design run shows the v3 tiers covering v5 reads, the layout should use unlisted barcode concentration, as v3/v4 does now.

## 4. Public runs

**ENA search.** ENA portal `read_run`.
- These patterns return 0 rows from ENA: `"*PIPseq*V*"`, `"*PIPseq*v5*"` and `"*Illumina*Single*Cell*3*"`. So do `"*PIPseq V*"`, `"*Single Cell 3*"` and `"*Illumina*Single*"`.
- Other patterns with inner wildcards work: `*RNA*Prep*` returns 10,687 runs. The cause of the failures is not established.
- PRJNA1269311's ENA protocol contains "Single Cell 3' RNA Prep", yet no ENA pull returned it; GEO found it. **All counts below are lower bounds.**
- I therefore pulled broad sets and filtered them locally: `*PIPseq*` (868 runs), `*PIP*seq*`, `*Fluent*` (12,641), `*RNA*Prep*` (10,687), and titles.
- I added ENA runs of the BioProjects that GEO names.

**GEO search.** E-utilities `gds`, series only. I fetched SOFT text (brief) for 208 series and classified it with the same regexes. Series that matched only the apostrophe-free phrase `"Single Cell 3 RNA Prep"` (571 hits, mostly 10x) were not fetched.

**Classifier** (`scripts/classify.py`):
- A study counts as v5 if it has a v5 token in PIPseq/Fluent context. Matching forms include "PIPseq V", "PIPseq™ V", "PIPseq-V", "PIPseqV5", "Fluent V", "v5.0" and "Pipseek V5".
- It also counts if it names the Illumina kit ("Single Cell 3′ RNA Prep").
- A study is flagged when it also names a v3/v4 kit (PIPseeker's own "v3.3" is excluded) or another vendor's words (the repo's `VENDORS`).

**22 studies:**

| category | studies | with paired ENA FASTQ |
|---|---|---|
| single kit (v5 / PIPseq V and/or Illumina name, nothing else) | 13 | 12 (PRJNA1304968 has no FASTQ yet) |
| mixed version (names V and v4.0PLUS, or Illumina and v4.0PLUS) | 3 (PRJNA1216192, PRJNA1309089, PRJNA1373237) | 3 |
| ambiguous: a v4 kit called by the Illumina name | 1 (PRJNA1472464) | 1 |
| other vendor named | 5 (PRJNA1235201, PRJNA1199038, PRJNA1255131, PRJNA1330391, PRJNA1308940) | 5 |

Notes on the other-vendor studies:
- PRJNA1199038's 10x word is a nuclei-isolation kit, not a library kit.
- PRJNA1235201 and PRJNA1255131 have run titles naming PIPseq V ("Illumina PIPseq V" and "Fluent V").

**The 13 single-kit studies,** by `kit_name`. The name is taken from the study's own ENA and GEO text by regex (`scripts/categorise.py`), never set by hand:
- "v5", with a numeral 5 (5): PRJNA1251639, PRJNA1252110, PRJNA1252747 (T2, "Kit V5"); PRJNA1268916 (T20, "Pipseek V5"); PRJNA1304968 (title "PIPseqV5", no FASTQ yet).
- "PIPseq V", no digit (3): PRJNA1284893 (T20, "PIP-Seq V"); PRJNA1445686 (T2); PRJNA1477052 (T10, "PIPseq-V").
- Illumina name only (4): PRJNA1269311 (T20); PRJNA1438756 (T2); PRJNA1481985 (T10); PRJNA1480533 (no size).
- PIPseq V and Illumina (1): PRJNA1377937 (T10).

By size: T2 5, T10 3, T20 3, none stated 2. Every size token was checked to sit next to the kit name. No T100 or M1 study was found. Protocol windows are in `studies_table.md` and `runs.tsv`. Illumina catalogue numbers seen in protocols: #20135689 (T2), #20135691 (T10).

Read shape (mean spot length = bases / spots):
- 2x150 (spot length about 300-302): the most common shape.
- 45x72 style (about 118): PRJNA1377937 and PRJNA1235201.
- 2x100 (about 200-202): PRJNA1284893 and PRJNA1480533.

## 5. Design proposal (metadata only)

**PRJNA1251639 / SRR33161101**:
- PIPseq T2 Kit V5, first public 18 Apr 2025.
- Spot length 302 (2x151), so R1 reaches the spacer and polyT.
- R1: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR331/001/SRR33161101/SRR33161101_1.fastq.gz
- R2: https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR331/001/SRR33161101/SRR33161101_2.fastq.gz

Two earlier studies were skipped:
- PRJNA1216192 (29 Jan 2025) names "PIPseq™ V T10" and a "T20 ... v4.0PLUS" kit.
- PRJNA1235201 (14 Mar 2025) names 10x GEM-X.

This is the same design row the repo already holds in `bench/pipseq/results/catalogue.tsv`. Under VALIDATION §5 rule 2, PRJNA1251639 is a design study. PRJNA1252110 and PRJNA1252747 were already called in `results.tsv` (as held-out v5), so a v5 held-out record must come from the other single-kit studies. Not proposed: the earliest study using only the Illumina name is PRJNA1269311 (T20, 18 Jul 2025). Naming it as a design study would spend it under rule 2.

## 6. Findings about the repo's records (read only; nothing changed)

1. `bench/pipseq/results/catalogue.tsv` lists **PRJNA1216192 and PRJNA1373237 as `pipseq_v4`**, but both name a V/Illumina kit as well as v4.0PLUS.
   - Cause: `select.py`'s token reads "v4.0PLUS" but not "PIPseq V" without a digit, nor the Illumina name.
   - PRJNA1373237: its R1 is 46 nt in `results.tsv`. That fits the v5 recipe (45), not v4.0PLUS (54). Its "no known layout" calls are not separate evidence, because a 46-nt R1 cannot meet the v3/v4 layout's length.
   - PRJNA1216192: its two catalogued runs were called v3-v4 HIGH, unlike every v5 run. That fits v4.
   - Under VALIDATION §5.1 these calls are not relabelled. A fixed rule applies only to studies not yet seen.
2. PRJNA1472464 (`pipseq_v4` in the repo) calls its v4 kit "now the Illumina Single Cell 3'RNA Prep Kit". Its 54-nt R1 run fits v4.
3. `select.py` misses these clean v5 studies:
   - "PIPseq V T2" / "PIP-Seq V T20" / "PIPseq-V T10" with no digit: PRJNA1445686, PRJNA1284893, PRJNA1477052.
   - "Pipseek V5": PRJNA1268916.
   - Protocols that use only the Illumina name, which its ENA QUERY never returns: PRJNA1269311, PRJNA1438756, PRJNA1481985, PRJNA1480533. PRJNA1377937 names PIPseq V only in a heading ("PIPseq V Cell Preparation").
4. Part of `select.py`'s QUERY, `"*Fluent*Bio*"`, returns only 4 runs today, while `"*Fluent*"` returns 12,641. The `*PIPseq*` part carries the draw.

## Reproduce

```sh
cd <this dir>
python3 -I -B scripts/ena_fetch.py ena          # writes ena/*.tsv (gzipped afterwards)
python3 -I -B scripts/geo_fetch.py geo
python3 -I -B scripts/geo_soft.py geo > geo/geo_hits.txt
python3 -I -B scripts/ena_projects.py geo/geo_hits.json runs.tsv ena/geo_projects.tsv.gz
python3 -I -B scripts/classify.py ena runs.tsv geo/geo_hits.json
python3 -I -B scripts/categorise.py runs.tsv ena geo/geo_hits.json > studies_table.md
```
