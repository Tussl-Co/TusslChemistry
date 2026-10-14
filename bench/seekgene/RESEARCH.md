# SeekGene SeekOne: research and design runs

## What the design runs showed (9 Oct 2026, read once each, now design data)

| kit | design run (study) | Read 1 as read | layout |
|---|---|---|---|
| DD 3' (`DDV2`) | SRR27533148 (PRJNA1064036) | 17-base barcode, 12-base UMI, polyT from base 29 | `SeekGene SeekOne DD 3'` |
| DD 5' (`DD5V1`) | SRR23931570 (PRJNA947149) | the same barcode and UMI, TSO `TTTCTTATATGGG` at bases 29-41 | `SeekGene SeekOne DD 5'` |
| MM (`MM`, `MM-D`, `DDV1`) | SRR22784277 (PRJNA913187) | 1-4 base anchor, three 8-base barcodes split by `ACTGGTGA` and `GGTAGTGACA` (94% of reads), 8-base UMI, polyT | not claimed |

- The 17-base DD barcode does not fit tchem's 16-base list limit, so it is checked in two parts:
  bases 8-17 against SeekGene's list (1,536 10-mers, `src/lists/seekgene.rs`), bases 1-7 by
  concentration on the 960 7-mers. As a list, those 960 7-mers matched any read within one
  substitution: five unrelated oracle reports turned into SeekGene DD 3' LOW calls, so the list
  part was dropped.
- **MM is not claimed.** A constant amplicon sequence at base 34 sits one substitution from one
  of the 479 8-mers, so the list vouched for tag libraries at LOW. MM waits on a vouching check
  that hits spread over many barcodes.
- Not claimed either: DD-Q (scFAST-seq, FFPE) and the SimpliPlex sample-barcode read.

> Research by a metadata-only agent on 9 Oct 2026, kept as written. Its tables are here; its
> scripts, downloads and vendor documents stay out of the repository (they read local downloads,
> and some documents may not be redistributed). Paths it names under `dl/`, `raw/`, `docs/`,
> `src/` or `scripts/` are not committed.


## Agent research: SeekGene (SeekOne) kits: research for tussl-chemistry

Date: 9 Oct 2026. Sources were metadata and code only. No FASTQ or sequence data was downloaded or opened, and tchem was not run.

## Files

| File | Content |
|---|---|
| `chemistries.tsv` | One row per chemistry as coded. Segment positions are 1-based on the barcode read. |
| `lists.tsv` | Barcode and linker files: pinned raw URL, length, count, sha256 (as downloaded and after the tussl `Lines` fetch recipe), users, licence. |
| `runs.tsv` | 519 ENA runs plus 4 GEO-only series that have no ENA runs. Each row has its stated chemistry, file count, `first_public`, spot length (base_count/read_count), any kit version stated verbatim, and the protocol window. Each window was checked to contain the token behind its label. |
| `design_proposal.tsv` | One design run per chemistry, chosen from metadata only. |
| `scripts/` | The scripts that made the tables (`runs.py` holds the classifier). |
| `dl/` | Downloaded code, lists, docs and metadata. These are untrusted copies. |

## Primary sources (pinned)

- **seekgenebio/seeksoultools** at `a9a508ea75e9dc41d7821f47f7475ffbfc0ab74a` (default branch `v1.3.0`, HEAD on 10 Jun 2025). Licence: MIT, "Copyright (c) 2023 seekgenebio".
  - `chemistry.py` last changed in `3c94763c` (v1.3.0, 5 Jun 2025). The files are identical at HEAD.
  - The 2023 version is `eb5a0f6f5b9dd18f15ccee7165239f5edc7f0121`.
- **seekgene/SeekArcTools** at `a49e76b0f443b1fd9568381ee56806d19873f51f` (MIT). This is the Multiome ATAC + RNA kit.
- **seekgene/SeekSoulMethyl** at `f67b12cb6acdc20e34d7f861e600a480ebfd401e` (branch `nf_rna_methy`, MIT). This is the Multiome Methylation + RNA kit.
- **seekgenebio/seekonetools** at `edd515bf481ffa75de671a69498f800e878cd029` (2022, MIT). This is the old pipeline. Its P3CB copy is the BOM/CRLF version, so do not use it.
- **seekgene/seekarc_to_10x** has no licence. Do not cite it as a list source.
- **seekgene/SimpliPlex_demux** at `c8fb8d44f6c14f8b18a25465cc74104deab11cfc` (no licence, README only). It gives the SeekOne DD SimpliPlex R1 layout: CB 17, anchor 7, linker 17, sample barcode 10 (bases 42-51), UMI 10, polyT 30, then cDNA. SeekSoul Tools runs it as `--chemistry custom --structure B17X34U10`. Its M1-M16 whitelist is exactly `P3CBGB/sample.barcode.txt.gz`, in the same order.
- **Live docs** at `seeksoul.online/cloudplatform-doc/en/document/Software/3_SeekSoul_Tools/1_v1.3.0a/` (rna `3_Use_guide/1_rna/2_Analysis.html`, fast `3_Use_guide/2_fast/2_Analysis.html`, release `3_release.html`; fetched 9 Oct 2026). Quotes:
  - "DDV2 corresponds to 3' Transcriptome Kit; DD5V1 corresponds to 5'+ Immune Profiling Kit".
  - "DD-Q corresponds to scFAST-seq Kit".
  - The fast module serves "scFAST-seq Kits and SeekOne® DD Single Cell 3' RNA-seq Kit for FFPE Samples".
- **SeekSoulTools docs.** `seeksoul.seekgene.com` does not resolve now (DNS failure on 9 Oct 2026). The Wayback copies were used; the URLs are in `chemistries.tsv`. Key quotes:
  - v1.3.0 rna/run: "DDV2 corresponds to the SeekOne™ DD Single Cell 3' Transcriptome-seq Kit. DD5V1 corresponds to the SeekOne™ DD Single Cell 5' Transcriptome-seq Kit. MM corresponds to the SeekOne™ MM Single Cell Transcriptome Kit. MM-D corresponds to the SeekOne™ MM Large-well Single Cell Transcriptome-seq Kit."
  - v1.2.2 fast/run: "DD-Q corresponds to the SeekOne™ DD Single Cell Full-length RNA Sequence Transcriptome-seq Kit." The FFPE example uses `--chemistry DD-Q` with `--scoremin 0.2 --matchnmin 0.33`.
  - v1.3.0 rna/description: "In MM design, to increase the base balance of the linker portion during sequencing,1-4 bp shifted bases, called anchors, were added. The anchor determines the starting position of the barcode." It also says the pipeline looks for the anchor "within the first 7 bases of the Read1".
  - v1.2.2 vdj/description: "The Read1 structure of the SeekOne™ DD Single Cell Immune Profiling Kit is B17U12".
  - v1.3.0 release notes: "In the FAST and FFPE modules, the 17 base with low-quality bases are removed following the extraction of Reads1 UMI." This explains why DD-Q changed from `B17U12` (2023) to `B17U12X17`.
- **Notation.** `B` barcode, `L` linker (checked against a list), `U` UMI, `X` skipped. The structure is read on R1 (`process_barcode` cuts R1, and R2 is the cDNA). Every listed barcode segment uses the first list unless a second list is given (`barcode_wl_dict.get(B, barcode_wl_dict[0])`). Matching allows 1 mismatch (`match_type` 1); the DDVS sample barcode allows up to 3.

## Findings that matter for implementation

1. **MM, MM-D and DDV1 are one layout.** They have the same R1 structure (`B8X8B8X10B8U8`), the same 479-entry list at all three positions and the same shift. They differ only in adapter trimming. By the `mod.rs` rule ("one entry per processing-equivalent layout") they need a single entry.
2. **The anchor offset is outside `OFFSET_MAX`. This blocks implementation.**
   - In code, the barcode starts one base after the first `A` in `R1[0:7]`, so the 0-based offset is 1 to 7, and offset 0 is impossible.
   - The docs say the anchor is 1 to 4 bp, which gives an offset of 1 to 4.
   - Both ranges exceed tussl's `OFFSET_MAX = 3`. The anchor sequences themselves were not found.
   - The design run must measure the real offset distribution.
3. **The MM linkers are not checked by any shipped chemistry (`X`).**
   - `ACTGGTGA` (8) and `GGTAGTGACA` (10) come only from the unused `__SO01V3` entry. Treat them as inferred until the design reads confirm them.
   - The list folder name `SO01V3` matches the MM catalogue number "SO01V3.1" in ENA texts (PRJNA913187 and others).
4. **R1 tails are inferred from the adapter lists.** `adapter1` is trimmed from R1 after the structure is cut.
   - MM, DDV1, DDV2 and DDVS: `TTTTTTTTTTTT` at end 5, which is polyT.
   - DD5V1: there is no R1 tail entry. `adapter2` has `CCCATATAAGAAA`, the reverse complement of `TTTCTTATATGGG`, which is tussl's `TSO`. The SeekSoulMethyl README gives the same 13-mer for the 5' DD-MET5 library.
   - DD-Q has no polyT. Its R1 from base 47 is cDNA, because FAST aligns both reads.
5. **P3CB (17-mer, 1,474,560 entries) must be pinned to `a9a508e`.** The 2023 copy has 959 entries with an embedded U+FEFF, which tussl's fetch recipe would not strip. Byte-identical fallbacks: SeekArcTools and SeekSoulMethyl at the SHAs above.
   - The list is an exact product of 960 7-mers (positions 1-7) and 1,536 10-mers (positions 8-17). This is an observation only.
   - tussl's fetch script already gunzips (`fetch_whitelists.sh`, `normalise`), so `Recipe::Lines` works for the `.gz` files.
6. **Multiome uses its own read-1 tails** (SeekArcTools):
   - `DD_AG`: `B17U12X7`, with 7F = `TTGCTGT`.
   - `DD5_AG`: `B17U12X13`, with `TTTATTATATGGG`. This differs from tussl's TSO at one base, so do not merge them.
   - `DD_AA` (ATAC): `B17U12`, then 43 bases are cut (`r1.sequence[43:]`). 43 = 7F (7) + 17L `CGTCCGTCGTTGCTCGT` (17) + ME `AGATGTGTATAAGAGACAG` (19). The 17L and ME sequences come from the SeekSoulMethyl README.
7. **DDVS and SimpliPlex conflict.**
   - DDVS (`B17U12X24B10`) reads a sample barcode from the same 16-entry list at R1 bases 54-63.
   - The SimpliPlex README puts it at 42-51, with a 10-base UMI.
   - Which layout the shipped kit uses is not found, and no public run names either.
8. **Licences and terms.** seeksoultools, SeekArcTools and SeekSoulMethyl are MIT. No terms were found that forbid fetching the lists at run time.

## Public runs

- **Search.** ENA portal API, `result=read_run`, on `library_construction_protocol`:
  - `*seekgene*` returned 324 runs, `*seekone*` 489, `*seek*one*` 489, `*seek*soul*` 9, `*seek*gene*` 324, `*seek*` 1,367, and "seek one" / "seek-one" 0.
  - Of the `*seek*` hits, 857 were other tools (Curio Seeker, MuSeek, PIPseeker, BS-Seeker, ChIPseeker). That left 510 runs in 77 studies.
- **GEO.** E-utilities `db=gds`, GSE entries, terms seekgene, seekone, seeksoul and "seek one": 31 series. 24 are already in ENA through their BioProject. The other 7 were joined to ENA by BioProject: 3 have ENA runs (9 runs), and 4 have none. GEO text was read from the GSE SOFT metadata.
- **Classification.** Each run was labelled from its protocol text. Only the text near "SeekOne", "SeekGene" or "K0xxxx" was used, and a library kit had to be named. SeekMate dissociation kits, cell counters and "sent to SeekGene" are not counted as kit evidence. Catalogue numbers seen: K00202 (DD 3'), K00501 (DD 5'), K00601/K00701 (V(D)J), K00104 and SO01V3.1 (MM), K00801 (whole-sequence / full-length).
- **Totals.** 84 distinct studies: 80 with ENA runs and 4 GEO-only with none.
  - 73 studies have at least one paired run (exactly 2 FASTQ files).
  - In 68 studies every run is paired.
  - Of the 519 runs, 463 have 2 files, 6 have 1 file and 50 have no FASTQ in ENA.

| Stated chemistry | Studies | Studies with paired runs | Runs | Paired runs |
|---|---|---|---|---|
| DD 3' (DDV1 or DDV2; never stated) | 37 (36 + 1 GEO-only) | 32 | 236 | 197 |
| MM | 16 (14 + 2 GEO-only) | 14 | 87 | 86 |
| DD 5' (DD5V1) | 7 (6 + 1 GEO-only) | 4 | 43 | 32 |
| DD, end not stated | 5 | 5 | 30 | 30 |
| DD full-length (DD-Q; includes K00801 "whole transcriptome") | 5 | 4 | 12 | 10 |
| DD multiome ATAC + RNA (SeekArcTools) | 1 | 1 | 4 (2 ATAC, 2 GEX) | 4 |
| Mixed: DD 3' and 10x in one study | 2 | 2 | 18 | 15 |
| Mixed: MM and 10x in one study | 1 | 1 | 4 | 4 |
| Mixed: MM and DD 3' in one study | 1 | 1 | 6 | 6 |
| Mixed: 10x Chromium and "SeekOne platform" both named | 1 | 1 | 48 | 48 |
| Not SeekOne (10x kit; SeekGene named for dissociation or protocol) | 3 | 3 | 14 | 14 |
| SeekGene named, no kit stated | 1 | 1 | 4 | 4 |
| SeekSpace (spatial, out of scope) | 4 | 4 | 13 | 13 |
| DDVS, SimpliPlex, MM-D, DDV1 by name, DD-MET3/5 | 0 | 0 | 0 | 0 |

## Design proposal (rule and choice)

- **Rule.** Take the earliest study by ENA run `first_public`. Within that study, take the paired run with the lowest run accession. Break a tie on `first_public` by the lowest run accession.
- **Read lengths.** All chosen runs have a spot length of 300 (base_count/read_count) and an R1/R2 gz byte ratio between 0.90 and 1.04. That suggests 150 + 150, but metadata does not split the length by mate.

| Layout | Run | Study | first_public |
|---|---|---|---|
| MM layout (also MM-D, DDV1) | SRR22784277 | PRJNA913187 | 2023-02-27 |
| DD 3' (expected DDV2: "kit V2.1") | SRR27533148 | PRJNA1064036 | 2024-01-18 |
| DD 5' (DD5V1) | SRR23931570 | PRJNA947149 | 2023-03-27 |
| DD-Q full-length (tie with PRJNA1332757, broken by run accession) | SRR32280356 | PRJNA1175628 | 2025-09-28 |
| Multiome ATAC (DD_AA) | SRR38275770 | PRJNA1458007 | 2026-05-03 |
| Multiome GEX (DD_AG or DD5_AG) | SRR38275772 | PRJNA1458007 | 2026-05-03 |

These studies become design studies, and none of their runs may be held out.

**Multiome will have no held-out runs.** PRJNA1458007 is the only multiome study, and it supplies both multiome design runs.

## Unresolved

- **DDV1 vs DDV2.**
  - Kit versions are sometimes stated, verbatim in the `version_stated` column of runs.tsv: "kit V2.1" (PRJNA1064036), "v2 protocol" (PRJNA1412541), "Kit v1.1" (PRJNA1279856, DD 3'), "V1.1" (PRJNA1164785, full-length), "Protocol v1.2" (PRJNA1447531), "v1.3" (PRJEB106495, 5').
  - The same v1.x numbers appear on 3', 5' and full-length kits, so they look like kit or protocol revisions. No text uses the code names DDV1 or DDV2, and a kit version is not evidence of a chemistry.
  - Only reads can separate the three-part 479-list barcode from P3CB.
- **Products for DDV1 and DDVS:** not found (DDVS probably relates to SimpliPlex; see finding 7). **DDVS:** the 24 skipped bases and the role of the second barcode are not documented. No public runs name it.
- **MM anchor.** The anchor sequences and the true offset range (code: 1-7; docs: 1-4) are not found. The X-linker sequences are inferred.
- **DD-Q.** The content of the 17 skipped bases is not found.
- **DD5V1 TSO.** It is inferred from the adapter2 reverse complement, not stated.
- **Multiome GEX.** Whether a run is DD_AG (3') or DD5_AG (5') is not stated. The 7F sequence in the ATAC library is not stated (TTGCTGT is stated only for the GEX library).
- **Methylation libraries.** The read structure appears only in README figures and was not extracted.
- **K00801 "Single Cell Whole Trancriptome Kit"** (PRJNA1332757, PRJDB39486) is put under full-length because K00801 is the "Whole-Sequence Microbead Kit" in PRJNA1175628's text. This is an inference.
- **Mixed studies.** The per-run kit is not stated in the run metadata. Sample titles or GEO sample pages might resolve this.
