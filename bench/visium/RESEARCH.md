> Research notes written before this kit's catalogue was drawn (8 Oct 2026), kept as the record of
> sources and of which runs' reads were examined (those runs' studies are never held out: see
> results/examined.tsv). Scripts these notes name were scratch work and are not kept; the
> catalogue's rules are the ones in scripts/select.py, which say where they differ from these notes.

# Visium catalogue research (10x Genomics spatial gene expression)

Researched 8 Oct 2026. Sources: Cell Ranger 10.1.0 tarball (`$TUSSL_CHEMISTRY_DATA/cellranger/cellranger-10.1.0.tar.gz`, files under `lib/python/cellranger/`), the Space Ranger source on GitHub (https://github.com/10XGenomics/spaceranger, `main`), the ENA portal API, and the first 64 KB of ENA FASTQ files. The 10x support site (10xgenomics.com/support/...) returned HTTP 429 to every request, so nothing below is taken from the Space Ranger web docs. Claims that would need those docs are marked **to verify**.

## 1. Read structure (chemistry_defs.json, Cell Ranger 10.1.0; the same keys are in spaceranger `main` `lib/python/cellranger/chemistry_defs.json`)

| Chemistry | Spatial barcode | UMI | cDNA / probe | Strand | List |
|---|---|---|---|---|---|
| SPATIAL3Pv1 … v5 ("Visium V1…V5 Slide") | R1 offset 0, 16 bp | R1 offset 16, 12 bp (min 10) | R2, full length | + | visium-v1 … visium-v5 |
| SPATIAL-HD-v1 ("Visium HD v1", probe/RTL) | R1, two 14-bp "spot_segment" parts: bc1 at offset 11, bc2 at offset 25; `barcode_extraction: joint_bc1_bc2`, min_offset 8, max_offset 12 | R1 offset 0, 9 bp | R2, 50 bp (min 30) | − | slide `visium_hd_v1`, parts bc1/bc2 |
| SPATIAL-HD-v1-3P (HD 3', poly-A) | same as above | same | R2 full length | + | `visium_hd_v1` |
| SPATIAL-HD-v1-11mm / -3P-11mm | same | same | as above | −/+ | `visium_hd_11mm_v1` |

- Probe-based (FFPE, CytAssist) standard Visium has no separate chemistry: https://raw.githubusercontent.com/10XGenomics/spaceranger/main/lib/rust/cr_types/src/chemistry/mod.rs `is_rtl()` returns `None` for SPATIAL3Pv1–v5 with the comment "Custom and spatial chemistries may be either 3'GEX or RTL". So R1 is the same for poly-A and probe libraries; only R2 content (probe ligation products vs cDNA) differs. Only HD has a dedicated RTL chemistry (`SpatialHdV1Rtl` = SPATIAL-HD-v1).
- Slide serial prefix → chemistry (https://raw.githubusercontent.com/10XGenomics/spaceranger/main/lib/rust/cr_wrap/src/bin/spaceranger/slide_id.rs): V1→SPATIAL3Pv1 … V5→SPATIAL3Pv5, H1→HD (RTL or 3'). "6.5 mm (V1-, V2-, V3-, V4-prefix) slides", "11 mm (V5-prefix) slides"; `--unknown-slide visium-1`→v1 "original slides", `visium-2`→v4 "new 6.5mm slides", `visium-2-large`→v5 "new XL slides". That V4/V5 slides are the CytAssist slides is **to verify** against 10x docs (consistent with the reads below: every run naming CytAssist that I sampled matched visium-v4 or visium-v5, none visium-v1).
- Observed in reads (first 64 KB, 10 + 8 + 8 runs, exact match of R1[0:16]): poly-A and FFPE-v1 runs match visium-v1 at 67–96%; CytAssist runs match visium-v4 (61–98%) or visium-v5 (35–94%); HD runs match no Visium list. On HD runs, an exact bc1 followed immediately by bc2 was found in 30–38% of R1s, the matched bc1 entry starting at R1 offset 7–12 (mode 9). chemistry_defs gives bc1 offset 11 (14 bp) within an 8–12 window; the listed bc1 entries carry a leading G and are 15–16 nt, which may explain the shift. This is an observation from exact matching, not a reproduction of Space Ranger's extraction.

## 2. Barcode lists (counted from the 10.1.0 tarball)

| List | Barcodes | Notes |
|---|---|---|
| visium-v1.txt | 4,992 × 16 bp | identical to visium-v2.txt |
| visium-v3.txt | 4,992 | identical to visium-v4.txt; 0 shared with v1/v2 |
| visium-v5.txt | 14,336 | 0 shared with v1–v4 |
| visium_hd_v1.slide | protobuf: 3,350 bc1 + 3,350 bc2 seqs; pitch field 2.0 (µm); grid 3,350 × 3,350 | every bc1 = "G" + the bc2 at the same index; bc1 15–16 nt, bc2 14–15 nt |
| visium_hd_11mm_v1.slide | 5,524 + 5,524; grid 5,524 × 5,524 | same G+bc2 relation |

The protobuf field meanings (field 7/8/9 = bc1/bc2 lists, field 5 = 2.0 pitch and 3350×3350 dims) are my reading of an undocumented binary; **to verify**.

Overlap by exact match with single-cell lists in `$TUSSL_CHEMISTRY_DATA/whitelists/`: 10x_v3_permit (6,794,880): v1 15, v3/v4 8, v5 36; 10x_v4_3p_permit (7,372,800): v1 14, v3/v4 18, v5 43; 10x_v2 (737K): 0/1/0; 10x_v3_5p: 7/8/15; Flex v1: 0/0/6. Visium lists are effectively disjoint from every 10x single-cell list (< 0.3%); conversely Visium runs put only 0.1–0.5% of R1 prefixes in the 3' v3 list.

Availability and terms: shipped only inside Cell Ranger / Space Ranger installs (`lib/python/cellranger/barcodes/`). Not in the public spaceranger GitHub repo (tree of `main` lists code only) and not in github.com/Teichlab/scg_lib_structs (https://api.github.com/repos/Teichlab/scg_lib_structs/git/trees/master?recursive=1: 349 entries, untruncated, no Visium file or page). The spaceranger LICENSE: use "solely for your internal use", "in connection with a 10x Genomics Product", and "You agree not to redistribute or sublicense the Software" — so the lists should be read from a user's install, not vendored (https://github.com/10XGenomics/spaceranger/blob/main/LICENSE).

## 3. ENA

Query `library_construction_protocol="*visium*"` (case-insensitive; same count with `*VISIUM*`): **14,035 runs, 919 studies**; fastq_ftp files: 2 → 8,807 runs (599 studies), 1 → 4,129, 0 → 1,095, 3–4 → 4. library_strategy: OTHER 9,391, RNA-Seq 4,373, ATAC-seq 165, others small. `*cytassist*` adds 175 runs not in the Visium set; `*space ranger*`/`*spaceranger*` adds 33; `*spatial gene expression*` returns 0 (wildcard with spaces), `*spatial*gene*expression*` 370 (count endpoint).

Read lengths (modal R1/R2, first 64 KB via `Range: bytes=0-65535`):
- visium-v1 list (poly-A FF and FFPE-v1): R1 28 in 9/9; R2 120 (SRR11616463, SRR11832355, SRR12185472, SRR12685794), 90 (SRR22863870, SRR19777667, SRR23499457), 50 for FFPE (SRR23770994, ERR10030051).
- CytAssist (v4/v5 lists): 28/151 SRR39596215, 28/50 SRR38762024, 28/90 SRR32084295, SRR35144583; 100/100 SRR34449864, 101/101 SRR33423943, 151/151 SRR34190588, 150/150 SRR37350534.
- HD: 43/50 in SRR31716668, SRR31118556, SRR32360183, SRR33823579; 48/50 SRR31520430; 151/151 SRR30635153, 75/75 SRR31590030, 150/150 SRR32655959 (all with bc1+bc2 hits).
- Mixed study: SRR12103471 (PRJNA642635, protocol opens with a Chromium single-cell prep) 26/166, no Visium hits — a Chromium 3' v2 run filed in a Visium study.

Query URL: https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=library_construction_protocol%3D%22*visium*%22&fields=run_accession,study_accession,first_public,fastq_ftp,library_construction_protocol,library_strategy&limit=0

## 4. Proposal

- **Query**: `library_construction_protocol="*visium*" OR library_construction_protocol="*cytassist*"`.
- **Classifier** (on the protocol of the study's first run, as `draw()` does; script `scratchpad/scripts/visium_classify2.py`):
  1. require `/visium|cytassist/i`; refuse `names_other_vendor(p, "10x")`;
  2. refuse 10x single-cell kit names `/chromium|gem-x|single[- ]cell (3|5|gene expression|immune)|next ?gem|multiome|\bflex\b/i` (studies that mix Chromium libraries, see SRR12103471). A broader `single[- ]cell|scrna|snrna` rule wrongly refused PRJNA1198312, a verified HD run (43/50, "Illumina Single cell RNA seq preparation"); the kit-name rule keeps it;
  3. refuse protocols citing both the FF guide CG000239 and an FFPE guide (CG000407; CG000495 also included, its being a CytAssist guide is **to verify**), e.g. PRJNA915021 ("CG000239 … or … CG000407");
  4. `/visium\s*hd|visium high[- ]definition/i` → `visium_hd`; else `/cytassist|\bffpe\b|formalin|paraffin|probe/i` → `visium_probe`; else `visium_v1_ff`. Optional: in `visium_probe`, refuse protocols that also say fresh/snap-frozen or OCT without "CytAssist" (23 of 184 studies name frozen/OCT; PRJNA851592 says "freshly frozen in OCT blocks or embedded in paraffin", and its run is 28/90, likely poly-A).
- **Shape**: standard Visium `r1 >= 28 and r2 >= 40`; HD `r1 >= 43 and r2 >= 40`. Lengths cannot catch a lost barcode read when both reads are 100–151 bp (several valid runs are, barcode intact).
- **Design runs** (earliest study, first run by accession): `visium_v1_ff` PRJNA628708 / SRR11616463 (2020-04-30; 28/120; 92% visium-v1 hits); `visium_probe` PRJNA851592 / SRR19777667 (2023-02-03) under rules 1–4, but it is FF-or-FFPE; with the optional frozen refusal it is PRJNA894440 / SRR23770994 (2023-03-17; FFPE, 28/50, 96% visium-v1 hits); `visium_hd` PRJNA1160022 / SRR30635153 (2024-09-24; 151/151; HD hits 30%; protocol cites CG000685).
- **Pool** (union query → 14,210 runs, 924 studies; 2-file runs; cap 2/study; before shape): `visium_v1_ff` 298 studies / 564 runs (563 held out); `visium_probe` 184 / 341 (340); `visium_hd` 54 / 95 (94). Refused: 62 studies for 10x single-cell kit names, 3 for the FF+FFPE guide pair, 2 for other vendors (incl. PRJNA975472, home of the SRR24710534 Parse+Visium negative control). Alternative split by list instead of by assay: `visium_v1slide` (no CytAssist, non-HD) 404 / 761, `visium_cytassist` 78 / 144, `visium_hd` 54 / 95. All 26 sampled runs would pass the shape rule, so final numbers should be close; **to verify** with the full draw.

## 5. Risks

- v1 vs 3' v3: R1 layout is identical (16 bp barcode + 12 bp UMI, 28 bp R1, + strand); only the list differs. Lists are disjoint (< 0.3% overlap), so list-hit rate separates them cleanly, but the tool must load visium lists (not redistributable).
- v1==v2 and v3==v4 lists are byte-identical sets: slides V1/V2 and V3/V4 cannot be told apart from reads; report "Visium 6.5 mm (v1/v2 list)" etc.
- The list is set by the slide, not by FF vs FFPE: FFPE on V1 slides (SRR23770994, ERR10030051) hits visium-v1; CytAssist fresh-frozen (SRR37350534) hits v4. So `visium_probe` holds two lists. Probe vs poly-A cannot be told from R1; R2 = 50 bp is a hint for probe (FFPE) but CytAssist runs were also seen at 90/151 bp. Needs probe-set matching on R2 if it matters.
- HD needs a new layout model: 9-bp UMI, two variable-position 14–16 nt segments from a 3,350 (or 5,524) list with an 8–12 bp start window; a single fixed 16-bp lookup will call it "no known layout". Exact-match hit rates are only ~30–38%, so detection needs tolerance.
- Mixed studies: Visium studies often include Chromium 3'/5'/ATAC runs under one protocol text; the single-cell-word refusal removes 74 studies but residual mixing (protocol mentioning only Visium) remains possible.
- 1,095 runs have no fastq_ftp and 4,129 have one file (often BAM-derived or R2-only); excluded by the 2-file rule.

## Runs examined during research

Each run below had the R1 content of its first 64 KB examined, not only its read lengths. All of them must be excluded from the held-out pool; at most one per kit may serve as that kit's design run. The examination was: R1[0:16] exact match against visium-v1, v3/v4 and v5 and against 10x_v3_permit (all runs); for HD runs also the bc1+bc2 joint match against visium_hd_v1.slide. The proposed design runs SRR11616463, SRR19777667 / SRR23770994 and SRR30635153 are among them.

SRR11616463	v1 lists + 3'v3 permit
SRR11832355	v1 lists + 3'v3 permit
SRR12103471	v1 lists + 3'v3 permit
SRR12185472	v1 lists + 3'v3 permit
SRR12685794	v1 lists + 3'v3 permit
SRR22863870	v1 lists + 3'v3 permit
SRR19777667	v1 lists + 3'v3 permit
SRR23770994	v1 lists + 3'v3 permit
ERR10030051	v1 lists + 3'v3 permit
SRR23499457	v1 lists + 3'v3 permit
SRR30635153	Visium lists + 3'v3 permit + HD bc1/bc2
SRR31590030	Visium lists + 3'v3 permit + HD bc1/bc2
SRR31520430	Visium lists + 3'v3 permit + HD bc1/bc2
SRR31716668	Visium lists + 3'v3 permit + HD bc1/bc2
SRR31118556	Visium lists + 3'v3 permit + HD bc1/bc2
SRR32360183	Visium lists + HD bc1/bc2
SRR33823579	Visium lists + HD bc1/bc2
SRR32655959	Visium lists + HD bc1/bc2
SRR39596215	Visium lists
SRR33423943	Visium lists
SRR34449864	Visium lists
SRR38762024	Visium lists
SRR34190588	Visium lists
SRR37350534	Visium lists
SRR32084295	Visium lists
SRR35144583	Visium lists
