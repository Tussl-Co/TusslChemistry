> Research notes written before this kit's catalogue was drawn (8 Oct 2026), kept as the record of
> sources and of which runs' reads were examined (those runs' studies are never held out: see
> results/examined.tsv). Scripts these notes name were scratch work and are not kept; the
> catalogue's rules are the ones in scripts/select.py, which say where they differ from these notes.

# ICELL8 / Shasta (Takara Bio, formerly WaferGen): held-out catalogue notes

Researched 2026-10-08. "Observed" means it was measured on the first 64 KB (or 256 KB for the overlap check) of ENA FASTQs. Unverified claims are marked **to verify**.

## 1. Read structure

**3' DE (WaferGen ICELL8 3' DE, Takara ICELL8 / ICELL8 cx 3' DE)**
- Chip: 5,184 nanowells in a 72 x 72 grid. Each well has a preprinted oligo-(dT30), an 11 bp well barcode and a 10 bp UMI. Sequencing in the paper was R1 25 bp (barcode + UMI) and R2 50 bp (cDNA). Source: Goldstein et al. 2017, https://europepmc.org/article/PMC/PMC5501953 (full text: https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5501953/fullTextXML).
- Takara manuals: "A total of 5,184 unique nanowell barcodes … each printed once into the chip." The oligo(dT) primers "also include random nucleotide sequences of fixed length". The manuals give no R1 layout and no cycle numbers. Sources: ICELL8 cx 3' DE manual (022720), https://www.takarabio.com/documents/User%20Manual/ICELL8%20cx%203%27%20DE%20User%20Manual/ICELL8%20cx%203%27%20DE%20User%20Manual_022720.pdf; https://www.takarabio.com/resourcedocument/x98945 (640164-6); https://www.takarabio.com/resourcedocument/x102348 (older 640005 "UMI" kit, which has a separate "anneal UMIs and nanowell barcodes" step).
- Third-party pipeline (auto_process_ngs): the barcode is R1[1-11] and barcode+UMI is R1[1-21]. Source: https://raw.githubusercontent.com/fls-bioinformatics-core/auto_process_ngs/devel/tests/icell8/icell8_barcodes.sh. The test R1 reads are 26 bp.
- One submitter wrote "positions 1-11 … 11bp cell barcode, positions 12-22 … 10bp UMI" (PRJNA489248, ENA protocol text).
- **Observed: two R1 variants**, told apart by where the polyT starts. The per-position T fraction changes from about 0.3 to 0.98 or more at that point:
  - **v1 (11+10)**: barcode 1-11, UMI 12-21, polyT from base 22. Seen in SRR7782797, SRR7841080, SRR10668933, SRR11050949, SRR17031723 and SRR6860528, all from WaferGen-era studies.
  - **v2 (11+14)**: barcode 1-11, 14 random bases 12-25, polyT from base 26. Seen in DRR195272, DRR249480, SRR21046957, SRR15425309 and SRR25533399 (R1 of 76-150 bp). The 25-28 bp runs SRR13296907, ERR4192134, SRR17334665, SRR17327624, SRR8924646, SRR12476705 and SRR10914862 have random bases through position 25 with no polyT at 22, so they are consistent with v2. In 26 bp NextSeq runs, base 26 is mostly G (dark cycle) rather than T. Takara does not state the 14 nt length in any document found (**to verify** against a CogentAP/mappa config).
  - Some R1s are cut to 21 bp (SRR13094236, SRR21595032). There the variant cannot be told.
- R2 is cDNA (Nextera XT fragments). A noticeable fraction start with `TCCCTACACGACGC` (observed; origin **to verify**).
- R1 lengths observed: 21, 25, 26, 28, 76, 100, 101, 150. R2 lengths: 50, 100, 101, 108, 130/131, 150.

**Full-length (SMART-Seq ICELL8, ICELL8 cx, SMART-Seq Pro) and Shasta Total RNA**
- In these kits the cell is identified only by the i5 x i7 pair: 72 forward x 72 reverse 8 nt indexes, sequenced as 2 x 8 nt index reads. R1 and R2 are tagmented cDNA with flexible lengths (2 x 75 to 2 x 250). Sources: https://takarabio.com/resourcedocument/x101924 §Final Sequencing Library Structure; https://takarabio.com/resourcedocument/x337712 (SMART-Seq Pro uses 72 forward + 144 reverse in sets A/B).
- Shasta Total RNA-Seq (on ICELL8 cx or the Shasta system) also uses 72 i5 x 72 i7 indexes, plus an 8 nt TSO sample barcode (96 options) that "is incorporated into Read2". Single-end sequencing is not supported. Sources: https://takarabio.com/resourcedocument/x348320 and https://www.takarabio.com/resourcedocument/x351230 (Appendix C/D).
- ICELL8 ATAC and CUT&Tag also put the barcode in the I1/I2 index reads (https://raw.githubusercontent.com/fls-bioinformatics-core/auto_process_ngs/devel/tests/icell8_atac/README).
- **Consequence**: the ENA R1/R2 FASTQs of full-length, Shasta, ATAC and CUT&Tag runs contain no cell barcode. These kits are out of scope until the tool supports index-read barcodes.

## 2. Barcode lists
- 3' DE: 5,184 well barcodes of 11 bp, each printed once per chip (Goldstein 2017; Takara manuals above). The list ships with the instrument software as `3'DE.xml` (cx) or `5184 BarcodeSequences 21Aug2015.gal` (older; x102348 p.~26). **No public download was found** (**to verify**; ask Takara). A per-run "well list" (TSV with a Barcode column) is the documented format (auto_process_ngs test data).
- Observed: the same barcode set is used across vendors' eras and both UMI variants. In 15 runs from 15 studies, 80-92% of each run's barcodes with 3 or more reads also appear in other studies, and their union is 4,636 distinct 11-mers. A 5,184 list could be mined from public runs, but whether that list is complete and correct is **to verify**.
- Full-length and Shasta: 72 + 72 (Pro: + 72) 8 nt index sequences are printed in the manual appendices (x101924 App. A, x337712 App. A, x348320 App. D). These are index reads only.

## 3. ENA
- Query `library_construction_protocol="*icell8*"` (the wildcard is case-insensitive): **20,942 runs, 71 studies**. File counts: 14,646 runs with 1 file, 6,242 with 2, 53 with 0 and 1 with 3.
- `"*wafergen*"`: 66,579 runs in 207 studies. These are mostly bulk PrepX/Apollo runs. After removing ICELL8 runs, none mention single cells, nanowells or MSND. `"*shasta*"`: **0 runs**. `"*icell 8*"`, `"*smartchip*"`, `"*5,184*"`: 0 runs beyond the ICELL8 set.
- **Demultiplexing per cell is common.** 12 studies have 460-6,539 runs each with about 30k-480k reads per run, mostly single-file (PRJNA552315 has 6,539 runs; PRJEB29049 1,015; PRJNA491097 1,417). Two per-cell studies kept R1: PRJNA438958 and PRJNA587639 each have one barcode per run (observed: top barcode = 100% of reads).
- Pooled 3' DE runs (1-18 per study, 2 files, 10^7-10^9 reads) keep R1 with barcode+UMI. Their read lengths are listed in §1.
- Read lengths for 10+ runs (R1/R2, observed): SRR7782797 26/101; SRR13094236 21/130; SRR13296907 26/131; ERR4192134 28/100; DRR195272 100/100; DRR249480 100/100; SRR17334665 25/150; SRR17327624 26/108; SRR21046957 150/150; SRR21595032 21/100. Full-length runs: SRR10324769 76/76 (per-cell); SRR17010753 150/150; SRR24548571 101/101; SRR26815712 76/76; SRR31174141 101/101. Non-ICELL8: SRR16101452 28/91 is 10x.

## 4. Proposal
- **Query**: `library_construction_protocol="*icell8*"` (a superset; WaferGen and Shasta add nothing usable).
- **Classifier** (all cases need `\bi-?cell ?8\b|wafergen` in the text):
  - Refuse if any of these match: ATAC, CUT&Tag, small RNA, PSCSR, TCR, 5-prime. Refuse if `names_other_vendor` matches. This refusal already drops PRJNA767076 (10x) and PRJNA1014780. Add an `"icell8"` key to `VENDORS` so other kits refuse ICELL8.
  - `icell8_3de`: matches `3\s*['′’]?\s*-?\s*DE\b|\b30 DE\b|3['′’]\s*(end|differential)|differential expression|3['′’]? transcriptome` and no full-length words.
  - `icell8_fl`: matches `full[- ]length|SMART-?Seq (Pro )?ICELL8|SMART-?Seq Pro|Shasta Total|i5 … i7`. Do not match a bare "SMART-seq": PRJNA880899 says "SMART-seq 3' DE".
- **Shape**: icell8_3de takes `21 <= r1 <= 30 and r2 >= 50`. Optionally also `r1 >= 76`, since three of the 10 named studies have long R1s. icell8_fl takes `r1 >= 50 and r2 >= 50`, but its runs have no barcode to detect. Use them only as negative ("no chemistry") controls.
- **Pool estimate (run of these rules on today's ENA)**:
  - icell8_3de: 10 studies with paired runs pass. With 2 runs per study that is 16 runs. The design run is SRR7782797 (PRJNA489248, 2018), which is the **v1 11+10** variant. That leaves **15 held-out runs in 9 studies**: v2 in 7 studies and 21 bp R1 in 2. With the short-R1-only shape it is 7 studies and about 11 runs.
  - An optional `icell8_unnamed` (ICELL8 named, no kit words, short-R1 shape) would add 8 studies and about 16 runs. By read evidence these are 3' DE (four v1, four v2), but the text does not name the kit. Two of these studies are per-cell with one barcode per run.
  - icell8_fl: 12 studies with about 23 paired runs held out, no in-read barcode.

## 5. Risks and open questions
- Two R1 layouts (11+10 and 11+14) under one kit name. A design built on the 11+10 run will mis-place the polyT for v2 runs. Model the kit as barcode 11 + UMI 10 and leave bases 22-25 unmodelled, or treat v2 as a separate layout. The 14 nt length is **to verify** from a Takara or CogentAP source.
- The barcode list is not public: either use "unlisted" (concentration) mode or mine a list from public data (**to verify**).
- 21 bp R1 runs have no polyT. 26 bp NextSeq runs show G at base 26.
- Protocol text is often generic ("ICELL8 platform … Nextera XT"). This has low recall for `icell8_3de` and needs the read shape to separate unnamed runs.
- Per-cell submissions with R1 kept give 1 barcode per run, which is degenerate for concentration-based detection.
- Full-length, Shasta, ATAC and CUT&Tag runs need index-read support. SRA usually drops I1/I2 from FASTQ. SRR20574028 headers carry a `_GTACCGT` suffix (meaning **to verify**).
- Pool sizes are small, about 9-17 studies.

## Runs examined during research (content, not just lengths; exclude from held-out)

- DRR195272: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- DRR249480: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- DRR249481: polyT-run position in both files (64 KB)
- DRR249482: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- ERR2797244: polyT-run position in both files (64 KB)
- ERR2797245: polyT-run position in both files (64 KB)
- ERR4192134: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- ERR4192135: polyT-run position in both files (64 KB)
- ERR4192137: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR10324769: polyT-run position in both files (64 KB)
- SRR10324770: polyT-run position in both files (64 KB)
- SRR10324776: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR10399350: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR10399351: polyT-run position in both files (64 KB)
- SRR10668933: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR10668934: polyT-run position in both files (64 KB)
- SRR10914862: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR10914863: polyT-run position in both files (64 KB)
- SRR11050949: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR11050950: polyT-run position in both files (64 KB)
- SRR12476705: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR12476706: polyT-run position in both files (64 KB)
- SRR13094236: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR13094237: polyT-run position in both files (64 KB)
- SRR13296907: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR13296908: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB)
- SRR15425309: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR15668721: polyT-run position in both files (64 KB)
- SRR15668722: polyT-run position in both files (64 KB)
- SRR16101452: polyT-run position in both files (64 KB)
- SRR16101453: polyT-run position in both files (64 KB)
- SRR16101459: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR17010753: polyT-run position in both files (64 KB)
- SRR17010754: polyT-run position in both files (64 KB)
- SRR17031723: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR17031724: polyT-run position in both files (64 KB)
- SRR17031729: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR17327624: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR17334665: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR20574028: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB)
- SRR20574029: polyT-run position in both files (64 KB)
- SRR21046957: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR21046958: polyT-run position in both files (64 KB)
- SRR21046959: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR21595032: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR21595033: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB)
- SRR22428154: polyT-run position in both files (64 KB)
- SRR22428155: polyT-run position in both files (64 KB)
- SRR22733270: polyT-run position in both files (64 KB)
- SRR22733275: polyT-run position in both files (64 KB)
- SRR23715205: polyT-run position in both files (64 KB)
- SRR23715206: polyT-run position in both files (64 KB)
- SRR24548571: polyT-run position in both files (64 KB)
- SRR24548572: polyT-run position in both files (64 KB)
- SRR24548574: sequences printed, polyT/base composition, top 11-mers (64 KB both files)
- SRR24985951: polyT-run position in both files (64 KB)
- SRR24985952: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB)
- SRR25533399: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR25533400: polyT-run position in both files (64 KB)
- SRR25994672: polyT-run position in both files (64 KB)
- SRR25994673: polyT-run position in both files (64 KB)
- SRR26271327: polyT-run position in both files (64 KB)
- SRR26271329: polyT-run position in both files (64 KB)
- SRR26815712: polyT-run position in both files (64 KB)
- SRR26815713: polyT-run position in both files (64 KB)
- SRR27082127: polyT-run position in both files (64 KB)
- SRR27082129: polyT-run position in both files (64 KB)
- SRR27441392: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB)
- SRR27441393: polyT-run position in both files (64 KB)
- SRR27995465: polyT-run position in both files (64 KB)
- SRR27995466: polyT-run position in both files (64 KB)
- SRR28101461: polyT-run position in both files (64 KB)
- SRR28101462: polyT-run position in both files (64 KB)
- SRR31174141: polyT-run position in both files (64 KB)
- SRR31174142: polyT-run position in both files (64 KB)
- SRR38072809: polyT-run position in both files (64 KB)
- SRR38072810: polyT-run position in both files (64 KB)
- SRR38072934: polyT-run position in both files (64 KB)
- SRR38072935: polyT-run position in both files (64 KB)
- SRR38335111: polyT-run position in both files (64 KB)
- SRR38335113: polyT-run position in both files (64 KB)
- SRR6860528: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR6860529: polyT-run position in both files (64 KB)
- SRR7782797: sequences printed, polyT/base composition, top 11-mers (64 KB both files); polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR7841080: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB); R1 11-mer barcodes cross-study overlap (256 KB)
- SRR7841081: polyT-run position in both files (64 KB)
- SRR8395253: polyT-run position in both files (64 KB)
- SRR8395254: polyT-run position in both files (64 KB)
- SRR8396163: polyT-run position in both files (64 KB)
- SRR8396164: polyT-run position in both files (64 KB)
- SRR8831274: polyT-run position in both files (64 KB)
- SRR8831275: polyT-run position in both files (64 KB)
- SRR8924646: polyT-run position in both files (64 KB); R1 per-position T fraction, 11-mer barcode counts (64 KB)
- SRR8924647: polyT-run position in both files (64 KB)
- SRR8924648: sequences printed, polyT/base composition, top 11-mers (64 KB both files)

(95 runs.) Note: the pool estimates in section 4 included some of these runs. After excluding every run listed above, the named icell8_3de studies still have these unexamined paired runs:
- PRJNA687537: SRR13296909
- PRJEB38660: ERR4192136, ERR4192138
- PRJDB10643: DRR249483
- PRJNA869450: 6 runs (SRR21046960 to SRR21046965)
- PRJNA880899: SRR21595034, SRR21595035, SRR21595036

PRJNA489248, PRJNA679842, PRJDB8870, PRJNA792549 and PRJNA792135 have no unexamined runs left.

The realistic clean 3' DE held-out pool is therefore 5 studies with up to 2 runs each, about 9 runs. All 5 are v2 (11+14) or have a 21 bp R1. The design run must come from an examined run, or from a vendor demo.
