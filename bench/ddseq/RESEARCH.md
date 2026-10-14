> Research notes written before this kit's catalogue was drawn (8 Oct 2026), kept as the record of
> sources and of which runs' reads were examined (those runs' studies are never held out: see
> results/examined.tsv). Scripts these notes name were scratch work and are not kept; the
> catalogue's rules are the ones in scripts/select.py, which say where they differ from these notes.

# Bio-Rad ddSEQ / Illumina Bio-Rad SureCell WTA 3' (kit `ddseq_wta3`)

Researched 8 Oct 2026. Nothing in the repo was changed. "Measured" means I measured it on ENA reads (first 64 KB of each file only).

## 1. Read structure

The bead oligo, from Teichlab: `AAGCAGTGGTATCAACGCAGAGTAC [BC1 6] TAGCCATCGCATTGC [BC2 6] TACCTCTGAGCTGAA [BC3 6] ACG [UMI 8] GAC (dT)V`. Read 1 starts after the custom primer `GCCTGTCCGCGG AAGCAGTGGTATCAACGCAGAGTAC`. Teichlab gives Read 1 = 68 cycles and Read 2 = 75 cycles (cDNA). The 8-bp sample index is on i7 (N7xx). Source: https://teichlab.github.io/scg_lib_structs/methods_html/SureCell.html

- **Phase block.** The ddSeeker paper (Romagnoli et al. 2018) is in *BMC Genomics* 19:960, not *Bioinformatics*. It says each oligo starts with a "phase block (PB)", and that "separation of cDNA from the beads can occur at different nucleotides within the PB, thus making the position of the two linkers variable". It checks that L2 starts exactly 21 nt after L1. Sources: https://doi.org/10.1186/s12864-018-5249-x and https://www.ebi.ac.uk/europepmc/webservices/rest/PMC6304778/fullTextXML. The authors of PRJEB33083 put a UMI-tools pattern in their ENA protocol: `(?P<discard_1>.{0,5})(?P<cell_1>.{6})(TAGCCATCGCATTGC){s<=2}(.{6})(TACCTCTGAGCTGAA){s<=2}(.{6})ACG(?P<umi_1>.{8})GAC T{0,6}`.
- **Measured phase.** I pooled R1 reads from 19 runs (36,011 reads, 19 studies). L1 starts at positions 6, 7, 8, 9, 10 and 11 in 5793, 5639, 6173, 5417, 4749 and 4625 reads. Positions 5 and 12–15 hold 54 reads, and 3563 reads (9.9%) have no exact L1. So BC1 starts at an offset of **0 to 5**, close to uniform. L2−L1 = 21 in 28,817 reads; 20 in 431. Where L2−L1 = 21, ACG…GAC sits exactly around the 8-bp UMI in 26,401 of 28,817 reads (91.6%).
- **Layout at phase 0.** BC1 0–5, L1 6–20, BC2 21–26, L2 27–41, BC3 42–47, ACG 48–50, UMI 51–58, GAC 59–61, then polyT. At phase 5 a 68-base R1 ends with GAC at 64–66 and only **one T**.
- **R2** is cDNA, 74–76 bases (measured).
- The **kallisto/kb `SURECELL`** definition matches phase 0 only: bc (0,0,6), (0,21,27), (0,42,48); umi (0,51,59); cDNA on file 1. It has no phase handling. Source: https://raw.githubusercontent.com/pachterlab/kallisto/master/src/main.cpp (around line 1361).
- **Illumina's Reference Guide** (doc 1000000021452 v01, marked "Illumina Proprietary") says only that Read 1 uses a custom primer and that analysis is done in the SureCell RNA Single-Cell App on BaseSpace. It gives no cycle counts and no barcode layout. Source: https://emea.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/surecell/surecell-wta3-library-prep-reference-guide-1000000021452-01.pdf. The 68/75 cycle numbers come from a search snippet of the Illumina product page; the page I fetched did not show them, so they are **to verify** from Illumina. ENA protocol PRJNA547819 states "68bp:Read1 and 75bp:Read2". I did not find the BaseSpace App documentation (**to verify**).

## 2. Barcode list

- There are 96 six-base blocks, and the same list serves BC1, BC2 and BC3. That gives 96³ ≈ 885k combinations. The minimum pairwise Hamming distance is 3 (measured).
- ddSeeker has the list hard-coded as `_barcodes` in `code/ddSeeker.py`, licensed GPL-3.0. Source: https://github.com/cgplab/ddSeeker/blob/master/code/ddSeeker.py
- The umis project has it as `examples/SureCell/barcodes.txt`, licensed MIT. Source: https://github.com/vals/umis/tree/master/examples/SureCell. The two lists are the same set (checked). Teichlab points to the umis copy.
- Illumina does not publish the list in the reference guide.
- **Prefer the MIT umis copy.**
- **Measured fit to the list.** I took 8,364 reads with both linkers at L2−L1 = 21 and BC1 complete, from 5 runs. All three blocks matched the list exactly in 88.2% of them. Per block: BC1 92.7%, BC2 97.8%, BC3 96.2%. BC1 often has an N at R1 base 5 (cycle-6 artefact on NextSeq).

## 3. ENA

The query `library_construction_protocol="*ddSEQ*" OR library_construction_protocol="*SureCell*"` (case-insensitive in practice) returns **3,741 runs from 69 studies**. Of these, **1,811 runs from 47 studies** have exactly two FASTQ files.

Wildcards with a space or hyphen (`*Bio-Rad*`, `*Illumina Bio-Rad*`, `*single-cell isolator*`) returned 0 runs. Use single tokens.

Many SureCell RNA studies have **one file only**: PRJNA438777 (507 runs), PRJNA699546 (732), PRJNA574184, PRJNA545723 and others. The rule that requires two files drops them.

ATAC studies are 22 of the 69. All of them are paired: R1/R2 lengths are 118/40, 122/36 or 150/150.

Read lengths, first run per study (R1, R2):

| study | run | lengths | L1 found |
|---|---|---|---|
| PRJNA547819 | SRR9311852 | 68, 75 | yes, pos 6–11 |
| PRJNA422868 | SRR8547723 | 68, 75 | yes |
| PRJEB33083 | ERR3379610 | 68, 75 | yes |
| PRJNA600520 (nuclei) | SRR10870711 | 68, 75 | yes |
| PRJNA520915 | SRR8528464 | 68, 74 | yes |
| PRJNA643512 | SRR12125440 | 68, 76 | yes |
| PRJNA705943 | SRR13826415 | 76, 76 | yes |
| PRJNA1108856 | SRR28960999 | 69, 76 | yes |
| PRJNA1269949 ("Isolator 2") | SRR33776330 | 68, 75 | yes |
| PRJNA594939 | SRR10669441 | 76, 76 | **no, in either file** |
| PRJNA433889, PRJNA657844 | SRR6713500, SRR12467642 | 26, 74 / 26, 98 | no (10x-like R1) |
| PRJNA738946/7 | SRR14856566 / SRR14856486 | 52, 75 | no; carries `TATGCATGAC` (ATAC/new-kit bead linker) |
| PRJNA640290 | SRR12043376 | 68, 75 | yes, but protocol names Chromium |

PRJNA509252 and PRJNA738945 are ATAC runs with 122-base R1. That R1 holds the SureCell RNA linkers and then `TCGTCGGCAGCGTC`. The other ATAC runs (118/150) use the `TATGCATGAC…AGTCACTGAG` linkers. **SureCell ATAC is a different read structure and out of scope.**

Bio-Rad's later **ddSEQ Single-Cell 3' RNA-Seq Kit** is "an entirely new assay". It uses a TruSeq R1 of at least 54 cycles, 7-bp blocks with linkers `TATGCATGAC` and `CAGTCAGTCA`, and Deconvolution-Oligo libraries. It is also out of scope. Source: https://teichlab.github.io/scg_lib_structs/methods_html/ddSEQ_3pRNA.html

## 4. Proposal

```python
QUERY = 'library_construction_protocol="*ddSEQ*" OR library_construction_protocol="*SureCell*"'
KIT = re.compile(r"SureCell|ddSEQ", re.I)
ATAC = re.compile(r"atac", re.I)
NEWKIT = re.compile(r"ddSEQ Single[- ]Cell 3.{0,2} RNA[- ]Seq|Omnition|deconvolution oligo", re.I)
def classify(p):  # kit ddseq_wta3
    if ATAC.search(p) or NEWKIT.search(p) or catalogue.names_other_vendor(p, OWN): return None
    m = KIT.search(p); return ("ddseq_wta3", snippet(p, m)) if m else None
def shape(r1, r2): return 66 <= r1 <= 80 and r2 >= 50
```

- `catalogue.VENDORS` has no ddSEQ/Bio-Rad key yet, so `names_other_vendor` needs a new entry. One option is `"biorad": re.compile(r"ddseq|surecell", re.I)`. In the simulation I instead refused a protocol if any existing vendor matched.
- **Simulated draw.** I ran it read-only against the repo's own `catalogue.draw` via a scratch script.
  - Result: **20 studies, 39–40 runs**, with design run PRJNA547819/SRR9311852 (2019-07-16). That leaves **38–39 held-out runs from 19 studies**.
  - The count changed by one run between two identical runs. A transient range-fetch failure makes `read_length` return 0, which drops a run.
  - The classifier refuses: PRJNA574650 (names 10X), PRJNA640290 (Chromium), every ATAC study, and the 52-base runs (by shape).

## 5. Risks and open questions

1. **The phase offset is 0–5, and the tool supports 0–3.** At offsets 4–5 about 30% of reads would fail. Detection needs a maximum offset of 5, or ddseq-specific handling. A fixed-offset model (kb) is wrong for about 5/6 of reads.
2. **PolyT is not a usable anchor.** At phase 5 only 1 T fits in a 68-base R1. Anchor on L1/L2 plus ACG…GAC instead.
3. **Linker variants.** Indels are common enough that ddSeeker models them: L2−L1 = 20 in 1.5% of reads, and the paper reports 1.2–1.4% linker events. No alternative linker sequences are known.
4. **Contamination.** PRJNA594939 passes the metadata and length filters but has no linker in either file, so its barcode read is lost or it is mislabelled. Tightening R1 to 66–70 removes it, but also drops the genuine 76/76 study PRJNA705943. Keep it and expect a correct "unknown/refuse", or exclude it by an explicit rule fixed now.
5. PRJEB87411 and PRJEB87408 are near-duplicate studies from one lab. PRJNA1269949 says "Isolator 2", but its reads are WTA-structured.
6. **To verify:** the cycle counts in Illumina's own text, and the documentation of the SureCell RNA Single-Cell App.

## Runs examined during research

I printed the first R1 read (first 70 bases) of the **first run, by accession, of every paired study** from the ENA query. That is the 47 runs below. Where R1 was 80 bases or shorter, I also counted where L1 starts in every read of the first 64 KB of R1. Other runs were seen only as read lengths, through `catalogue.read_length` in the simulated draw.

Codes:
- **A** = in the 19-run pool for the phase/linker/ACG–GAC measurement
- **B** = in the 5-run barcode-list fit
- **S** = first R1 read printed, plus L1 position scan on the first 64 KB of R1 (when R1 ≤ 80)
- **P** = first R1 read printed only (R1 > 80)
- **X** = both files: first two reads printed and L1 counted

SRR9311852 (PRJNA547819, design run): S, A, B
SRR8547723 (PRJNA422868): S, A
ERR3379610 (PRJEB33083): S, A, B
SRR10870711 (PRJNA600520): S, A
SRR8528464 (PRJNA520915): S, A
SRR12125440 (PRJNA643512): S, A
SRR11201417 (PRJNA609328): S, A
SRR10197365 (PRJNA574650): S, A
SRR13394502 (PRJNA690953): S, A
SRR13808894 (PRJNA705511): S, A
SRR18788302 (PRJNA827611): S, A
SRR17499020 (PRJNA795624): S, A, B
SRR10168924 (PRJNA573844): S, A
SRR13826415 (PRJNA705943): S, A
ERR14764724 (PRJEB87411): S, A
SRR28960999 (PRJNA1108856): S, A, B
SRR33776330 (PRJNA1269949): S, A, B
SRR27297071 (PRJNA1055016): S, A
SRR39621274 (PRJNA1495244): S, A
SRR10669441 (PRJNA594939): S, X
SRR6713500 (PRJNA433889): S
SRR12043376 (PRJNA640290): S
SRR12467642 (PRJNA657844): S
SRR14856486 (PRJNA738947): S
SRR14856566 (PRJNA738946): S
ERR14764693 (PRJEB87408): S
SRR8310244 (PRJNA509251, ATAC): P
SRR8310494 (PRJNA509253, ATAC): P
SRR8310525 (PRJNA509252, ATAC): P
SRR8310637 (PRJNA509255, ATAC): P
SRR8310661 (PRJNA509254, ATAC): P
SRR8994087 (PRJNA540583, ATAC): P
SRR8994111 (PRJNA540584, ATAC): P
SRR9331546 (PRJNA549875, ATAC): P
SRR11074305 (PRJNA623130, ATAC): P
SRR11321529 (PRJNA613046, ATAC): P
SRR14333205 (PRJNA713904, ATAC): P
SRR14855152 (PRJNA738945, ATAC): P
SRR13783625 (PRJNA704903, ATAC): P
SRR13386918 (PRJNA690721, ATAC): P
SRR23908682 (PRJNA798749, ATAC): P
SRR12817011 (PRJNA668908, ATAC): P
SRR23860411 (PRJNA944632, ATAC): P
SRR25915891 (PRJNA1013101, ATAC): P
SRR22444736 (PRJNA906529, ATAC): P
SRR28886066 (PRJNA714108, ATAC): P
SRR29712896 (PRJNA1132135, ATAC): P

**Consequence.** All 20 studies in the proposed draw had their first run's reads examined. These include 19 held-out runs: every one in group A except SRR9311852 (the design run) and SRR10197365 (PRJNA574650, which the classifier refuses), plus SRR10669441 and ERR14764693. The second-run held-out rows (lengths only) were not examined.

Not from ENA: the vals/umis `examples/SureCell/K562_R1.fastq` (first two reads).
