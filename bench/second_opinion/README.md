# Second opinions on the disputed calls

A HIGH call is **disputed** when the run's label names another kit and only the reads' own evidence
supports the call (`bench/coverage/results/adjudication.tsv`, rule of 8 Oct 2026). At coverage
a62c4ed, 999 of the 8,928 HIGH calls on held-out runs were disputed. Before any tool ran, three
methods other than tchem were set to read the same files and say which side they take
(`PLAN.md`, fixed on 10 Oct 2026):
1. **Sanger's `infer_10x_chemistry.py`.** The Wellcome Sanger Institute's Cellular Genetics checker
   (cellgeni/nf-reprocessing-public-10x @ a8611621). It covers 10x only and matches barcode lists
   and read geometry.
2. **A strand referee for 3' vs 5'.** It checks whether the cDNA read is sense or antisense on the
   spliced transcripts of 40 genes. The genes are disjoint from those tchem uses, and any k-mer
   shared with tchem's `data/reference/housekeeping_human.fa` is dropped. It uses GENCODE v44 on
   GRCh38 for human runs and GENCODE M35 on GRCm39 for mouse runs.
3. **Cell Ranger 10.1.0 `count --chemistry auto`.** It runs on a seeded sample of 40 human disputes,
   against refdata-gex-GRCh38-2024-A, one run at a time.

**None of this changes a level.** All three read the reads. Under the independent-label rule, they
corroborate a call and do not stand in for an independent label. The record keeps "disputed" as its own column.
This bench shows what other methods say about those calls.

## Results

Every number below is recounted from `results/` by `scripts/summarize.py results`.

| Method, on which disputes | Disputes | Sides with tchem | Sides with the label | No verdict |
|---|---:|---:|---:|---|
| Sanger, 10x label and 10x call | 914 | **851** | **0** | 31 cannot split (one list); 32 no call (21 variable read lengths, 6 read geometry, 5 Flex) |
| Sanger, 10x reads under another vendor's label (BD, Parse, MGI, Singleron) | 12 | **12** (finds a 10x list) | **0** | 0 |
| Sanger, another vendor's reads under a 10x label (Singleron, SeekGene, BD, MGI, inDrops, Parse) | 22 | **22** (finds no 10x list) | **0** | 0 |
| Strand referee, label and call differ in 3' vs 5' | 127 | **77** | **0** | 33 undetermined; 17 no referee (11 species not stated, 6 other species) |
| Strand referee, on the 31 Sanger cannot split | 31 | **21** | **0** | 6 undetermined; 4 no referee |
| Cell Ranger, seeded sample of human 10x disputes | 40 | **40** | **0** | 0 |

Cell Ranger's chemistries on the 40 were "Single Cell 3' v3" 24, "3' v2" 6, "5' PE" 5,
"5' R2-only" 3 and "3' v4" 2, each the layout tchem called. The sample includes 8 runs labelled
3' that tchem called 5', a run labelled Multiome that tchem called 3' v2, and a Flex-labelled run
tchem called 3' v3. The median wall time was 69 s per run (48 to 144 s).

The other 51 disputes have no label claimed, or neither side is a 10x kit:
- 30 Flex calls whose label is GEM-X Flex v2 (not claimed). Sanger has no Flex list.
- 6 3' v3 calls on runs whose label is not claimed (4 GEM-X Flex v2, 2 SeekGene MM). Sanger finds
  the v3 list on 4 and refuses 2 for variable read lengths.
- 8 SeekGene DD calls on MM-labelled runs. Sanger finds no 10x list.
- 7 within a vendor: BD bead version, and PIPseq v3 called V. Sanger finds no 10x list, which is
  consistent with both sides and judges neither.

## What the numbers mean

**Not one of the 999 disputes was settled for the label.** Wherever a method could tell the two
kits apart, it agreed with tchem.

**"Cannot split" is a limit of the method, not a dissent.** 10x 3' v2 and 5' v1/v2 use the same
barcode list (737K-august-2016). A tool that only matches lists can say no more than "3' v2 or
5' v1/v2". That is consistent with both tchem and the label, so it settles nothing. tchem goes
further and reads the TSO and the cDNA strand. For these runs, the strand referee is the method
that can decide.

**Sanger's no-calls are its own strictness, not disagreement.**
- It refuses a file whose cDNA reads were trimmed to varying lengths before upload, which is common
  in public archives.
- It refuses a file whose read lengths do not match the geometry it expects.
- It has no Flex barcode list, so it never calls Flex.

It named a different kit in none of them. tchem calls such files anyway, which is part of what it
is for. These runs are left with one witness fewer.

**"Undetermined" means too few reads fell on the referee's genes to cross its bar.** The bar is at
least 0.2% of reads and a 3-to-1 strand ratio. **"No referee"** means a species other than human or
mouse, or a run whose species is not stated.

**Cell Ranger is slow here by design.** It has no standalone chemistry detector. `--chemistry auto`
is a step inside `cellranger count`, which needs a transcriptome reference and aligns reads to it.
3' vs 5' is decided from where they land. That takes 11 GB of reference, about 32 GB of memory and
about a minute per 200,000-read sample, one sample at a time. It is efficient for its own users,
who have the reference anyway. tchem answers the same question without a reference, in about a
second.

**What this does not show.**
- The methods read the same reads as tchem, so they show that other methods agree. They do not show
  what the lab did.
- 50 of the 127 3' vs 5' disputes (33 undetermined, 17 with no referee), and the 32 Sanger could not
  call, are left with no second opinion.
- Free-text labels in public archives are wrong often enough that 11% of confident calls disagree
  with them. That limits how high a level can climb on public data. It does not limit the tool.

## Found on the way: damaged read heads

Sanger refused two files for malformed records. Both were clean at ENA and damaged in our copy.

A head is a byte-range prefix of a single-member gzip. It is decoded without its CRC, which sits at
the file's end, so bytes lost in transfer can decode into plausible text. A scan of every held-out
head (`scripts/scan_heads.py`) found 6 of 10,732 such files, listed in `results/damaged_heads.tsv`:
- tchem had refused two of them: "record 1 does not start with '@'" and "record 97139 has no '+'
  separator".
- The other four were damaged beyond where tchem reads.

The fetcher now checks every record and fetches again before failing (a62c4ed). The six were fetched
again and called fresh. The two refusals became right HIGH calls (10x 3' v3 and SeekGene DD 3'), and
the four others were called as before. The coverage table at a62c4ed carries this.

## Files

| File | What |
|---|---|
| `PLAN.md` | the scope, mappings and outcome rules, fixed before any tool ran (amended once before running, for species) |
| `scripts/second_opinions.py` | Sanger and the strand referee on every disputed HIGH call |
| `scripts/cellranger_sample.py`, `fill_chemistry.py` | the Cell Ranger sample; its chemistry is read from each run's `web_summary.html` |
| `scripts/summarize.py` | every number above, from `results/` |
| `scripts/scan_heads.py` | the FASTQ structure check over every head file |
| `results/second_opinions.tsv` | per disputed run: label, call, Sanger's id or error, strand shares and end |
| `results/cellranger.tsv` | per sampled run: label, call, Cell Ranger's chemistry, wall seconds |
| `results/species.tsv` | ENA `scientific_name` per disputed run |
| `results/damaged_heads.tsv` | the 6 damaged head files and their first bad record |

The referee FASTAs are built by `bench/head_to_head/scripts/build_nuclear_ref.py` from the gene lists
in `PLAN.md`. Their sha256 prefixes are 1c24964746d16cad (human, 426 transcripts) and
07b3d235d8f72bfd (mouse, 234 transcripts).
