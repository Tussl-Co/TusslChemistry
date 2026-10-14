# How the validation record is made

The README section [How well it works](../README.md#how-well-it-works) gives the record that every
match prints. It also gives recall and what each level means. This page gives the rest:

- what the numbers show, and what they do not show;
- each benchmark, as it was when we first ran it.

[bench/VALIDATION.md](../bench/VALIDATION.md) gives the procedure that every kit follows.
[bench/second_opinion](../bench/second_opinion/README.md) gives the second opinions on the
disputed calls. [glossary.md](glossary.md) defines each term.

Short commit hashes in `bench/` and in the docs name commits in the development history. That
history is not published. The public history starts at v0.1.0.

## What the numbers show, and what they do not show

- **Independent labels decide every level.** An independent label is one of these:
  - the run's own protocol text;
  - the lab's own record. For 10x, this is 805 runs from 501 studies whose ENA or GEO record names
    exactly one kit (`bench/head_to_head/LAB_LABELS.md` and `bench/tenx_lab`);
  - the kit that the vendor states;
  - a BD part number that BD's own documents tie to one bead type (`bench/bd/BEAD_PARTS.md`).
- **Read-derived labels do not count.** For the other 10x runs, the label is Arc Institute's
  barcode-list matching. A right call against it shows that two list matchers agree. It does not
  show that tchem found the kit that the lab recorded. We report these calls, but a level does
  not count them.
- **Most disputes come from the labels.** Free-text records often name the wrong kit version.
  - Where the two kits use different barcode lists, the reads match the called kit's list at 74%
    or more, usually 88 to 100%. They match the labelled kit's list only weakly.
  - Where the two kits share one list (10x 3' v2 and 5' v1/v2), the call depends on the tail or
    the cDNA strand.
  - Two PIPseq runs with the label v3 show the ACG spacer of PIPseq V.
  - Four SeekGene studies state the MM kit, but show the barcode read of the DD kit.
  - Of the runs whose protocol names a 10x 3' or 5' kit, 15 have Singleron or SeekGene reads.
    Their lists match at 86 to 99%.
  - Of the same runs, 16 have Flex, Multiome, BD, MGI, Parse or inDrops reads (88 to 99%).
  - `bench/coverage/results/adjudication.tsv` gives the evidence for each dispute.
  - Two of the labelling rules are not certain, and the adjudications record this. The rules read
    "GEM-X Flex" as Flex v2 (34 runs). They also catalogued a study that names two assays as one
    assay (14 runs).
- **The second opinions agree with tchem, and never with the label.** Three methods other than
  tchem read the files of the 999 disputed calls (`bench/second_opinion`). We wrote the plan
  before they ran.
  - Sanger's `infer_10x_chemistry.py` agreed with tchem's barcode list in 851 of the 914 disputes
    where the label and the call are both 10x kits. It could not separate 31, because 3' v2 and 5'
    v1/v2 share one list. It made no call on 32 (trimmed reads, an unexpected geometry, or Flex).
  - A strand referee agreed with tchem in 77 of the 127 disputes where the label and the call
    differ in 3' or 5'. It gave no result on the others, or it had no genes for their species.
  - Cell Ranger's own `--chemistry auto` named tchem's layout in all 40 disputes of a seeded
    random set of human runs.
  - None of the three agreed with the label in any dispute. They read the same reads, so they
    support a call but do not change a level.

## Each benchmark, as first run

Each paragraph below gives a benchmark as it was when we first ran it, with the numbers of that
day. Since then, one coverage pass (`bench/coverage`) includes every held-out run. The README's
table shows that record, and it gives the current numbers for each layout.

### Long reads

`bench/long_read` used 59 public long-read 10x runs. We catalogued each run with its stated kit,
and froze the catalogue before tchem ran. The catalogue had:

- all five 10x gene-expression kits;
- Oxford Nanopore runs from 2019 to 2025;
- PacBio runs: raw Kinnex arrays, segmented reads, MAS-Seq and subreads;
- targeted amplicons and nuclei.

51 runs had a stated kit. The first pass named 48 of them at HIGH, made no call on 3, and got
none wrong. We then fixed the three causes that it showed: routing on read length, shared
barcodes and junk reads. After the fixes, tchem named all 51 (50 HIGH, 1 MEDIUM). Because those 51
runs then became design runs, quote the held-out test below.

### Long reads, held-out test

A second catalogue had 25 runs from 21 studies (`catalogue_heldout.tsv`). It had no run, study or
project in common with the first catalogue. We froze the catalogue and the tool before the test,
and ran the test once.

- 16 runs with a stated kit: tchem named all 16 (15 HIGH, 1 MEDIUM), none wrong.
- 3 runs whose labels name a 10x 5' kit, but not its version. tchem called each one a 10x 5'
  layout at HIGH (two 5' v3, one 5' v1/v2). They have no exact label, so they do not count.
- 5 runs of other kinds (10x Multiome, Parse, Visium, PIPseq, bulk Iso-Seq): tchem did not call
  any of them a 10x gene-expression kit. This includes PIPseq, whose reads have the 10x adapter.
- 1 run is available only as an NCBI SRA archive. It must be converted to FASTQ
  (`fasterq-dump`) before tchem can read it.

### Singleron

We used 20 public runs from 20 studies, and verified the kit of each run from its own reads. Of
223 studies that named Singleron, 17 did not contain Singleron libraries.

- 16 runs were Singleron GEXSCOPE: V1 6, V2 6 and V3 4.
- One run of each version was a design run, so 13 runs were held-out runs.
- tchem named all 16 at HIGH with the correct version, and all 13 held-out runs.
- The other 4 runs were GEXSCOPE MicroBead. tchem reports them as Drop-seq at MEDIUM, because the
  two kits have the same structure.

The commit that added the layouts gives these counts (2d0140e). The record did not keep the names
of the three design runs. Thus the coverage record counts only runs that it can trace
to held-out rows. These are the 3 runs that it names as held-out, and Singleron's own
demonstration runs.

### 10x Multiome GEX

This is the gene-expression library of the 10x Multiome ATAC + Gene Expression kit. Its reads look
like 3' v3, and tchem tells the two apart by the Multiome barcode list.

- One vendor run was the design run. The other 10 runs were held-out runs: 5 more from 10x, and 5
  archive runs from other labs whose protocols name the kit.
- tchem named all 10 correctly at HIGH. No other input in the parity oracle changed its verdict.
- The coverage record added one archive run. Arc's pipeline recorded its list as 3' v3, but the
  run's own label says Multiome.
- This gave 11 held-out HIGH calls from 7 studies, and the level `checked`.

### Parse Evercode WT

`bench/parse` drew 147 public runs from 78 studies by metadata alone. The protocol of each run
names the version of the WT kit. We froze the catalogue before we read any run.

The version of a Parse kit does not set the read layout:

- early v3 kits have the same layout as v1 and v2;
- later v3 kits use shorter linkers;
- some runs with the label v2 have the v3 layout.

The first held-out pass did not know the v3 layout: 37 runs had it, and all were LOW. We made the
v3 layout from 3 of those runs. Their 3 studies (6 runs) then became design studies, so 138
held-out runs remained (commit 7f48a4f).

On the 138 held-out runs:

- 125 were called Parse, 101 of them at HIGH.
- 9 had no call: R2 was trimmed short, or extra bases came before the structure.
- 4 were called 10x or BD at HIGH. The titles and protocols of all 4 name a 10x or BD library
  beside Parse, so the draw had the wrong label for them.
- tchem called the v3 layout on 36 runs, 22 of them at HIGH. The GGGG in its first linker reads
  badly on two-colour sequencers, so tchem finds that linker missing in many runs.

We know a fix: judge only the reliable bases of the linker. We will test it on runs that we have
not seen, not on these runs.

### MGI DNBelab C4 / C-series HT

`bench/mgi` used public runs whose protocols name the kit. It also used a rule on read shape.
Archives often drop the barcode read, or file the oligo library of the kit beside the cDNA. Each kit had one design run. We ran the held-out runs once:

- 52 runs from 26 studies.
- 49 were called MGI, 33 of them at HIGH. 1 had no call.
- 2 runs, from one study, were called 10x 3' v3 at HIGH. The study's metadata names only the MGI
  kit, but its 28-base barcode reads match the 10x v3 list at 97%. Chance gives about 3%.
- 13 of the MGI calls are LOW, all for one reason. A newer design puts the UMI at the start of R2,
  and a fixed GGCAGCGTC after the barcodes on R1. The layout does not know this design yet.
- Every run with the layout that MGI documents was HIGH or MEDIUM.

### BD Rhapsody, 2024 onwards

`bench/bd` used 511 held-out runs from 278 studies. The protocols of these runs name Rhapsody and
no other vendor. We set the bead versions on BD's own demonstration runs, so none of these runs is
a design run.

- 480 were called BD: 246 original V1 beads, 121 Enhanced and 113 Enhanced V2/V3 (the WTA Next
  generation). 466 of them were HIGH: 238, 118 and 110.
- 24 had no call.
- 7 were called 10x. The 6 at HIGH match a 10x list at 96 to 98%. They are 10x libraries in
  studies that name BD. The one at MEDIUM matches the 10x 3' v3 list at 52%.

### Illumina Single Cell 3' / PIPseq

`bench/pipseq`: v3 and v4 have the same layout. Each version had one design run. We ran the
held-out runs once:

- v3: 11 of 15 called. v4: 29 of 44 called. All were at HIGH, and none was wrong.
- The 19 runs with no call have an R1 trimmed below the 51 structured bases. Or they have a fixed
  sequence after the UMI where polyT must be, as v5 has.
- tchem does not claim v5 (5 runs), and it called none of them.

### ScaleBio QuantumScale

tchem names this kit from R2 alone: the UMI, then ScaleBio's RT barcode. It was right at HIGH on
both of ScaleBio's public test sets, and the second set was held-out. No archive run names the
kit yet.

### 10x Flex v1

This covers Fixed RNA Profiling, Next GEM and GEM-X. Its R1 looks like 3' v3, and tchem tells Flex
apart by the Flex barcode list. That list is available only inside Cell Ranger.

- One vendor run was the design run.
- 49 runs were held-out runs, and we ran them once. They were the 3 other vendor Flex runs, and 46
  runs from 24 studies by other labs whose own metadata names the kit (`bench/flex`). We froze
  the catalogue before the pass.
- R1 had 28 to 151 bases.
- tchem named all 49 correctly at HIGH, none wrong. No other input in the parity oracle changed
  its verdict.
- The 49 runs come from 25 studies. By studies, the level was `checked`, with an error bound of
  12% at 95% confidence.
- By runs, the bound would be less than 7%. But runs in one study are not independent.

Without the Flex list, tchem does not guess Flex. It reports that it could not test for Flex.

### Vendor gold set

`bench/vendor_gold` catalogued 164 datasets that the vendors published, each with the kit that the
vendor states. We froze the catalogue before tchem ran. We could download 128 of them, because
Parse publishes no FASTQ.

tchem named every supported kit correctly: 105 of 105 at HIGH, none wrong.

| Kit | Right at HIGH |
|---|---|
| 10x 5' v3 (GEM-X) | 4 of 4 |
| 10x 3' v4 (GEM-X) | 8 of 8 |
| 10x 3' v3/v3.1 (standard, HT, LT, nuclei) | 33 of 33 |
| 10x 3' v2 | 10 of 10 |
| 10x 5' v1 to v2 (with HT) | 21 of 21 |
| BD Rhapsody (10 Enhanced, 6 original beads) | 16 of 16 |
| ScaleBio RNA | 5 of 5 |
| Singleron | 8 of 8 |

- tchem did not support Multiome and Flex at that time. It called none of their 10 datasets at
  HIGH.
- After we added 10x Multiome GEX, the same set gave 111 of 111 at HIGH, with the 6 Multiome
  datasets (`results_multiome.tsv`).
- After we added Flex v1, it gave 115 of 115 (`results_flex.tsv`).

tchem declined, and did not guess, on kits that it did not know yet:

- BD's newer bead kits (WTA Next and intracellular AbSeq, 5 datasets);
- ScaleBio QuantumScale (2 datasets);
- three old 5' runs whose downloaded cDNA file had no transcript signal.

We fixed the scorer's mapping of vendor wording before the run. An audit of the scores found
three errors in how the scorer read that wording, and we corrected them separately. We keep both
scorings.

### Chaos testing

`bench/chaos` applies 17 kinds of damage to a real run of each of the 25 layouts. This gives 425
inputs, and each input runs twice. The kinds of damage are:

- swapped files;
- a constant cDNA read that matches a barcode list;
- a barcode read that is reverse-complemented;
- a barcode read that is trimmed;
- two kits mixed;
- only 500 reads;
- Windows line endings;
- lower-case bases;
- 30% Ns;
- 3% substitutions;
- R2 a copy of R1;
- no barcode read;
- files of unequal length;
- FASTA instead of FASTQ;
- empty files;
- truncated gzip;
- quality lines of the wrong length.

Each input must give the right verdict, a lower confidence, or a typed refusal, and the same
result both times. It must never cause a crash or a wrong HIGH call. `scripts/gate.sh` runs these
cases on every commit, with 50 hostile-input probes (names, paths, symlinks, BGZF and BAM). We
also fuzz the readers and detection (`fuzz/README.md`).
