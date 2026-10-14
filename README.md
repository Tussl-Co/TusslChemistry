# TusslChemistry

TusslChemistry tells you which single-cell kit made a set of FASTQ or BAM files. Its command is
`tchem`. It reads the first 100,000 reads and needs no genome index. It also tells you where the
barcode and UMI are, which barcode list to use, and which strand to count.

```
$ tchem --whitelists ~/whitelists --r1 S1_R1.fastq.gz --r2 S1_R2.fastq.gz
...
VERDICT  10x 3' v3/v3.1  (confidence HIGH)
         geometry 1{b[16]u[12]x:}2{r:}
```

Each single-cell pipeline (Cell Ranger, STARsolo, alevin-fry or simpleaf, kallisto|bustools) needs
four facts about the reads. These are the position of the cell barcode, the UMI and the cDNA, the
barcode list, and the strand. If one fact is wrong, the pipeline does not stop. It reports fewer
cells, fewer genes or the wrong genes. Public metadata often does not give these facts, or gives
them incorrectly.

We could check 295 experiments in Arc Institute's scBaseCount. For 21, it records a 3' kit, but the
reads are from a 5' kit. The ENA title shows this for 5 of them, and a strand test on nuclear genes
for the other 16.

The terms in this README have one meaning each. [docs/glossary.md](docs/glossary.md) defines them.

## Quick start

```sh
git clone https://github.com/Tussl-Co/TusslChemistry && cd TusslChemistry
cargo install --path .                    # installs tchem into ~/.cargo/bin
scripts/fetch_whitelists.sh ~/whitelists  # barcode lists (10x's are ~315 MB of text)
# add --cellranger ~/cellranger-10.1.0 to the line above to also copy Flex's list from Cell Ranger
tchem --whitelists ~/whitelists --pack    # optional: faster loads
tchem --whitelists ~/whitelists --r1 R1.fastq.gz --r2 R2.fastq.gz
tchem --whitelists ~/whitelists S1_L001_*.fastq.gz      # a lane's R1, R2, I1, I2 in any order
tchem --whitelists ~/whitelists nanopore_reads.bam       # long reads: FASTQ or BAM
tchem --whitelists ~/whitelists SRR19374511              # a public run, read directly from ENA
```

Each GitHub release also has prebuilt archives with their sha256. The Linux builds (x86_64 and
aarch64) are static. The macOS builds are for Apple silicon and Intel. CI builds the Intel one,
but does not run it. An archive holds `bin/tchem` and the fetcher.

You can also set the barcode list directory with `TUSSL_CHEMISTRY_WHITELISTS`. Other options are
`--reads N` (default 100,000, maximum 1,000,000), `--all` (the evidence for every candidate) and
`--help`.

## What tchem reports

For each set of files, tchem reports:

- **The kit and its version.** tchem knows 25 layouts from 10x, BD Rhapsody, Parse Evercode,
  SPLiT-seq, Singleron, MGI, Illumina (PIPseq), ScaleBio, SeekGene, inDrops and Drop-seq. The
  report gives a confidence and the evidence for it.
- **The geometry** for simpleaf, piscem and alevin-fry, and the barcode list or lists to use.
- **The barcode read and the cDNA read.** tchem finds them when the files are swapped.
- **The strand of the cDNA read**: forward for 3' kits and reverse for 5' kits. tchem measures
  it without a reference genome.
- **The library type**: gene expression, immune receptor (VDJ), CRISPR guide or a tag library.
  These library types can share beads, but they need different pipelines.
- **The instrument**, from the read names or a BAM read group. For a public run, the report also
  shows the instrument that ENA records. When the read names do not identify the instrument, the
  report says "unknown".
- **Long reads.** For 10x cDNA on Oxford Nanopore or PacBio, tchem gives the kit. It also gives
  the `--kit` value that Oxford Nanopore's `wf-single-cell` needs, for example `3prime:v3`.
- **Index reads.** Give the files of a lane in any order. tchem finds the index reads (I1 and I2)
  and sets them aside.
- **The cause of a no match.** For example, the barcode read is missing, a common fault in ENA's
  two-file exports. Or the reads have a fixed structure that no known layout has.

A run takes about 1.3 seconds of CPU time for short reads (median of 72 public runs). It takes
about 0.1 seconds for long reads.

For long reads, the verdict gives the setting for the long-read pipeline:

```
VERDICT  10x 5' v3 (GEM-X), long reads  (confidence HIGH)
         wf-single-cell --kit 5prime:v3
```

### The JSON report

`--format json` prints one object for a pipeline to read. This is the report for a 10x 3' v3 run
from the gate's regression inputs:

```json
{
  "tool": "TusslChemistry",
  "version": "0.1.0",
  "schema_version": "1.1",
  "verdict": "match",
  "confidence": "high",
  "layout": {
    "name": "10x 3' v3/v3.1",
    "geometry": "1{b[16]u[12]x:}2{r:}",
    "barcode_lists": [
      "whitelists/10x_v3_permit.txt"
    ],
    "note": "",
    "long_read_kit": "3prime:v3",
    "validation": {
      "applies_to": "high",
      "level": "error_under_3pct",
      "held_out_high_calls": 4348,
      "studies": 2282,
      "family_only_calls": 0,
      "disputed": 637,
      "wrong": 0,
      "error_bound_95": 0.0013146362839614374
    },
    "family_validation": null
  },
  "candidates": [],
  "files": {
    "barcode_read": "R1.fastq.gz",
    "cdna_read": "R2.fastq.gz"
  },
  "cdna_strand": "forward",
  "library": {
    "kind": "gene expression",
    "gene_expression": true,
    "constant_region_share": 0.00212,
    "housekeeping_share": 0.14774,
    "scaffold_share": 0.00048,
    "tag_start_cover": 0.09316
  },
  "swapped": false,
  "strand_evidence": {
    "poly_a": 0.02551,
    "poly_t": 0.00329,
    "chrm_sense": 0.0581,
    "chrm_antisense": 0.00063,
    "nuclear_sense": 0.11542,
    "nuclear_antisense": 0.02696
  },
  "contradictions": [],
  "notes": [],
  "hint": null,
  "not_checked": [],
  "reads_sampled": 100000,
  "read_kind": "short",
  "instrument": {
    "platform": "illumina",
    "model": "NovaSeq 6000",
    "id": "A00228",
    "from": "read_name",
    "archive_model": null
  },
  "long_read_evidence": null,
  "fingerprint": null
}
```

`verdict` is `match`, `ambiguous` or `no_match`. An ambiguous verdict lists its `candidates`. A
no match gives a `hint`.

The object is a versioned contract. A pipeline must read `schema_version` and ignore keys that it
does not know. A new key increases the minor version. A renamed or removed key increases the major
version. [docs/json.md](docs/json.md) describes each
key, and [schema/report.v1.json](schema/report.v1.json) is the JSON Schema.

## Public runs by accession

You can give a run accession (SRR, ERR or DRR) or an experiment accession (SRX, ERX or DRX)
instead of files. tchem gets the FASTQ addresses of the run from ENA. It reads only the start of
each file, then stops the download. Thus the size of the run is not important. In `bench/remote`, a
sample of 100,000 reads used 4.8 to 13.5 MB per run (median 11.0 MB). Each of the 20 runs took 2.6
to 3.8 seconds, but the time depends on your network.

```
$ tchem --whitelists ~/whitelists SRR19374511
Run SRR19374511 (SRR19374511), first reads streamed from ENA (5.3 MB): SRR19374511_1.fastq.gz, SRR19374511_2.fastq.gz
...
VERDICT  10x Multiome GEX  (confidence HIGH)
```

The first line of the report gives the instrument from the read names. ENA's `instrument_model`
is next to it.

For an experiment with more than one run, tchem reads the first run. An experiment is one
library, so it has one kit. tchem sets index reads aside, as it does for local files.

Some runs are only in SRA format, and ENA has no FASTQ files for them. For these runs, tchem
tells you so and gives the `fasterq-dump` command to use.

This mode needs `curl`. tchem tests all data from the network before it uses it:

- It tests the accession pattern before it puts the accession in a URL.
- It limits the size of ENA's answer and parses the answer strictly.
- It accepts only file URLs on ENA's FASTQ host.
- It limits the bytes of each download and uses only HTTPS.

`bench/remote` reads 20 held-out runs in this mode, from five kits. Each call is the same as the
call on the cached reads of the same run. The results also record the bytes that each run used.

## Commands for your pipeline

`--emit cellranger|starsolo|simpleaf|kb` prints a command for your pipeline. The command contains
the barcode and UMI positions, the barcode list, the strand and the file order. You set the
reference: the command contains a named placeholder for it.

```sh
tchem --whitelists ~/whitelists --emit starsolo S1_L001_R1_001.fastq.gz S1_L001_R2_001.fastq.gz > star.sh
STAR_INDEX=/refs/star/GRCh38 sh star.sh
```

`--emit` does not guess. If it cannot give a correct command, it gives no command and exits with
code 3. This occurs in these conditions:

- The confidence is less than HIGH. To get a command at MEDIUM or LOW, use `--emit-anyway`.
- The library type is not gene expression.
- The pipeline cannot use the kit.
- The pipeline cannot read the files correctly. For example, Cell Ranger finds R1 and R2 by
  their names, so it cannot use swapped files.

We tested each value that `--emit` prints. STARsolo, kb and simpleaf ran each command on real
reads. Cell Ranger's argument parser accepted each Cell Ranger `--chemistry` value that tchem
gives. We did not yet run Cell Ranger on these commands. [docs/emit.md](docs/emit.md) gives the
pipelines, the kits and the tests.

## Use as a Rust library

`tchem` is a small program around the `tussl_chemistry` crate. A Rust mapper or pipeline can call
the crate directly, on reads that it already holds. It needs no subprocess and no temporary files:

```rust
use tussl_chemistry::{Reads, Verdict, Whitelists, detect};

let whitelists = Whitelists::load("whitelists".as_ref())?;   // once, shared by every run
let r1 = Reads::from_file("R1.fastq.gz".as_ref(), 100_000)?;  // FASTQ or BAM, or push reads
let r2 = Reads::from_file("R2.fastq.gz".as_ref(), 100_000)?;

let report = detect(&r1, Some(&r2), &whitelists)?;
if let Verdict::Match { layout, confidence } = report.verdict() {
    // layout.geometry, layout.barcode_lists, layout.barcode_input (which file holds barcodes)
}
let strand = report.cdna_strand();   // Forward (3') / Reverse (5') for the mapper
let library = report.library();      // gene expression, VDJ, CRISPR guide, tag library
let kind = report.read_kind();       // Short or Long; long-read layouts carry long_read_kit
```

The public API is small: `Whitelists`, `Reads`, `detect`, `Report` and the types that they
return. All other items are private and can change. The crate version stays at 0.x until the API
is stable.

## Supported layouts

| Layout | Barcode read | Barcode list |
|---|---|---|
| 10x 3' v2 | R1: 16 bp barcode, 10 bp UMI | 737K-august-2016 |
| 10x 3' v3/v3.1 | R1: 16 bp barcode, 12 bp UMI | 3M-february-2018 |
| 10x 3' v4 (GEM-X) | R1: 16 bp barcode, 12 bp UMI | 3M-3pgex-may-2023 |
| 10x 5' v1/v2 | R1: 16 bp barcode, 10 bp UMI, TSO | 737K-august-2016 |
| 10x 5' v3 (GEM-X) | R1: 16 bp barcode, 12 bp UMI, TSO | 3M-5pgex-jan-2023 |
| 10x Multiome GEX | R1: 16 bp barcode, 12 bp UMI | 737K-arc-v1 (gene expression) |
| 10x Flex v1 | R1: 16 bp barcode, 12 bp UMI; R2 reads probes | 737K-fixed-rna-profiling (from Cell Ranger) |
| BD Rhapsody (original V1 beads) | R1: three 9 bp cell labels split by linkers | BD CLS1-3 |
| BD Rhapsody (Enhanced Beads) | R1: as V1, starting 0 to 3 bases in | BD CLS1-3 |
| BD Rhapsody (Enhanced V2/V3 beads) | R1: as Enhanced, 384 sequences per cell label | BD V2 CLS1-3 |
| SPLiT-seq (Rosenberg 2018) | R2: UMI, three 8 bp round barcodes, two linkers | rounds 1-3 |
| Parse Evercode WT | R2: as SPLiT-seq, Parse's second linker (WT v1, v2, early v3) | rounds 2-3 |
| Parse Evercode WT v3 | R2: UMI, three 8 bp rounds, 12 bp linkers ATGAGGGGTCAG and TCCAACCACCTC | Parse's v3 tables, rounds 1-3 |
| Singleron GEXSCOPE V1 | R1: three 8 bp barcodes split by 16 bp linkers (4 alternatives each), C, 12 bp UMI | one list of 96 |
| Singleron GEXSCOPE V2 | R1: as V1 with 9 bp barcodes | bc1-3 (96, 96, 192) |
| Singleron GEXSCOPE V3 | R1: 0 to 3 bases in, 9 bp barcodes split by ACGATG and CATAGT, C, 12 bp UMI | bc1-3 |
| MGI DNBelab C4 / C-series HT | R1: two 10 bp barcodes, 10 bp UMI | MGI C-series (1,536 per barcode) |
| Illumina Single Cell 3' / PIPseq v3-v4 | R1: four barcode parts split by ATG, GAG, TCGAG; 12 bp UMI, polyT; offset 0 to 3 | none (counted) |
| Illumina Single Cell 3' RNA Prep / PIPseq V | R1: four barcode parts split by ATG, GAG, TCGAG; 3 bp binning index, ACG; offset 0 to 3 | none (counted) |
| ScaleBio RNA | R1: barcode, linker, UMI, barcode, offset 0 to 1 | none (counted) |
| ScaleBio QuantumScale RNA | R2: 8 bp UMI, 8 bp RT barcode | QuantumScale RT |
| SeekGene SeekOne DD 3' | R1: 17 bp barcode (7 counted + 10 listed), 12 bp UMI, polyT | SeekGene P3CB, bases 8-17 |
| SeekGene SeekOne DD 5' | R1: as DD 3', TSO | SeekGene P3CB, bases 8-17 |
| inDrops | R1: two barcode halves around the W1 adapter, offset 0 to 3 | none (counted) |
| Drop-seq | R1: 12 bp barcode, 8 bp UMI | none (counted) |

"None (counted)" means the kit has no public barcode list, so tchem counts the barcodes (see
[How tchem decides](#how-tchem-decides)). On long reads, tchem detects the five 10x
gene-expression layouts: 3' v2, 3' v3/v3.1, 3' v4, 5' v1/v2 and 5' v3. It uses the same table.

The layouts are data in `src/layouts/`, with one file for each vendor. The compiler tests the
table. To add a kit, follow [the kit recipe](ROADMAP.md#the-kit-recipe). The level of each layout
is in [How well it works](#how-well-it-works).

## How tchem decides

tchem reads the first reads of each file. It compares each known layout with them, in both file
orders. For each layout, it measures these signals:

1. **Listed barcodes**: the share of reads with a barcode on the vendor's list. A barcode can
   have one different base. An `N` matches any base. Cell Ranger and alevin-fry use the same
   tolerance.
2. **Linkers**: the share of reads with each fixed sequence at its position. tchem permits one
   different base in each five bases.
3. **Unlisted barcodes**, for kits without a public list. Real cell barcodes put most reads on a
   few thousand sequences. Random sequence does not. tchem measures the share of reads on the most
   frequent sequences.
4. **UMI randomness**: each UMI position must have a Shannon entropy of at least 1.5 bits.
5. **The tail after the UMI**: polyT for 3' kits, or the template-switch oligo (TSO) for 5' kits.
   tchem tests the tail only when the read is long enough.
6. **Start offsets**, for beads whose barcode starts a variable number of bases into the read.

The support of a layout is the mean of its listed-barcode, linker and unlisted-barcode signals.
tchem calculates each signal over all reads in the sample. A read that is too short for a segment
does not show that segment. A contradiction is a thing in the reads that the layout does not
permit. Examples are a missing linker, a read that is too short, or a UMI that is not random. The
report lists the contradictions for each layout.

tchem then gives one of these results:

- **HIGH** needs these four conditions:
  - The support is at least 70%.
  - There is no contradiction.
  - No other layout has half as much support.
  - A barcode list or a linker gives evidence. Unlisted barcodes alone never give HIGH.

  A kit with no barcode list also needs each of its linkers in at least 70% of reads. This stops
  a HIGH call on a file that mixes two libraries. List evidence counts only when it is spread
  over many barcodes, as in a real library. No single sequence can hold half of the hits. A
  constant read, such as the amplicon of a tag library, can be near a list entry by chance.
- **MEDIUM** is a clear lead over the other layouts, but less than HIGH needs.
- **LOW** is not a weaker MEDIUM. At least 80% of the barcodes are on the list, spread as above.
  But the reads contradict another part of the layout. The report shows
  the contradictions.
- **Ambiguous** means that two or more layouts fit the reads equally well, and a pipeline must
  process them differently.

10x 3' and 5' kits share a barcode list. A 26 bp R1 ends at the UMI, so it cannot show polyT or
the TSO. In this case, tchem measures the strand of the cDNA read without a reference:

- polyA runs against polyT runs. A 3' cDNA read continues into the polyA tail.
- the human and mouse mitochondrial genomes, which are in tchem. Their transcripts come mostly
  from one strand.
- the spliced transcripts of 40 highly expressed human nuclear genes.

Each signal that prefers a strand must agree with the others. If they agree, tchem tells 3' from
5'. If they disagree, or no signal prefers a strand, the verdict is ambiguous. The report gives
the cause.

The library type comes from the cDNA read:

- VDJ: reads from the constant regions of TCR and BCR genes.
- CRISPR guide: reads from the sgRNA scaffold.
- Tag library (antibody, hashtag, guide or barcode amplicon): a small set of read starts that
  covers almost all reads.

We set these thresholds on 480 public experiments whose ENA titles give the library type.

When no layout matches, the report tells you what the files show. The reads can have a fixed
structure, which indicates a kit that tchem does not know. The reads can all have cDNA length,
which indicates that the barcode read is missing. In that case, get the SRA copy with
`fasterq-dump --include-technical`. Or the reads can be long reads.

### Long reads

A long read contains the full molecule, in either direction. Nanopore errors are mostly
insertions and deletions, so no part is at a fixed position. tchem treats a sample as long reads
when at least one fifth of its reads are longer than 301 bases. No short read is that long. tchem
then searches each read in three steps:

1. Find 10x's Read 1 adapter at either end, in either direction, with up to 4 edits.
2. Find the 16-base barcode after the adapter in each 10x list. Permit a shift of two bases.
3. Find polyT (3') or the TSO (5') after the UMI. Measure the UMI length (10 or 12 bases).

The edit distance uses Myers' bit-parallel algorithm. A test compares it with the textbook
algorithm. The 10x lists share some barcodes. A shared barcode is not evidence for one list, so
tchem ignores it.

The kit is the one whose list, tail and UMI length agree. It comes from the same layout table as
for short reads. A read that is too short to hold the start of a molecule is not evidence
against a kit.

## How well it works

**Summary.** We ran tchem on 10,130 held-out public runs from 5,212 studies. A held-out run is a
run whose reads we did not use to make or change tchem. tchem made 8,923 HIGH calls. None of them
is wrong. The 8,923 HIGH calls are in four groups:

| Group | HIGH calls | Meaning |
|---|---|---|
| Right, exact | 6,694 (3,563 studies) | An independent label names the same layout. |
| Right, family | 995 | An independent label names only the kit family, and the call is in that family. |
| Read-derived label only | 235 | Only a label from the reads agrees. A level does not count these calls. |
| Disputed | 999 | The label names a different kit, but the reads' evidence supports the call. |

An independent label comes from a source that did not read the reads. Examples are the run's
protocol text, the lab's record, the vendor's statement or a BD part number. A disputed call is
not right and not wrong. [Disputes](#disputes) gives the evidence.

### What a level means for you

Each match prints the level of its layout. A level counts held-out studies, not runs. Runs from
one study have the same lab and protocol, so they are not independent trials. The error rate
comes from the rule of three. If n independent studies have no wrong call, the error rate is less
than 3/n at 95% confidence.

| Level | The record | What to do with a HIGH call |
|---|---|---|
| `error_under_3pct` | 100 or more held-out studies, none wrong | Use it. You can use it without a person. |
| `validated` | 30 to 99 studies, none wrong | Use it. Compare the first runs from a new source with the kit you used. |
| `checked` | 1 to 29 studies, none wrong | Confirm the kit before a large run. The error bound is too wide to claim a rate. |
| `unvalidated` | No held-out HIGH call of this exact layout | Confirm the kit before you use the call. |
| `contradicted` | A held-out HIGH call was wrong | Do not use the call. |

Below HIGH (MEDIUM, LOW or ambiguous), the report shows the contradictions and the other
candidates. Read them, and confirm the kit before use. `--emit` gives no command below HIGH,
unless you use `--emit-anyway`.

### The record for each layout

`bench/coverage` made this table from one coverage pass of the released binary over all held-out
runs. [bench/VALIDATION.md](bench/VALIDATION.md) gives the procedure. For a `checked` layout, the
table shows the error bound, but the bound is too wide to claim. The bound has a maximum of 100%.
A layout that is not in the table has no held-out HIGH call.

<!-- coverage:begin -->
| Layout | Reads | Level | Held-out HIGH calls (studies) | Error bound (95%) | Kit family only | Read-derived label only | Disputed | Wrong |
|---|---|---|---|---|---|---|---|---|
| 10x 3' v2 | short | error_under_3pct | 626 (335) | 0.9% | 0 | 27 | 85 | 0 |
| 10x 3' v3/v3.1 | long | checked | 11 (11) | 27.3% | 0 | 0 | 0 | 0 |
| 10x 3' v3/v3.1 | short | error_under_3pct | 4346 (2282) | 0.1% | 0 | 140 | 637 | 0 |
| 10x 3' v4 (GEM-X) | short | validated | 172 (87) | 3.4% | 0 | 2 | 77 | 0 |
| 10x 5' v1/v2 | long | checked | 3 (2) | 100.0% | 0 | 0 | 0 | 0 |
| 10x 5' v1/v2 | short | error_under_3pct | 747 (416) | 0.7% | 0 | 64 | 110 | 0 |
| 10x 5' v3 (GEM-X) | short | checked | 36 (17) | 17.6% | 0 | 1 | 11 | 0 |
| 10x Flex v1 | short | error_under_3pct | 217 (118) | 2.5% | 0 | 0 | 35 | 0 |
| 10x Multiome GEX | short | error_under_3pct | 218 (135) | 2.2% | 0 | 1 | 7 | 0 |
| BD Rhapsody (Enhanced Beads) | short | checked | 9 (1) | 100.0% | 114 | 0 | 1 | 0 |
| BD Rhapsody (Enhanced V2/V3 beads) | short | checked | 2 (1) | 100.0% | 104 | 0 | 6 | 0 |
| BD Rhapsody (original V1 beads) | short | checked | 18 (8) | 37.5% | 221 | 0 | 0 | 0 |
| Illumina Single Cell 3' / PIPseq v3-v4 | short | checked | 45 (24) | 12.5% | 0 | 0 | 0 | 0 |
| Illumina Single Cell 3' RNA Prep / PIPseq V | short | checked | 10 (5) | 60.0% | 0 | 0 | 2 | 0 |
| MGI DNBelab C4 / C-series HT | short | validated | 72 (46) | 6.5% | 0 | 0 | 2 | 0 |
| Parse Evercode WT | short | unvalidated | 0 (0) | - | 152 | 0 | 1 | 0 |
| Parse Evercode WT v3 | short | unvalidated | 0 (0) | - | 50 | 0 | 0 | 0 |
| ScaleBio RNA | short | checked | 4 (1) | 100.0% | 0 | 0 | 0 | 0 |
| SeekGene SeekOne DD 3' | short | checked | 51 (26) | 11.5% | 15 | 0 | 15 | 0 |
| SeekGene SeekOne DD 5' | short | checked | 6 (3) | 100.0% | 0 | 0 | 0 | 0 |
| Singleron GEXSCOPE V1 | short | unvalidated | 0 (0) | - | 97 | 0 | 1 | 0 |
| Singleron GEXSCOPE V2 | short | checked | 16 (8) | 37.5% | 238 | 0 | 7 | 0 |
| Singleron GEXSCOPE V3 | short | checked | 8 (2) | 100.0% | 4 | 0 | 0 | 0 |
| inDrops | short | validated | 77 (42) | 7.1% | 0 | 0 | 2 | 0 |
<!-- coverage:end -->

### The record for each kit family

For BD Rhapsody, Parse Evercode, SeekGene SeekOne DD and Singleron GEXSCOPE, tchem finds the bead
type or version from the reads. But labs seldom record it. Thus most of their runs have only a
family label: the "Kit family only" column above.

The family record answers this question. When tchem names a layout of the family at HIGH, is the
run from that family? The rules are the same as for layouts. A call that names the wrong layout in
the correct family is right here. The table above shows that call as disputed or wrong.

<!-- family:begin -->
| Kit family | Layouts | Reads | Level | Held-out HIGH calls (studies) | Error bound (95%) | Read-derived label only | Disputed | Wrong |
|---|---|---|---|---|---|---|---|---|
| BD Rhapsody | 3 | short | error_under_3pct | 473 (252) | 1.2% | 0 | 2 | 0 |
| Parse Evercode WT | 2 | short | error_under_3pct | 202 (109) | 2.8% | 0 | 1 | 0 |
| SeekGene SeekOne DD | 2 | short | validated | 72 (38) | 7.9% | 0 | 15 | 0 |
| Singleron GEXSCOPE | 3 | short | error_under_3pct | 363 (187) | 1.6% | 0 | 8 | 0 |
<!-- family:end -->

### Disputes

999 HIGH calls disagree with the label of their run. Most of these labels name the wrong 10x
version. Where the two kits use different barcode lists, the reads match the list of the called
kit. The match is at least 74%, and usually 88 to 100%. The reads match the list of the labelled
kit only weakly.

Three other methods read the files of each dispute. We wrote the plan before they ran. No method
agreed with the label in any dispute (`bench/second_opinion`):

- Sanger's `infer_10x_chemistry.py` agreed with tchem's barcode list in 851 of 914 disputes.
- A strand referee used genes that tchem does not use. It agreed with tchem in 77 of the 127
  disputes between 3' and 5'.
- Cell Ranger's `--chemistry auto` agreed with tchem in all 40 disputes of a random selection of
  human runs.

These methods read the same reads as tchem. Thus they support a call, but they do not change a
level. `bench/coverage/results/adjudication.tsv` gives the evidence for each dispute.

### Recall

Recall shows what tchem said on the held-out runs, for each independent label. A level does not
use recall. The label is a layout, or a kit family when the label names only the family. The
table shows labels with 20 runs or more. `bench/coverage/results/recall.tsv` has all labels.
`bench/coverage/scripts/recall.py` makes this table, and the gate compares it with the coverage
pass.

<!-- recall:begin -->
| Label | Read | Runs (studies) | HIGH, as labelled | HIGH, another kit (disputed) | MEDIUM or LOW | Ambiguous | No call |
|---|---|---|---|---|---|---|---|
| 10x 3' v3/v3.1 | short | 4776 (2484) | 91% | 5% | 2% | 0% | 2% |
| 10x 3' v2 | short | 1308 (686) | 48% | 46% | 2% | 1% | 3% |
| 10x 5' v1/v2 | short | 1014 (527) | 74% | 4% | 2% | 18% | 2% |
| 10x Multiome GEX | short | 555 (284) | 41% | 1% | 0% | 0% | 58% |
| BD Rhapsody (family) | short | 480 (260) | 91% | 1% | 3% | 0% | 5% |
| Singleron GEXSCOPE (family) | short | 413 (217) | 82% | 0% | 14% | 0% | 4% |
| Parse Evercode WT (family) | short | 263 (139) | 78% | 1% | 18% | 0% | 4% |
| 10x Flex v1 | short | 247 (131) | 88% | 6% | 0% | 0% | 6% |
| 10x 3' v4 (GEM-X) | short | 197 (99) | 87% | 10% | 1% | 0% | 2% |
| inDrops | short | 128 (65) | 60% | 0% | 2% | 0% | 38% |
| MGI DNBelab C4 / C-series HT | short | 98 (50) | 73% | 2% | 23% | 0% | 1% |
| Illumina Single Cell 3' / PIPseq v3-v4 | short | 66 (34) | 68% | 3% | 0% | 0% | 29% |
| 10x 5' v3 (GEM-X) | short | 58 (31) | 52% | 36% | 0% | 9% | 3% |
| SeekGene SeekOne DD 3' | short | 57 (29) | 89% | 0% | 4% | 0% | 7% |
| BD Rhapsody (original V1 beads) | short | 24 (10) | 79% | 17% | 4% | 0% | 0% |
<!-- recall:end -->

"HIGH, another kit" is the disputed column of the layout record. None of these calls is wrong,
and the second opinions agree with tchem. A high value shows that the archive labels of that kit
often name the wrong version. It does not measure the accuracy of tchem.

A MEDIUM or LOW call also names a kit:

- For MGI and Parse, each MEDIUM or LOW call names the labelled kit. The confidence is LOW because
  the barcode list agrees, but the reads contradict part of the layout. The report shows both.
- For Singleron, 20 runs get the labelled kit at LOW. 37 runs get Drop-seq at MEDIUM. Singleron's
  GEXSCOPE MicroBead kit has the Drop-seq structure ([docs/validation.md](docs/validation.md)).

### Known gaps

We found these gaps on the held-out runs on 10 Oct 2026. They are not yet fixed:

- **PIPseq v3-v4:** 19 of 66 labelled runs get no call. Their barcode parts and linkers fit at
  96 to 99%. But the polyT starts some bases later than the layout permits. A missing polyT is a
  contradiction for a layout with no barcode list. The report says: "the reads contradict every
  candidate".
- **inDrops:** 49 of 128 labelled runs get no call. In the run that we examined, the barcode read
  is in the second file. It is also shorter than the layout needs.
- **10x Multiome GEX:** 322 of 555 labelled runs get no call. The label comes from protocol text,
  which describes both libraries of a Multiome experiment. ENA records 267 of the 555 runs as
  ATAC-seq, and tchem gives no call on 265 of them. ENA records 194 runs as RNA-Seq or ssRNA-seq,
  and tchem gives a match on 182 of them (`bench/coverage/results/multiome_strategy.tsv`).
- **10x 5' v1/v2:** 18% of the runs are ambiguous. We know the read lengths of 131 of these runs.
  In 120 of them, read 1 has 26 to 28 bases. This read ends at or just after the UMI, before the
  TSO. 3' v2 has the same barcode list. When the cDNA read shows no strand, tchem cannot tell
  these kits apart.

A fix for PIPseq or inDrops changes detection. Thus the fix waits for a release with a new
coverage pass. The design-data rule applies: the runs that we use for the fix become design runs.

### Comparison with Sanger's tool

We compared tchem with the Wellcome Sanger Institute's `infer_10x_chemistry.py` (Cellular
Genetics). It is the most similar existing tool. The test used 483 held-out public experiments from
Arc Institute's scBaseCount (`bench/head_to_head`). We froze the selection plan before either tool
ran. Both tools can process 445 of the gene-expression experiments:

| | tchem | Sanger's tool |
|---|---|---|
| Barcode list the same as scBaseCount's | 445 / 445 | 426 / 445 (19 no call, 0 different) |
| 3' and 5' told apart | 387 (87%) | 226 (51%) |
| Correct, of the tool's calls where the truth is known | 323 / 323 | 181 / 181 |
| CPU time per experiment (median of 72, one at a time) | 1.35 s | 7.08 s |

scBaseCount's barcode list is itself found from the reads, so the first row shows agreement, not
proof. The truth for 3' and 5' comes from ENA titles and a strand referee on nuclear genes. Each
tool told 3' from 5' on a different set of experiments, so the totals in the third row differ.

No HIGH call from tchem was wrong. The library type was correct for all 42 CRISPR guide libraries
with a title, and for 53 of 57 VDJ libraries with a title. We set the library-type thresholds on
these same experiments (`bench/library_type`), so these two numbers are not from held-out data. We measured the CPU time of each tool with
`wait4` on cached files (`bench/head_to_head/scripts/timing.py`, `results/timing.tsv`).

[docs/validation.md](docs/validation.md) tells you what these numbers show and do not show. It
also describes each benchmark when we first ran it. [bench/VALIDATION.md](bench/VALIDATION.md)
gives the procedure for each kit.

## Engineering

tchem is written in Rust to [TigerStyle](https://tigerstyle.dev) rules. A test (`tests/tidy.rs`)
requires at least two assertions in each function that run in a release build. 60 functions are
exempt, each for a stated reason. 16 functions on the per-read path have checks that run only
in the debug build. The gate runs that build on all parity oracle inputs.

Release builds keep their assertions and stop at the first failure. tchem gives a typed refusal for bad input. It does not
panic.

Each commit must pass the gate. The gate runs the tests, 281 parity oracle inputs, 425 chaos
cases and 50 probes. We also fuzz the readers and detection (`fuzz/`). tchem has one dependency:
`flate2`, for gzip. [CONTRIBUTING.md](CONTRIBUTING.md) gives the rules, the gate and the steps to
repeat each test.

## Data

- `data/reference/human_chrM.seq` and `mouse_chrM.seq`: NCBI RefSeq NC_012920.1 (rCRS) and
  NC_005089.1.
- `data/reference/housekeeping_human.fa`: 607 spliced transcripts of 40 genes. Also
  `tcr_bcr_constant_human.fa`: TCR and BCR constant regions. Both come from GENCODE v44 on GRCh38
  (`bench/head_to_head/scripts/build_nuclear_ref.py`).
- TusslChemistry does not include barcode lists. `data/lists.tsv` is the register of all lists.
  For each list, it gives the source (at a fixed upstream commit) and a second source. It also
  gives the barcode count and length, the sha256 of the file, the owner and the terms.
- `scripts/fetch_whitelists.sh` downloads each list from the register. It refuses a file that
  does not match.
- Flex's list is only in Cell Ranger. You download Cell Ranger from 10x under 10x's licence.
  `--cellranger <install dir>` copies the list from there. Without it, the script gets all other
  lists, and tchem reports Flex as not tested.
- 10x's barcode lists belong to 10x, under 10x's terms. Do not distribute them. Do not distribute
  the `.packed` files that `--pack` makes from them.

## Questions, citing and security

Ask a question in [Discussions](https://github.com/Tussl-Co/TusslChemistry/discussions).
Report a wrong answer or a bug as an issue. Report a security problem privately, as
[SECURITY.md](SECURITY.md) describes.

If you use tchem in published work, please cite it. [CITATION.cff](CITATION.cff) has the details.
GitHub's "Cite this repository" gives them as APA or BibTeX.

## Licence

Copyright 2026 Tussl Ltd. Licensed under the [Apache License, Version 2.0](LICENSE); see
[NOTICE](NOTICE).
