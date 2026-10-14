# Ready commands: `--emit`

`tchem --emit <pipeline>` prints a command for a downstream pipeline instead of the report. The
command holds only what the reads decide:

- where the barcode and the UMI are,
- which barcode list to correct against,
- which strand to count,
- which file is the barcode read and which is the cDNA.

You supply the reference. Each value that only you know is a named placeholder. In the shell
form, a placeholder reads an environment variable, and the command stops with a sentence if the
variable is not set:

```sh
$ tchem --whitelists ~/whitelists --emit starsolo S1_L001_R1_001.fastq.gz S1_L001_R2_001.fastq.gz
# tchem 0.1.0: 10x 3' v3/v3.1 (confidence HIGH)
# validation: error_under_3pct: 4348 held-out HIGH calls from 2282 studies, none wrong (error under 0.1% at 95%); ...
# set STAR_INDEX: a STAR genome index built with a gene annotation (GTF)
STAR \
  --genomeDir "${STAR_INDEX:?set STAR_INDEX to a STAR genome index built with a gene annotation (GTF)}" \
  --readFilesIn S1_L001_R2_001.fastq.gz S1_L001_R1_001.fastq.gz \
  --readFilesCommand 'gunzip -c' \
  --soloType CB_UMI_Simple \
  ...
$ tchem --whitelists ~/whitelists --emit starsolo S1_L001_R1_001.fastq.gz S1_L001_R2_001.fastq.gz > command.sh
$ export STAR_INDEX=/refs/star/GRCh38 && sh command.sh
```

With `--format json`, the command is one object (see [schema/emit.v1.json](../schema/emit.v1.json)).
`argv` is the argument list, and each `${NAME}` item stands for the matching entry in
`placeholders`. `shell` is the script above.

## Pipelines and kits

| Pipeline | Kits it gets a command for | Fills in |
|---|---|---|
| `cellranger` | 10x 3' v2/v3/v4, 5' v1/v2/v3, Multiome GEX | `--chemistry`, `--fastqs`, `--sample`, `--id`; `CELLRANGER_REF` |
| `starsolo` | 10x 3' and 5' kits, Multiome GEX, Drop-seq | `CB_UMI_Simple` positions, the list (`None` for Drop-seq), `--soloStrand`, file order; `STAR_INDEX` |
| `simpleaf` | the same | the geometry, `--unfiltered-pl` or `--knee`, `--expected-ori`; `SIMPLEAF_INDEX` |
| `kb` | the same | the `-x` technology string, `-w`, `--strand`; `KB_INDEX`, `KB_T2G` |

Multiome GEX is always nuclei. The RNA of nuclei is mostly unspliced, and Cell Ranger counts
its intronic reads. Thus, for Multiome:

- The STARsolo command counts `--soloFeatures GeneFull`.
- The kb command uses the `nac` workflow (`--workflow nac`, `--sum=total`). This workflow also
  needs `KB_C1` and `KB_C2` from `kb ref --workflow nac`.
- The simpleaf command does not change. It has a note: add its spliced, unspliced and ambiguous
  counts.

`--emit` refuses the other kits for now: BD, Parse, SPLiT-seq, Singleron, MGI, PIPseq, ScaleBio
and inDrops. Their barcodes are in several segments. These four pipelines accept such barcodes
only in forms that tchem does not write yet (STARsolo `CB_UMI_Complex`, combined lists). Use the
vendor's own pipeline for those kits.

## When no command is written

A wrong kit does not make a pipeline fail. The run finishes with fewer cells or the wrong
genes. So `--emit` refuses rather than guesses, and exits with **3**. In each of these cases,
stderr says why:

- The verdict is not a HIGH match. `--emit-anyway` emits a MEDIUM or LOW match all the same, and
  adds a note to the command saying so.
- The library is not gene expression (VDJ, CRISPR guides, antibody tags).
- The reads are long reads. Use the report's `wf-single-cell --kit` instead.
- The kit is 10x Flex, whose second read is a probe and not cDNA. Use `cellranger multi` with its
  probe set.
- The cDNA strand could not be measured.
- Only one file was given.
- For Cell Ranger, the files are swapped, or they do not have Illumina names
  (`<sample>_S1_L001_R1_001.fastq.gz`). Cell Ranger finds R1 and R2 by name. Archive downloads
  frequently do not have Illumina names. The other pipelines take files in an explicit
  sequence, so they accept swapped files.

## How the commands were tested

`bench/emit` runs every command on real reads (200,000 per case), on a small reference built
from 40 housekeeping genes. The results are in
[bench/emit/results/results.tsv](../bench/emit/results/results.tsv). What each check does and
does not cover:

- **STARsolo 2.7.11b, kb 0.30.2 and simpleaf 0.30.0** run to the end on six kits. The kits are
  10x 3' v3 (also with swapped files), 3' v4, 5' v1/v2, 5' v3, Multiome GEX and Drop-seq. kb
  finds 95 to 100% of barcodes on the list.
  - The bench ran STARsolo with `--readFilesCommand zcat`. tchem now writes `gunzip -c`,
    because macOS's `zcat` reads only `.Z` files. On 12 Oct 2026, STARsolo 2.7.11b ran both on
    a whole 10x 3' v3 lane (33,436,697 reads). The two `Solo.out` directories were identical.
  - On two of the 10,000-cell runs, simpleaf needs alevin-fry's read floor lowered to 1. This is
    because only about 5% of reads land on the small reference. The bench lowers it and records
    that it did; the emitted command never does.
- **Strand.** Each STARsolo run is repeated with the opposite strand. With the emitted strand,
  2 to 40 times more reads are assigned to genes on six of the seven runs.
  - The seventh is the Multiome nuclei run, which is inconclusive on this reference (3.1% against
    3.4%). Its emitted strand is forward. Both Cell Ranger's own definition of `ARC-v1` and
    tchem's strand measurement on that run agree.
- **Cell Ranger 10.1.0's argument parser accepts all six `--chemistry` values** in the layout
  table. The parser refuses an unknown value before anything runs. Cell Ranger itself was not
  run: that needs its full install and a reference.
- **Not tested yet.** Only STARsolo has a strand control. A wrong UMI length would still let a
  pipeline finish and assign reads, so no check here would catch it. The Multiome commands'
  intron counting (STARsolo `GeneFull`, kb `--workflow nac`) has had only one test: kb's
  argument parser (`kb count --dry-run`). `bench/emit` has not run it yet.
