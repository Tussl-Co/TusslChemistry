# Roadmap

This file tells where tussl-chemistry goes, in which sequence, and what "done" means at each
step. We update it when a milestone moves, not for each commit. The terms in this file have the
meanings in [docs/glossary.md](docs/glossary.md).

**The aim.** One cross-vendor step tells a pipeline what a pair of FASTQ files is. It gives the
kit, the barcode list, the strand, the library type and the instrument. It says how far to trust
the answer, and it gives a ready command.

One engine measures evidence in the reads. The kits are data, so a new kit is mostly knowledge
and not new code. That knowledge is the read structure, the barcode lists, and runs that prove
the kit.

Structures marked "to verify" are not yet tested on real reads. Vendor documents and archive
metadata are sometimes wrong. For example, 17 of 223 public studies that name Singleron did not
contain Singleron libraries. For this reason, we claim nothing until the reads confirm it.

## Where things stand (10 Oct 2026)

- 25 short-read layouts from 10x, BD, Parse, SPLiT-seq, MGI, PIPseq (v3-v4 and V), ScaleBio,
  Singleron, SeekGene (DD 3' and 5'), inDrops and Drop-seq.
- Long reads for the five 10x gene-expression kits.
- The held-out record: 8,923 HIGH calls on 10,130 public runs from 5,212 studies. None is wrong,
  and 999 are disputed (README, "How well it works").
- At `error_under_3pct`: 10x 3' v3, 5' v1/v2, 3' v2, Multiome GEX and Flex v1. Also the BD
  Rhapsody, Singleron GEXSCOPE and Parse Evercode WT kit families.
- At `validated`: 10x 3' v4, MGI, inDrops and the SeekGene DD kit family.
- Every commit passes the gate: the tests, the parity oracle, `validate.sh`, 425 chaos cases and
  50 probes. The readers and detection are fuzzed.

## The kit recipe

This is the one procedure for a kit. Every kit goes through these steps before tchem names it:

1. Get the read structure from the primary source: the vendor's own pipeline or user guide.
2. Get the barcode lists from a durable source.
3. Enter each list in the list registry (`src/lists/`, rendered to `data/lists.tsv`).
4. Give each list an exact count, a checksum and its terms.
5. Make sure that `scripts/fetch_whitelists.sh` fetches each list. Do not ship a list.
6. Find public runs. Verify each run from its own reads, never from metadata only.
7. Use about one run for each kit version as a design run. The other runs are held-out runs.
8. Add the layout, with the compile-time table checks.
9. Make sure that the gate passes and the parity oracle explains every change.
10. Do a coverage pass on every held-out run (`bench/coverage`), as `bench/VALIDATION.md` says.
11. Make sure that every match prints the level of the layout from that pass.
12. Add a permanent check to `scripts/validate.sh`.

A held-out run that you use to change the code becomes a design run (the design-data rule). Do
this before the next coverage pass.

Each level comes from the rule of three. If n independent trials have no HIGH call wrong, the
error rate is less than 3/n, at 95% confidence. We count the trials in **studies**, not runs.
Runs from one study share a lab, a protocol and often a library. Every match prints the level of
its layout (`layout.validation`):

| Held-out studies | Level |
|---|---|
| 100 or more, no HIGH call wrong | `error_under_3pct` |
| 30 or more, no HIGH call wrong | `validated`, error under 10% |
| 1 to 29, no HIGH call wrong | `checked`: the report gives the counts and the error bound, but claims no rate |
| none | `unvalidated` |
| any, with a HIGH call wrong | `contradicted` |

## Next

### Trust: validate every kit the same way

- **0.1.1: PIPseq v3-v4 and inDrops recall.** 19 of 66 PIPseq runs and 49 of 128 inDrops runs
  get no call (README, "Recall").
  - In some PIPseq runs, the polyT starts a few bases late.
  - In some inDrops runs, the barcode read is in the second file. It is shorter than the layout
    needs.
  - A fix to either kit changes detection. The runs that we use become design runs (the
    design-data rule), and we do the coverage pass again.
- **Assertions in detection, with the same coverage pass.** A change to a detection file makes
  the held-out record old, so these changes wait for that pass:
  - Remove the assertions that only repeat the line above them. A code review on 10 Oct 2026
    found them in `input/bam.rs`, `cli/mod.rs`, `lib.rs`, `signals/structure.rs`,
    `layouts/mod.rs` and `lists/mod.rs`.
  - Give a release-build check to each entry in `DEBUG_ONLY` (`tests/tidy.rs`) where the check
    costs little.
  - Add postconditions to the parsers of untrusted names (`signals/instrument.rs`,
    `input/fastq.rs` `decoded`, `input/bam.rs` `record_name`).
  - Put the last printed words onto the glossary: the words `CeleScope chemistry` in the
    Singleron notes (`layouts/singleron.rs`).
  - Make the unit tests of `input/fastq.rs` and `input/whitelist.rs` remove their temporary
    files, as `tests/api.rs` does.
- **More held-out studies.** A third draw from the ENA catalogue (`bench/tenx_ena`) is possible.
  2,618 studies remain, in seeded sequence. This draw would take 10x 3' v4 from 87 studies to
  near 100. MGI and PIPseq have little more public data. ScaleBio has 6 studies at most, and
  QuantumScale has none.
- **Labels.** These labelling rules are in doubt. We must correct them before the next draw:
  - We read "GEM-X Flex" as Flex v2.
  - We catalogue a study that names two assays as one assay.
  - Multiome labels from protocol text also cover the ATAC library.
- Do the coverage pass again for each release. `scripts/release_check.sh` already refuses a
  release with an old record.

### Coverage: more kits

Some kits need a larger layout model first:

| Enabler | Needed for |
|---|---|
| Barcodes in an index read (I1/I2) | 10x ATAC and Multiome ATAC (barcode in i5), inDrops v3, some plate methods |
| Barcodes in R2, probe barcodes | 10x Flex v2, some combinatorial kits |
| Variable-length barcodes | inDrops (8 to 11 bases), ScaleBio, PIPseq tiers (to verify) |
| Kit families for identical structures | Drop-seq, Seq-Well, Singleron MicroBead reported as one family, members named |

Frozen catalogues are already available for five kits (`bench/{visium,celseq2,atac,ddseq,icell8}`).
Each kit needs a layout, and then a held-out record:

- Barcodes from an index read. This makes 10x scATAC and Multiome ATAC possible.
- Visium (v1 and v2 slides, CytAssist). Its lists come from a Space Ranger install.
- CEL-seq2, SORT-seq and MARS-seq. Some "CEL-seq2" studies use a 6+8 / 384-barcode layout. MARS-seq
  has its barcode on R2.
- Bio-Rad ddSEQ (start offsets 0 to 5), then Visium HD (a two-part barcode layout).
- Takara ICELL8. It has no public barcode list, so we can detect only its structure.
- ddSEQ and ICELL8 have no held-out runs yet. They have only design runs until new studies appear.

Then, approximately in the sequence of how frequent each kit is in public data:

- 10x GEM-X Flex v2. It has its own list and a 10-base probe barcode.
- 10x 3' v1. It is old, but it is still in the archives.
- 10x On-Chip Multiplexing and CellPlex: the library type and the sample tags.
- Parse Evercode WT v4, when public runs are available. The reads cannot show the kit size (Mini,
  WT or Mega).
- SeekGene MM. Its list possibly counts now, because a list counts only when its hits are spread.
- Honeycomb HIVE (no list).
- BD Rhapsody targeted panels, AbSeq and Sample Tag library types. TCR/BCR Next (no public runs).
- Seq-Well and Seq-Well S3, which we expect in the Drop-seq family. All to verify.
- sci-RNA-seq3, Smart-seq3 and Smart-seq3xpress, VASA-seq, Quartz-seq2, Microwell-seq, M20
  Genomics and MobiDrop. All to verify.
- Spatial sequencing, where barcodes give a position: 10x Visium HD, BGI Stereo-seq, Curio Seeker
  (Slide-seq) and Open-ST.

### Other modalities on the same engine

- scATAC: 10x ATAC v1 and v2, and Multiome ATAC. These need barcodes in an index read.
- Feature libraries in addition to the library typing of today: CITE-seq antibodies, cell
  hashing, Perturb-seq guides and VDJ.
- Long reads from more vendors, when public data appears: Parse, ScaleBio and BD cDNA on Oxford
  Nanopore or PacBio.

### Acting on the answer

- Run `bench/emit` again with STARsolo, kb and simpleaf installed. Include the Multiome intron
  counting (GeneFull, kb nac). Until now, only the kb parser has tested it.
- Add a test that finds a wrong UMI length. Add strand controls for kb and simpleaf.
- Run Cell Ranger itself on a small reference, not only its argument parser.
- Add more pipelines: CeleScope, Parse's split-pipe and Oxford Nanopore's `wf-single-cell`.
- seqspec is the IGVF standard that describes the read structures of assays.
  - Import a spec as a candidate layout. It must go through the kit recipe before we claim it.
  - Export a seqspec for each detected layout.

### Distribution and upkeep

- Publish: tag `v0.1.0`, then submit the bioconda recipe and the nf-core module. Both are in
  `packaging/`, and we built and tested both locally.
- Library type: tchem typed reads that start with an adapter as a tag library. In ERR15993554,
  70% of reads start with the SMART oligo. tchem must skip known adapter starts when it counts
  read starts.
- Add chaos probes and fuzz seeds for each new input path or kit.

## Milestones

| Milestone | Contents | Exit criteria |
|---|---|---|
| v0.1.0 | today's layouts and record | Public repository, tag, release, packages submitted |
| v0.1.1 | PIPseq and inDrops recall, labelling rules | Both kits run again on held-out runs; a new coverage pass |
| v0.2 | index-read barcodes (ATAC), more held-out studies, `--emit` tests | scATAC and Multiome ATAC detected, with a held-out record; at least four more kits at `validated` |
| v0.3 | Visium, the CEL-seq family, ddSEQ | Detected, with a held-out record |

## Open decisions

- Spatial: we propose it for 1.x. Then 1.0 stays single-cell RNA, broad and proved.
- seqspec as the interchange format: we propose yes.
- The checks on every commit use five held-out runs (three Singleron, two GEM-X 10x demos). We
  can move them to design. Or we can change the rule of `validate_kits.tsv` (`bench/VALIDATION.md`,
  section 7).

## Working rules

These rules apply to people and agents who change the code:

- Work on a branch.
- Run `scripts/gate.sh` before every commit, and examine its exit code.
- Freeze a catalogue before you run on it.
- Do a coverage pass only with a frozen binary.
- A held-out run that you use to change the code becomes a design run.

CONTRIBUTING.md gives the other rules.

## Done

- Short reads, 10x: 3' v2, v3/v3.1, v4 (GEM-X); 5' v1/v2, v3 (GEM-X); Multiome GEX; Flex v1 (list
  from a Cell Ranger install).
- Short reads, BD Rhapsody: V1, Enhanced and Enhanced V2/V3 beads.
- Short reads, other vendors: SPLiT-seq; Parse Evercode WT and WT v3; inDrops; Drop-seq.
- ScaleBio RNA and QuantumScale (from R2).
- Singleron GEXSCOPE V1, V2, V3; MGI DNBelab C4/HT; SeekGene SeekOne DD 3' and 5'.
- Illumina Single Cell 3' / PIPseq v3-v4, and Illumina RNA Prep / PIPseq V.
- Long reads: the five 10x kits on Oxford Nanopore and PacBio, FASTQ or BAM.
- Library type (gene expression, VDJ, CRISPR guide, tag library).
- The 3'/5' strand without a reference.
- The files of a lane, identified by position, with index reads set aside.
- A held-out record of 10,130 public runs, recall for each kit, and second opinions on every
  disputed call.
- `--emit` commands for Cell Ranger, STARsolo, simpleaf and kb.
