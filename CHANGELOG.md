# Changelog

The JSON report has its own version (`schema_version`, see [docs/json.md](docs/json.md)). This
file follows the version of the tool.

## 0.1.0 (2026-10-14)

The first public release.

**Detection**

- tchem names the layout of a single-cell RNA library from its first reads, for all supported
  vendors. These are the kits:
  - 10x 3' v2/v3/v4, 5' v1/v2/v3, Multiome GEX and Flex v1;
  - BD Rhapsody (original, Enhanced and Enhanced V2/V3 beads);
  - Parse Evercode WT (both layouts) and SPLiT-seq;
  - MGI DNBelab C4 / C-series HT;
  - Illumina Single Cell 3' / PIPseq v3-v4 and Illumina Single Cell 3' RNA Prep (PIPseq V);
  - SeekGene SeekOne DD 3' and 5';
  - ScaleBio RNA and QuantumScale;
  - Singleron GEXSCOPE V1-V3, inDrops and Drop-seq.
- The report gives the geometry and the barcode lists. It also gives which file is which, even
  when the files are swapped.
- The report gives the cDNA strand without a reference. It also gives the library type: gene
  expression, VDJ, CRISPR guides, or a tag library (for example, antibody or hashtag tags).
- tchem measures the strand and the library type on the cDNA read, in the file where the kit puts
  it.
- tchem reads long reads (Oxford Nanopore, PacBio; FASTQ or BAM) for the five 10x
  gene-expression kits. A long-read run in two files (`_1` and `_2`) goes to the long-read
  detector, which reads the first file.
- A barcode list counts as evidence only when its hits are spread over many barcodes. Thus a
  constant read (the amplicon of a tag library) near a list entry cannot support a kit.

**Trust in each answer**

- Every match gives the held-out record of its layout (`layout.validation`). The record counts
  the HIGH calls on held-out runs, and the studies of those runs.
- The record also gives a level, from the number of studies and the rule of three. It also counts
  the calls right only for the kit family, and the disputed calls.
- `bench/coverage` measures the record on the binary that we ship. `scripts/release_check.sh`
  refuses a release with an old record.
- On the held-out runs: 8,923 HIGH calls on 10,130 public runs from 5,212 studies. None is wrong.
  999 are disputed, each with its evidence.
- A level counts only independent labels: the protocol of the run, the lab's own record, the
  vendor, a BD part number. The README shows calls with only a read-derived label, but a level
  does not count them (README, "How well it works"; bench/VALIDATION.md).
- The README section "How well it works" tells what each level means when you act on a HIGH
  call. It also shows recall. Recall is how frequently tchem calls each labelled kit HIGH, below
  HIGH, ambiguous or not at all.
- The recall section also gives the known gaps: no call on some PIPseq and inDrops runs, and
  Multiome labels on ATAC libraries. bench/coverage/scripts/recall.py makes the recall table, and
  the gate compares it.

**Acting on the answer**

- `--emit cellranger|starsolo|simpleaf|kb` prints a ready command. The reference is a named
  placeholder.
- `--emit` refuses (exit 3) when a pipeline would read the files wrongly. This includes paths
  that a pipeline would split or read as an option, pipes, BAM files and SRA files.
- For Multiome (nuclei), the STARsolo and kb commands count intronic reads, as Cell Ranger does.
- The STARsolo command reads gzipped files with `gunzip -c`, so it runs on Linux and macOS.
- `bench/emit` runs the STARsolo, kb and simpleaf commands to the end on real reads. STARsolo has
  a strand control.
- Only the argument parser of Cell Ranger tested its commands. No test finds a wrong UMI length
  yet. docs/emit.md tells what each test covers.
- `tchem SRR12345678` reads a public run by accession. It streams the first reads of each file
  from ENA through curl, and does not download the full run.
- curl uses HTTPS only, and ignores the `~/.curlrc` of the user.

**Asking for a kit**

- A no-match tells why. Either the reads contradict every candidate and no vendor list supports
  one, or the best support is less than 30%. Only the second suggests a kit that tchem does not
  support yet.
- In the first case, the report points to the contradictions of each candidate. The reads of a
  variant of a known kit, or trimmed reads, look the same. The text and the JSON `hint` both say
  this.
- A no-match with the structure of an unknown kit prints a link. The link opens the unknown-kit
  issue form, with the `fingerprint` of the run filled in.
- The fingerprint holds no reads. For a sample of fewer than 1,000 reads, it leaves out the
  consensus for each position and the shared k-mer. Over so few reads, the consensus would be
  the reads themselves. CONTRIBUTING.md explains this.

**Instrument**

- Every report names the instrument (`instrument`). The answer comes from the name of the first
  read: the Illumina model by instrument prefix, Element, MGI DNBSEQ, Ultima, PacBio or Oxford
  Nanopore.
- For BAM, the answer comes from the platform tags of the read group. For a run read by
  accession, the report shows the model from ENA next to it. When the names are stripped, the
  report says the instrument is unknown.

**Robustness**

- Untrusted input gets a refusal, never an abort. `bench/chaos` covers read names with control
  characters, IUPAC bases, hostile paths and pipes. It has 425 damaged files and 50 probes.
- `scripts/gate.sh` runs `bench/chaos` on every change.
- Fuzzing (`fuzz/`) ran on 10 Oct 2026, 2 hours for each target. It gave 44 million inputs to
  the readers and 4.9 million paired inputs through detection and every output.
- The fuzzing found no crash, no hang and no out-of-memory failure.

**Contract**

- The JSON report has a version (`schema_version` 1.1). Version 1.1 adds
  `layout.family_validation`. [schema/report.v1.json](schema/report.v1.json) describes the report.
- `--emit --format json` has its own schema, [schema/emit.v1.json](schema/emit.v1.json).
- Exit codes: 0 the reads were read (any verdict), 1 input refused, 2 command line refused, 3 no
  command fits.
- The Rust library stays at 0.x until its API is stable. Pipelines should build on the JSON
  contract.

**Documentation**

- The documentation follows Simplified Technical English rules, with one meaning for each term
  ([docs/glossary.md](docs/glossary.md)). `scripts/ste_check.py` checks the rules on every
  commit.
