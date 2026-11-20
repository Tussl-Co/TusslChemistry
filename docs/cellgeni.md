# The gaps that cellgeni found, and what we changed

The Cellular Genetics team of the Wellcome Sanger Institute (cellgeni) measured tchem 0.1.0 on
1,227 runs of their pipeline, on 19 Oct 2026. Their report is `docs/tchem-benchmark.md` in
[cellgeni/nf-reprocessing-public-10x](https://github.com/cellgeni/nf-reprocessing-public-10x).
This page gives each gap that they found, what we changed, and how we tested the change.

**Status:** the changes are in development. They are not in a release, and they are not in this
repository yet. A change to detection makes the held-out record old. Thus the changes wait for
version 0.1.1 and a new coverage pass. Until then, the record in the README is the record of 0.1.0.

## Each gap

| Gap that cellgeni found | What we changed | Status |
|---|---|---|
| 5' feature libraries typed as gene expression (0 of 23 typed as tag) | The tag test counts read starts at base 1 and at base 11 | Changed |
| A tag library's polyA runs called 3' v2 | A tag or guide library's strand does not decide 3' or 5'. The ambiguity gives this reason | Changed |
| 24 VDJ runs ambiguous between 3' v2 and 5' v1/v2 | A VDJ library decides 5'. tchem types VDJ before the tag test | Changed |
| HIGH calls on files of fewer than 1,000 reads | A HIGH call needs at least 1,000 reads. Below that, the call is MEDIUM, with a note | Changed |
| Index reads of 16 bases refused | tchem sets aside index reads of up to 19 bases | Changed |
| ENA's unpaired third file refused in local files | tchem sets aside an unnumbered file beside its own `_1` and `_2` | Changed |
| Every barcode list needed under tchem's own file names | `scripts/fetch_whitelists.sh --from <dir>` takes the lists from a directory that a pipeline already has | Changed |
| 10x 3' v1, ATAC and bzip2 input | None | Not supported |

### 5' feature libraries

In a 5' antibody, CITE-seq or hashtag library, the tag starts at base 11 of the cDNA read. In a
3' library, it starts at base 1. tchem 0.1.0 looked only at base 1. Now the tag test counts the
read starts at both positions and uses the higher count.

### Strand of a tag library

The reads of a tag or guide library do not come from a transcript. Their polyA runs and their
strand do not show 3' or 5'. Now tchem does not use them to break a tie between 3' v2 and 5' v1/v2.
The result is an ambiguity, with the reason "the cDNA read is a tag or guide library, not a
transcript". Their barcode list is the same, so the barcode list in the report stays correct.

### VDJ on the 737K list

10x makes VDJ kits only as 5' kits. Now a VDJ library breaks a tie between 3' v2 and 5' v1/v2 in
favour of 5' v1/v2, at HIGH. tchem now types a library as VDJ before the tag test. In a VDJ
library, many reads start with the same constant-region bases, and the tag test could count them.

### Small files

A HIGH call now needs at least 1,000 reads in each file. A call on fewer reads is MEDIUM at most.
The report gives a note, for example "only 999 reads in a file: a HIGH call needs at least 1000".

This change does not stop all calls on very small files. On a file of 4 reads, tchem can still
make a LOW call, and the call can be wrong. `--emit` refuses any call that is not HIGH and exits
with 3, so a pipeline that uses `--emit` stops on these files (see [emit.md](emit.md)).

### Inputs

tchem now sets aside index reads of up to 19 bases, so it accepts 10x dual-index runs with 16-base
I2 reads. It also sets aside ENA's unpaired third file in local files, as it already did for a run
that it reads by accession. The rule is strict: there must be exactly three files, and the
unnumbered file must have the same stem as `_1` and `_2`. The report says which files it set aside.

### Barcode lists

A pipeline that already holds the 10x lists does not need a second copy. Run
`scripts/fetch_whitelists.sh --from <dir>`. The script finds each list in that directory by the
sha256 of its lines, under any file name. It makes a link to a plain file. It makes a copy of a
gzip file or a file with Windows line ends.

A missing list still stops tchem. We did not change this. With a list missing, a kit that has no
list (Drop-seq, for example) could match a run and give a wrong call.

## How we tested the changes

| Claim | Check | Result |
|---|---|---|
| The changes correct the runs that cellgeni named | 0.1.0 and the new build on the same reads | Correct. Details below |
| The changes do not make other calls worse | 0.1.0 and the new build on 82 held-out runs, a seeded draw across all short-read layouts | 73 identical. The 9 changes are all corrections by their labels |
| Each new test can fail | We put each bug back and ran its test | Each test failed |
| Hostile input stays safe | The 50 chaos probes | 50 of 50 pass, with the same table as 0.1.0 |
| The record of each layout stays true | The coverage pass on the held-out runs | Not done. It comes with 0.1.1 |

The runs that cellgeni named:

- **SRR27765505** (5' CITE-seq): 0.1.0 called 3' v2 at MEDIUM, typed as gene expression. The new
  build gives an ambiguity between 3' v2 and 5' v1/v2, typed as a tag library. The tag start
  cover went from 0.055 to 0.953.
- **SRR27765486** (hashtag): the new build types it as a tag library. The call stays an ambiguity
  between 3' v2 and 5' v1/v2.
- **SRR22190228** (VDJ, PRJNA898288): 0.1.0 gave an ambiguity. The new build calls 5' v1/v2 at
  HIGH.
- **A local run with ENA's unpaired file and 16-base I2 reads:** 0.1.0 refused it. The new build
  sets aside both files and calls HIGH.
- **The same reads, cut to 999 and to 1,000 reads:** the new build calls MEDIUM with the note,
  then HIGH.

In the seeded draw, 4 TCR runs went from an ambiguity to 5' v1/v2 at HIGH. 5 libraries with HTO,
ADT, CITE-seq, hashtag or BD tag labels are now typed as tag libraries. We read the reads of the BD
run (SRR25629551) and found a tag library. No change replaced one kit or barcode list with
another. One run of the 83 failed to download from ENA.

## The design-data rule

We used some of cellgeni's runs to make and test these changes. These runs become design runs, and
the held-out record does not count them. They are SRR27765505, the VDJ runs of PRJNA898288 and
the runs of GSE274751.

## What we did not change

- **10x 3' v1, ATAC and bzip2 input:** not supported.
- **Drop-seq at MEDIUM on 20-base CRISPR runs:** we could not find the runs, so we could not
  examine them.
- **Short read 1 and the 3'/5' ambiguity:** sometimes read 1 ends at the UMI and the cDNA read
  shows no strand. Then tchem still cannot tell 3' v2 from 5' v1/v2. We do not know a fix that we
  can prove.
