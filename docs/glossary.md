# Glossary

Each term in this file has one meaning. The README, CONTRIBUTING.md, ROADMAP.md, CHANGELOG.md,
SECURITY.md, `bench/VALIDATION.md` and the files in `docs/` use each term only with that meaning.
The text in these files follows Simplified Technical English rules. `scripts/ste_check.py` checks
these rules on every commit:

- A sentence that describes has 25 words or fewer. A step in a procedure has 20 words or fewer.
- A paragraph has six sentences or fewer.
- A word that this glossary replaces does not occur. For example, write "barcode list", not
  "whitelist".

The check does not test grammar, voice or the meaning of a word. A reviewer checks those.

Terms in `code format` are names that tchem prints, such as a level or a JSON key.

## The tool

| Term | Meaning |
|---|---|
| tussl-chemistry | The project: the repository, the Rust crate (`tussl_chemistry` in code) and this documentation. |
| tchem | The command that tussl-chemistry installs. |
| report | What tchem prints for one set of files, as text or as JSON (`--format json`). |
| gate | The checks that every commit must pass: `scripts/gate.sh`. |
| parity oracle | 281 recorded inputs. The gate runs tchem on each and compares the report, byte for byte, with the recorded report. It does this for text and for JSON. |
| chaos case | A real run that the gate damages in one of 17 ways before it runs tchem on it. There are 425 chaos cases. |
| chaos base | The real run of one layout that the chaos cases damage. Each layout has one. |
| probe | An input with a hostile file name, path or base. There are 50 probes. |
| regression input | A run with a known kit that `scripts/validate.sh` gives to tchem on each commit. |
| fuzzing | A test that gives tchem's readers and detection millions of generated inputs (`fuzz/`). |

## Reads and files

| Term | Meaning |
|---|---|
| read | One sequence from a sequencing instrument, with its name and qualities. |
| barcode read | The read that holds the cell barcode and the UMI. On Illumina this is usually R1. |
| cDNA read | The read that holds the transcript sequence. On Illumina this is usually R2. |
| index read | A read (I1 or I2) that holds the sample index. tchem sets index reads aside. |
| short reads | Reads from Illumina, MGI, Element or Ultima, with the barcode at a fixed position. |
| long reads | Reads from Oxford Nanopore or PacBio. One read holds the full molecule, in either direction. |
| swapped | The barcode read is in the file named R2, and the cDNA read is in the file named R1. |
| sample | The reads that tchem reads from the start of each file. The default is 100,000 reads. A "GEO sample" is a different thing: a GEO record with a GSM accession. |
| read type | Short reads or long reads (`read_kind` in JSON). A level applies to one layout and one read type. |

## Kits and layouts

| Term | Meaning |
|---|---|
| kit | A product that a vendor sells to make single-cell libraries, such as 10x Chromium 3' v3.1. |
| layout | The position of the barcode, the UMI, the linkers and the cDNA in the reads of one kit. tchem knows 25 layouts. On long reads, it detects the five 10x gene-expression layouts. One layout can cover more than one kit, such as "10x 3' v3/v3.1". |
| kit family | Two or more layouts of one vendor product line that labs usually do not tell apart: BD Rhapsody, Parse Evercode WT, SeekGene SeekOne DD and Singleron GEXSCOPE. |
| barcode list | The vendor's list of all cell barcodes for a kit. tchem does not include barcode lists. `scripts/fetch_whitelists.sh` downloads them. |
| listed barcode | A barcode on the barcode list, exactly or with one different base. |
| unlisted barcode | A barcode of a kit that has no public barcode list. |
| linker | A fixed sequence between the parts of a barcode. |
| geometry | The text that simpleaf and alevin-fry use to describe a layout, such as `1{b[16]u[12]x:}2{r:}`. |
| strand | The direction of the cDNA read relative to the transcript: forward for 3' kits and reverse for 5' kits. |
| library type | Gene expression, immune receptor (VDJ), CRISPR guide, or a tag library (antibody, hashtag or amplicon). |
| chemistry | Use this word only for Cell Ranger's `--chemistry` option. Elsewhere, write "kit" or "layout". |

## What tchem says

| Term | Meaning |
|---|---|
| verdict | The answer of a report: `match`, `ambiguous` or `no_match`. |
| match | The verdict when one layout fits the reads. A match has a confidence. |
| ambiguous | The verdict when two or more layouts fit the reads equally well, and a pipeline must process them differently. |
| no match | The verdict when no layout fits the reads (`no_match` in JSON). The text report prints it as "no known layout". |
| confidence | HIGH, MEDIUM or LOW. Only a match has a confidence. |
| HIGH | A match with at least 70% support, no contradiction, and evidence from a barcode list or a linker. |
| MEDIUM | A match with a clear lead over the other layouts, but less than HIGH needs. |
| LOW | A match where the barcode list agrees with the layout, but the reads contradict another part of it. LOW is not a weaker MEDIUM. |
| call | The layout that a match names, with its confidence. A "HIGH call" is a match at HIGH. |
| support | The mean share of reads that show a layout's barcode, linker and unlisted-barcode evidence. |
| contradiction | A thing in the reads that a layout does not permit, such as a missing linker or a read that is too short. |
| candidate | A layout that the report shows with its evidence. An ambiguous verdict names its candidates. |
| no call | A run where tchem gives no match: the verdict is `no_match`, or tchem refuses the input. |
| hint | The lines under a no-match verdict that tell the most probable cause and the next step. |
| fingerprint | What each file shows when the verdict is no match, without reads. It goes with a request for a new kit. |

## The held-out record

| Term | Meaning |
|---|---|
| run | One sequencing run in a public archive, with an accession that starts SRR, ERR or DRR. A vendor's demonstration run is also a run. |
| experiment | The archive record for one library, with an accession that starts SRX, ERX or DRX. An experiment has one or more runs. |
| study | The project that a run belongs to, with an accession that starts PRJNA, PRJEB or PRJDB. All demonstration runs of one vendor are one study. |
| label | What a source other than the reads says the kit of a run is. |
| independent label | A label from a source that did not read the reads: the run's protocol text, the lab's own record, the vendor's statement, or a BD part number. |
| read-derived label | A label from a source that read the reads, such as Arc Institute's barcode-list matching. A level does not count it. |
| exact label | A label that names one layout. |
| family label | A label that names a kit family but not the layout in it. |
| design run | A run whose reads we used to make or change a layout. A held-out run that we use to change the code becomes a design run (the design-data rule). |
| design study | A study that contains a design run. |
| examined study | A study that contains a run whose reads we looked at during research. |
| held-out run | A run that is not a design run and is not in a design study or an examined study. Only held-out runs go into the held-out record. |
| held-out record | The calls of one frozen tchem binary on all held-out runs, with their labels (`bench/coverage`). |
| catalogue | The list of runs that one benchmark selects from archive metadata, with the label and the split (design or held-out) of each run. We commit it before tchem reads any of the runs. |
| draw | One selection of runs from archive metadata into a catalogue, by rules that we commit first. |
| read head | The first part of a run's files, which a benchmark downloads and keeps. |
| coverage pass | One pass of a frozen tchem binary over all held-out runs (`bench/coverage/scripts/run.py`). |
| referee | The source of a run's label, such as the protocol text or the lab's record (the `referee` column of `calls.tsv`). |
| precision | The share of HIGH calls that are not wrong. A level measures precision. |
| right | A HIGH call that agrees with the run's independent label. |
| wrong | A HIGH call that the label does not permit, and that adjudication does not make right or disputed. |
| disputed | A HIGH call whose label names a different kit, where the reads' evidence supports the call. The evidence is the barcode list or, where two kits share a list, the strand. A disputed call is not right and not wrong. |
| adjudication | The decision, by a rule written before the coverage pass, that a call that disagrees with its label is right, disputed or wrong (`bench/coverage/results/adjudication.tsv`). |
| second opinion | A different method that reads the files of a disputed call: Sanger's `infer_10x_chemistry.py`, a strand referee, or Cell Ranger (`bench/second_opinion`). |
| strand referee | A test of the cDNA strand on 40 genes that tchem does not use. |
| recall | The share of runs with a label that tchem calls, by outcome. A level does not use recall. |

## Levels

| Term | Meaning |
|---|---|
| level | The validation level of one layout and one read type, from its held-out record. Every match prints the level of its layout. |
| `unvalidated` | No held-out HIGH call with an independent exact label. |
| `checked` | HIGH calls in 1 to 29 held-out studies, none wrong. The report gives the counts and the error bound, but the bound is too wide to claim a rate. |
| `validated` | HIGH calls in 30 to 99 held-out studies, none wrong. The error rate of a HIGH call is less than 10%, at 95% confidence. |
| `error_under_3pct` | HIGH calls in 100 or more held-out studies, none wrong. The error rate of a HIGH call is less than 3%, at 95% confidence. |
| `contradicted` | At least one held-out HIGH call is wrong. |
| family level | The level of a kit family. A call that names the wrong layout in the correct family is right for the family level. |
| rule of three | If n independent trials have no error, the true error rate is less than 3/n at 95% confidence. A level counts studies as trials, because runs from one study are not independent. |
| error bound | 3 divided by the number of studies, from the rule of three. |

Do not use "checked" or "validated" as ordinary verbs. Write "tested", "verified" or "compared".

## The rules of the held-out record

`bench/VALIDATION.md` gives the full procedure. Three of its rules have names:

| Term | Meaning |
|---|---|
| the design-study rule | A design study cannot give a held-out run for the layout that its design run changed. Earlier records call this "rule 2". |
| the independent-label rule | A level counts only independent labels. Earlier records call this "A4". |
| the design-data rule | A held-out run that we use to change the code becomes a design run, before the next coverage pass. Earlier records call this "rule 7". Its list is `bench/coverage/results/rule7.tsv`. |
