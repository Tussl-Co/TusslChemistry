# The JSON report

`tchem --format json` prints one JSON object per run. The object is the output contract that
pipelines can build on. This page explains it, and
[`schema/report.v1.json`](../schema/report.v1.json) is its machine-readable form. The gate
checks every recorded report against that schema, so the page, the schema and the binary agree.

## Versioning

`schema_version` is the contract's version. It is separate from `version`, which is the tool's
version. Check `schema_version`, not `version`.

| Change | Version step | Example |
|---|---|---|
| A key is added | minor: `1.0` to `1.1` | a new evidence share |
| A key is renamed, removed or changes type | major: `1.x` to `2.0` | `files` becomes an array |

Rules for consumers:

- Ignore keys that you do not know. A minor version can add them.
- Refuse a major version that you were not written for.

Rules for the tool:

- Every key is always present. A key that does not apply is `null`, `[]` or `false`. It is
  never left out.
- A major version gets a new schema file (`report.v2.json`). The old file stays, so pipelines
  pinned to the old version can still validate.

## Exit codes

| Code | Meaning | stdout |
|---|---|---|
| 0 | The reads were read. The verdict can be `match`, `ambiguous` or `no_match`. | the report |
| 1 | An input was refused: an unreadable or malformed file, or a missing list directory | empty; the reason is on stderr |
| 2 | The command line was refused | empty; the reason is on stderr |
| 3 | `--emit` found no command that fits (see [emit.md](emit.md)) | empty; the reason is on stderr |

`no_match` is a result, not an error, so it exits 0.

## Keys

| Key | Type | Meaning |
|---|---|---|
| `tool` | `"tussl-chemistry"` | Always this value. |
| `version` | string | The tool's version. |
| `schema_version` | string `1.N` | The contract's version (see above). |
| `verdict` | `match`, `ambiguous`, `no_match` | The answer. |
| `confidence` | `high`, `medium`, `low` or null | Set only for a match. |
| `layout` | object or null | Set only for a match: `name`, `geometry` (simpleaf / piscem / alevin-fry), `barcode_lists` (paths), `note` (what processing must know beyond the geometry), `long_read_kit` (Oxford Nanopore `wf-single-cell --kit`), `validation` (see below). |
| `candidates` | array of layout names | Set only when the verdict is ambiguous. These are the layouts that the reads cannot separate. |
| `files` | object or null | `barcode_read` and `cdna_read`: which input file is which, after swap detection. For a run read by accession, these are the files' HTTPS URLs on ENA. For long reads, both name the one file that was read. A long-read run in two files is read from the first: the second is set aside and not named here (the text report names it). |
| `cdna_strand` | `forward`, `reverse` or null | The strand to give the mapper: `forward` for 3' kits, `reverse` for 5' kits. |
| `library` | object or null | The library type (`kind`) and the shares it was judged on. |
| `swapped` | boolean | True when the file labelled R1 holds what the layout puts on R2. |
| `strand_evidence` | object or null | The reference-free strand signals, as shares of reads. |
| `contradictions` | array of strings | Evidence against the chosen layout. |
| `notes` | array of strings | Facts about the match that do not count against it. |
| `hint` | string or null | Set only when the verdict is no_match: why nothing matched. |
| `not_checked` | array of layout names | Layouts that were not scored because their barcode list was not fetched. |
| `reads_sampled` | integer | The number of reads that were read. |
| `read_kind` | `short` or `long` | Which detector answered. |
| `instrument` | object | Which instrument read the files (see below). |
| `long_read_evidence` | object or null | Set only for long reads: adapter, barcode list, tail and UMI evidence. |
| `fingerprint` | object or null | Set only when the verdict is no_match: what each file shows, to share when asking for a kit (see below). |

All shares are fractions of reads from 0 to 1.

## `layout.validation`

This object gives the held-out record of the matched layout. A held-out run is a public run
that we did not use to make or change the layout. The record describes the layout for this read
kind, not the run that you gave. It applies to HIGH calls only, as `applies_to` says.

| Key | Meaning |
|---|---|
| `level` | `unvalidated`, `checked`, `validated`, `error_under_3pct` or `contradicted` |
| `held_out_high_calls` | HIGH calls on held-out runs whose own record names exactly this layout |
| `studies` | The studies those runs come from |
| `family_only_calls` | HIGH calls that are right only at the kit-family level, because the run's record does not name the version (BD bead versions, Parse kit versions) |
| `disputed` | HIGH calls whose label names a different kit, where the reads' evidence supports the call: the barcode list or, where two kits share a list, the strand |
| `wrong` | Held-out HIGH calls that the run's own record contradicts |
| `error_bound_95` | 3 / `studies`: the rule-of-three upper bound on the error rate of a HIGH call, at 95% confidence. It is null when no call has an exact independent label, or when a call is wrong. |

The level counts studies, not runs. Runs from one study share a lab, a protocol and often a
library. Thus they are not independent trials, and the rule of three needs independent trials.

| Level | Studies, with no wrong call |
|---|---|
| `unvalidated` | 0 |
| `checked` | 1 to 29 (counts are stated, but no error rate is claimed) |
| `validated` | 30 to 99 (error rate under 10%) |
| `error_under_3pct` | 100 or more (error rate under 3%) |
| `contradicted` | any wrong call, whatever the count |

`held_out_high_calls` and `studies` count only independent labels: the protocol of the run, the
lab's own record, the vendor's statement, or a BD part number. They do not count calls with only
a read-derived label (bench/VALIDATION.md, section 4).

## `layout.family_validation`

Added in schema 1.1. Labels of some kit families rarely name the version or bead type:
`BD Rhapsody`, `Parse Evercode WT`, `SeekGene SeekOne DD` and `Singleron GEXSCOPE`. For a layout in
such a family, this object gives the record of the family. It counts the held-out HIGH calls that
name a layout of the family, and how many of those runs are that family.

It has the keys of `validation`, without `family_only_calls`, and with `family`, the name of the
family. It is null for every other layout. A call that names the wrong layout in the correct
family is right here. `validation` counts the same call as disputed or wrong.

The numbers come from `data/validation.tsv` and `data/validation_family.tsv`, which `bench/coverage` writes. It runs the binary
over every held-out run that the benchmarks have and then applies the adjudications in
`bench/coverage/results/adjudication.tsv`. Each adjudication is a confident disagreement with
a run's label, decided once and recorded with its evidence.

## `fingerprint`

When nothing matches, `fingerprint` describes each file so that someone can ask for the kit.
It holds no reads. The object is `{"fingerprint_version": 1, "files": [...]}`, and each file has
these keys:

| Key | Meaning |
|---|---|
| `file` | `R1` or `R2` |
| `reads`, `length_min`, `length_mode`, `length_max` | how many reads were sampled, and their lengths |
| `structured` | whether the file shows a barcode read's positional structure |
| `fixed_columns` | how many positions one base dominates |
| `sketch` | A consensus over every read sampled, one character per position. `.` means random, an upper-case base means one base holds over 70%, a lower-case base means it leads, and `-` means mixed. It is set only for structured files sampled with at least 1,000 reads (over fewer, a consensus is the reads themselves); a cDNA file gets null. |
| `anchor` | `{kmer, share}`: the commonest 12-mer near the read start, set only when at least 15% of at least 1,000 reads carry it. A sequence that common is a linker, not a transcript. |

[CONTRIBUTING.md](../CONTRIBUTING.md) explains how to share it.

## `instrument`

This object records which instrument read the files. It comes from the read names or from a BAM
read group. Detection never uses it.

| Key | Meaning |
|---|---|
| `platform` | `illumina`, `element`, `mgi`, `ultima`, `pacbio`, `oxford_nanopore`, or null when nothing in the source says |
| `model` | The model, when the instrument prefix of the read name or the BAM `PM` tag fixes it. For example, `A…` is NovaSeq 6000, `LH…` is NovaSeq X, `m84…` is Revio. |
| `id` | The instrument, flowcell or read ID that the name carries |
| `from` | `read_name` or `read_group` (a BAM `@RG` line's `PL`/`PM` tags) |
| `archive_model` | For a run read by accession: the `instrument_model` that ENA records |

Archives often strip read names down to `SRR123.1`. In that case `platform` is null, which is
"unknown" and not a guess. ENA's model and the read names sometimes disagree. Some runs labelled
"Element AVITI" have NextSeq names, and some runs labelled "UG 100" have NovaSeq names. For this
reason, the report gives both, side by side.
