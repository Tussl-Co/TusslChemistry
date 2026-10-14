# Contributing

The most useful contribution is a kit that tchem does not know yet. A wrong kit does not make a
pipeline fail. The pipeline finishes with fewer cells or the wrong genes. For this reason, every
kit goes through the same steps before tchem names it. The terms in this file have the meanings
in [docs/glossary.md](docs/glossary.md).

## Questions and concerns

- Ask a question or give an idea in
  [Discussions](https://github.com/Tussl-Co/TusslChemistry/discussions).
- Report a wrong answer or a bug as an issue.
- Report a security problem privately, as [SECURITY.md](SECURITY.md) tells you.

Be civil and specific. Everyone in the project follows the [code of conduct](CODE_OF_CONDUCT.md).
To report how a person behaves, send an email to opensource@tussl.co. Do not open an issue for it. If
something else in the project concerns you, open an issue or start a discussion. We answer it
there.

## Asking for a kit

When tchem finds a structured read that matches no known layout, its report says so. Do these
steps:

1. Run tchem again with `--format json`.
2. Copy the `"fingerprint"` object from the report.
3. Paste it into the [unknown-kit issue form](https://github.com/Tussl-Co/TusslChemistry/issues/new?template=unknown-kit.yml).

The fingerprint holds no reads:

- For each file, it records the number of reads in the sample and their lengths.
- For a file with the structure of a barcode read, it adds a sketch with one character for each
  position. The character is a consensus over all reads in the sample. It shows where the fixed
  sequence, the random sequence and the polyT are.
- It includes a 12-base sequence only when at least 15% of reads have it. Only a linker is that
  frequent. The transcript of one person is not.

The fingerprint does not contain file paths or reads. The other parts of the JSON report contain
file paths, so examine the report before you share it.

With the fingerprint, please tell us:

- the name and version of the kit, as your lab knows them,
- where the read structure of the kit is documented, if possible in the vendor's own pipeline or
  user guide,
- one or more public runs (SRA, ENA or GEO accessions) that use the kit.

## Adding a kit

Every kit follows one procedure: [the kit recipe in ROADMAP.md](ROADMAP.md#the-kit-recipe). The
held-out record of the kit follows [bench/VALIDATION.md](bench/VALIDATION.md), section 5.

## Issues

Each issue has a label:

- `kit request`: a kit that tchem does not know. The unknown-kit form adds this label.
- `wrong answer`: a wrong kit, or a confidence or refusal that seems wrong. The wrong-answer
  form adds this label.
- `good first issue`: a small change with a clear result.
- `help wanted`: a change that we want a contributor to do.

## Taking an issue

1. Find an issue that has no assignee.
2. Write a comment on the issue to say that you will do it.
3. A maintainer assigns the issue to you.
4. For more than a small fix, agree the method on the issue before you write code.

A new kit always follows the kit recipe. If the issue has no pull request or update in 14 days, it
becomes available again. You can ask for more time on the issue.

## Pull requests

- Make one change in each pull request.
- CI runs `cargo fmt --check`, clippy and the tests on each pull request. Your change must pass
  them.
- `scripts/gate.sh` needs the benchmark data ([Reproducing the checks](#reproducing-the-checks)).
  A maintainer runs it before the merge.
- If a change alters a report, the commit message must give the reason. The parity oracle shows
  each altered report.
- A change to detection needs a new coverage pass. A maintainer runs it.
- Write the commit message as `area: what is now true`. An example is `emit: STARsolo reads
  gzipped files with gunzip -c`.
- You submit your contribution under the Apache License 2.0, section 5. There is no contributor
  licence agreement.

## The bar for code

Every change passes the gate, `scripts/gate.sh`. The gate runs these checks:

- clippy, with the lint table of the crate,
- the tests,
- the parity oracle: the output must not change, unless the commit says why it changes,
- the JSON schema,
- the validation tables,
- `scripts/validate.sh`,
- the chaos cases and the probes: damaged files and hostile names and paths. tchem must not abort
  or hang, and no HIGH call can be wrong.

Input from outside (files, the command line) gets a typed refusal, never a panic.

## Engineering rules

The code follows [TigerStyle](https://tigerstyle.dev):

- Functions have 70 lines or fewer. Lines have 100 columns or fewer.
- Every input has a limit. The passes over each read do not allocate memory. Allocation-count
  tests prove this.
- Untrusted input (FASTQ, barcode lists, packed lists, BAM, the command line) gets a typed
  refusal, not a panic.
- There is one dependency: `flate2`, for gzip.

The assertion rule, as `tests/tidy.rs` applies it in `cargo test`:

- Each function must have at least two assertions that run in a release build. The test counts
  `assert!`, `assert_eq!`, `assert_ne!` and `unreachable!`. It does not count the `debug_` forms
  or a statement under `#[cfg(debug_assertions)]`.
- 60 functions are exempt (`EXEMPT`), each with a reason. Most are lookups, display code or
  validators whose refusal is the check.
- 16 functions on the per-read path have their checks as `debug_assert` (`DEBUG_ONLY`). Each must
  have at least two checks in the debug build. The gate runs the debug build on all 281 parity
  oracle inputs.
- The test fails if an entry in either list is no longer necessary.
- Release builds keep `assert!` and `assert_eq!` on. A failed assertion stops the program.
- `src/emit/simpleaf.rs` tests that the geometry text of each layout agrees with its segments.
- The debug build is an optimised build with overflow checks and more cross-checks.

The parity oracle has 281 recorded inputs. Their text and JSON reports must stay the same, byte
for byte, in release and in debug builds. A change to a report must be intentional and reviewed.

The fuzzing (`fuzz/`) is coverage-guided. Its seeds are real reads of each layout that tchem
calls. On 10 Oct 2026, it gave these inputs, with the real barcode lists:

- 44 million inputs to the readers (FASTQ, gzip, BGZF, BAM, SRA),
- 4.9 million paired inputs through detection and every output.

It found no crash, no hang and no out-of-memory failure.

## Reproducing the checks

The benchmark scripts read their inputs from one directory, `TUSSL_CHEMISTRY_DATA`. This
directory holds public 10x, BD, Parse, SPLiT-seq and Drop-seq data. It also holds the ENA and SRA
downloads that the scripts fetch. Recorded result tables show this directory as
`$TUSSL_CHEMISTRY_DATA/...`.

```sh
export TUSSL_CHEMISTRY_DATA=/path/to/benchmark-data    # barcode lists in $TUSSL_CHEMISTRY_DATA/whitelists
scripts/gate.sh          # clippy, tests (with tidy), parity oracle in release and debug
scripts/validate.sh      # known runs and damaged copies, verdict by verdict
python3 bench/chaos/scripts/chaos.py target/release/tchem $TUSSL_CHEMISTRY_DATA/whitelists work/ chaos.tsv
```

Each benchmark under `bench/` keeps its code in `scripts/` and its recorded data in `results/`.

- `bench/head_to_head` holds the frozen selection (`select.py`, `selection.tsv`) and the runner.
  It also holds the adjudication evidence and the result of each experiment.
- `bench/long_read` holds the long-read catalogues (design and held-out) and the byte-range head
  fetcher. It also holds the runner and the scorer.
