# Fuzzing

Coverage-guided fuzzing (libFuzzer, through `cargo fuzz`) of what tchem reads and decides. It is not
part of the gate. Run it before a release, or after changing a reader or the detection code.

| Target | What it is given | What must hold |
|---|---|---|
| `reads` | any bytes, as the file or stream tchem is given: FASTQ (plain, gzip, BGZF), BAM or an SRA archive, sniffed from the first bytes | a typed refusal or reads, never a panic |
| `detect` | any bytes, split at the first 0xFF into R1 and R2, read as above and run through `detect`, the text report, the JSON report and every `--emit` pipeline, with the real barcode lists | a report, never a panic or a broken assertion |

Release builds keep their assertions on, so a broken invariant in detection is a crash here.

## Running

```sh
cargo install cargo-fuzz          # needs a nightly toolchain: rustup toolchain install nightly
TUSSL_CHEMISTRY_DATA=<data> python3 -I fuzz/scripts/seed_corpus.py <data> fuzz/corpus
cd fuzz
cargo +nightly fuzz run -O reads corpus/reads -- -max_total_time=7200 -jobs=6 -workers=6
TCHEM_FUZZ_WHITELISTS=<whitelist dir> cargo +nightly fuzz run -O detect corpus/detect -- \
  -max_len=400000 -max_total_time=7200 -jobs=6 -workers=6 -timeout=60 -rss_limit_mb=6144
```

The seeds are held-out read heads: up to two runs for each layout tchem called HIGH, as plain FASTQ,
gzip, BGZF, BAM and R1/R2 pairs (`scripts/seed_corpus.py`). A crash is saved under `artifacts/` as
the input that reproduces it: `cargo +nightly fuzz run <target> <file>`.

## Record

| Date | Commit | Target | Time | Jobs | Inputs run | Crashes, hangs, out of memory |
|---|---|---|---|---|---|---|
| 10 Oct 2026 | 28aad54 | `reads` | 2 h | 6 | 44,455,108 | none |
| 10 Oct 2026 | 28aad54 | `detect` | 2 h | 6 | 4,893,966 | none |

Seeds: 40 runs from 21 layouts. The corpora grew to 26,462 (`reads`) and 12,093 (`detect`) inputs.
Coverage reached 1,387 edges for `reads` and 5,954 for `detect`.
