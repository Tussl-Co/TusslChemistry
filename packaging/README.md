# Packaging

`bioconda/` is the recipe for
[bioconda-recipes](https://github.com/bioconda/bioconda-recipes) (`recipes/tussl-chemistry/`).
`nf-core/modules/tchem/` holds two modules for
[nf-core/modules](https://github.com/nf-core/modules) (`modules/nf-core/tchem/`):
`TCHEM_FETCHLISTS` and `TCHEM_DETECT`. Both were built and tested here, locally. Nothing was
published: the checks below were made on local builds. The steps to publish are in
"Steps to publish" below, with a box for each one.

## What the package holds

The package installs `tchem` and `tchem-fetch-whitelists`. The fetcher reads its registry from
`share/tussl-chemistry/data/lists.tsv`. The barcode lists belong to their vendors, so they are
never shipped: the fetcher downloads each list and refuses any list whose barcode count, length
or sha256 differs from the registry. The recipe's tests check that `data/` holds only
`lists.tsv`.

## What was built and tested here (2026-10-08, source at 56a3eac)

### Release workflow (checked earlier)

- The release workflow's steps were run locally: the tag and changelog check, the release notes
  taken from the changelog, and the packaging.
- The static musl build (stripped, 2.1 MB) passes both parity oracles byte for byte. The aarch64
  cross-build links here, but it was not run: there was no emulator here. The workflow runs its
  tests under qemu.
- The CI workflow's steps pass locally as written.

### Bioconda recipe

- **Builder:** bioconda-utils 5.0.0 (conda-build with bioconda's pinnings, conda 26.7.2) in a
  Docker image built from `condaforge/miniforge3`. The command was `bioconda-utils build recipes
  config.yml --packages tussl-chemistry`, run in a scratch bioconda-recipes tree that used
  bioconda's `config.yml`.
- **Source:** `git archive --prefix=tussl-chemistry-0.1.0/ HEAD`, which is the same layout as
  GitHub's tag tarball, served as a `file://` URL. The tested recipe differed from
  `bioconda/meta.yaml` only in `url` and `sha256`.
- **Toolchain picked by the pinnings:** `rust_linux-64` 1.98.1 (Cargo.toml's rust-version is
  1.88, and cargo refuses an older rustc), `gcc_linux-64` 16.2.0, `sysroot_linux-64` 2.17, and
  cargo-bundle-licenses 4.2.0. conda's cargo ignores `rust-toolchain.toml`, so this Rust, not
  1.96.1, built the package.
- **Result:** `tussl-chemistry-0.1.0-hfa8f182_0.conda` (linux-64, 0.5 MB). The package holds
  exactly `bin/tchem`, `bin/tchem-fetch-whitelists`,
  `share/tussl-chemistry/scripts/fetch_whitelists.sh` and
  `share/tussl-chemistry/data/lists.tsv`. `info/licenses/` holds `LICENSE` and `THIRDPARTY.yml`
  (flate2 and zlib-rs, with their licence texts). The run requirements are libgcc, bash, curl,
  gzip, coreutils, grep and gawk.
- **Recipe tests:** these passed in conda-build's test environment: `tchem --version` (and that
  it prints `tchem 0.1.0`), `tchem --help`, the fetcher's usage message, the registry's presence,
  and the check that no barcode list is shipped. The usage test is written so that it also
  passes under `pipefail`, because the fetcher exits 1 after it prints its usage.
- **Fix found by the build:** with Rust alone, conda-build warned that `libgcc_s.so.1` was not in
  any package. Rust's standard library links it. `compiler('c')` now adds `libgcc` to the run
  requirements through its run_exports.
- **Lint:** `bioconda-utils lint` on the recipe as committed (the sha256 placeholder replaced by
  a well-formed dummy) reports all checks OK. The only finding was `long_summary`, which is
  fixed: the summary is now short and the detail is in `description`. The lint does not fetch
  the URL.
- **Stand-in images:** the package was installed into a fresh environment (`mamba create` from
  the build's output channel, with conda-forge and bioconda for the dependencies). That
  environment went into two images. One is on `condaforge/miniforge3`. The other copies the
  environment onto `quay.io/bioconda/base-glibc-busybox-bash`, the base BioContainers uses. That
  base has no distribution userland and no system `libgcc_s`, so only the declared run
  requirements are present. The busybox image is the one tagged
  `quay.io/biocontainers/tussl-chemistry:0.1.0--FILL_IN` and used by nf-test. The recipe's test
  commands pass in it. The fetcher's tools (bash, curl, gunzip, sha256sum, awk, grep, od,
  mktemp, sort) all resolve to the package's environment. Bioconda's own mulled test
  (`--mulled-test`) was not run.
- **Byte-for-byte check:** each image ran with the data mounted read-only at its host path and no
  network. tchem's stdout was compared with the repository's release binary
  (`target/release/tchem`: release profile, toolchain 1.96.1, glibc, built at 56a3eac; not the
  musl release archive). Both got the same arguments. The outputs from both images are
  identical to the release binary's, in JSON and in text, for three samples: `pbmc_1k_v3_S1_L001` R1/R2 with `--whitelists
  $TUSSL_CHEMISTRY_DATA/whitelists` (10x 3' v3/v3.1, HIGH), and the two nf-core
  test-datasets heads below (10x 3' v3/v3.1, HIGH; 10x 5' v1/v2, HIGH).
- **Not built here:** linux-aarch64 and osx-arm64, which the recipe lists as additional
  platforms. Bioconda's CI builds them.

### nf-core modules

- **Template:** nf-core/tools 4.1.0, `nf-core modules create`, in a clone of nf-core/modules at
  3388c6a (2026-10-08). The modules follow the current template: versions go out as an `eval`
  tuple on the `versions` topic (`versions_tchem`), not as `versions.yml`. Containers use the
  `workflow.containerEngine in ['singularity', 'apptainer']` form.
- **`TCHEM_FETCHLISTS`:** takes an optional Cell Ranger install directory (`[]` for none). Only
  the Flex list needs it. The module emits the `whitelists` directory, with every list checked
  and packed.
- **`TCHEM_DETECT`:** takes `tuple val(meta), path(reads)` (one lane: R1 and R2 in any order,
  with index reads set aside, or one long-read FASTQ or BAM) and `path whitelists`. It emits
  `*.chemistry.json` (the report, read by `schema_version`). It never passes `--emit`, which
  exits 3 when no command fits and would fail the task. The report exits 0 whatever the verdict.
- **Tests:** the nf-test tests use existing files on nf-core/test-datasets' `modules` branch.
  They are subsampled 10x heads of 10,000 read pairs:
  `10k_pbmc_cmo/fastqs/gex_1/subsampled_SC3_v3_NextGem_DI_CellPlex_Human_PBMC_10K_1_gex_S2_L001_R{1,2}_001.fastq.gz`
  (3' v3) and
  `5k_cmvpos_tcells/fastqs/gex_1/subsampled_5k_human_antiCMV_T_TBNK_connect_GEX_1_S1_L001_R{1,2}_001.fastq.gz`
  (5' v2). The detect tests assert the verdict, the confidence, the layout and the list, and they
  snapshot the report. The report's paths are relative to the task directory, so its md5 is
  stable.
- **Runs:** nf-test 0.9.5 with Nextflow 26.04.6 (`--profile docker`, the clone's
  `tests/config/nf-test.config`). Both ran in a container built from `condaforge/miniforge3`,
  with the host's Docker socket. The process container was the busybox stand-in image above,
  under the modules' placeholder name. Earlier runs with the miniforge image gave the same
  results. All 5 tests pass: fetch, fetch stub,
  3' v3, 5' v2, and the detect stub. The run with real data fetched every list over the network,
  each checked against its sha256. A second run with the snapshots committed passed and left
  them unchanged. Two extra checks agree with the snapshot: the release binary, run in the task
  directory on the same relative inputs, prints a report with the same md5 (`1402f633…`, 3' v3).
  The verdicts are the ones the release binary gives on the same heads.
- **Lint:** `nf-core modules lint` passes 59 checks for detect and 54 for fetchlists. Each module
  fails one: `bioconda_version`. That check looks `bioconda::tussl-chemistry=0.1.0` up on
  anaconda.org, where it does not exist yet. The `FILL_IN` container tags were not flagged, but
  they must still be filled.
- **Pre-commit checks** (from the clone's `.pre-commit-config.yaml`) all pass on the two
  modules: prettier 3.9.7 with the clone's `.prettierrc.yml`, check-jsonschema 0.38.0 against
  `modules/yaml-schema.json` and `modules/environment-schema.json`, trailing whitespace and
  final newline, and `nextflow lint`. These were run directly, not through `pre-commit`.

## Steps to publish

The maintainer does these steps in this order, and ticks each box when the step is done.

- [ ] **Make the repository public** (Tussl-Co/TusslChemistry). Bioconda, BioContainers and
      nf-core all fetch from it.
- [ ] **Tag `v0.1.0`** and publish the draft release that `.github/workflows/release.yml`
      creates. The workflow checks that the coverage is fresh and builds static Linux and macOS
      archives with their sha256. Read the draft before you publish it.
- [ ] **Re-check at the tagged commit.** The report embeds the validation counts, which the
      bench/coverage work changes. The byte-check references and the detect snapshots' md5s
      therefore hold for 56a3eac only. At the tag, rebuild the package, repeat the byte check,
      and run nf-test again. Regenerate the snapshots if the report changed.
- [ ] **Fill the sha256** in `bioconda/meta.yaml` from
      `https://github.com/Tussl-Co/TusslChemistry/archive/refs/tags/v0.1.0.tar.gz`. Compute
      it with `curl -sL <url> | sha256sum` after the tag exists. It may differ from the sha256
      of the local archive. Never reuse the local value.
- [ ] **Bioconda PR:** copy `bioconda/` to `recipes/tussl-chemistry/` in bioconda-recipes and
      open the pull request. Its CI also builds linux-aarch64 and osx-arm64, which were not built
      here.
- [ ] **BioContainers tags:** after the bioconda merge, BioContainers builds
      `quay.io/biocontainers/tussl-chemistry:0.1.0--<build>` and the Galaxy Singularity image.
      Put both tags in place of `0.1.0--FILL_IN` in the two `main.nf` files. nf-core also
      accepts Seqera containers built from `environment.yml`.
- [ ] **nf-core/test-datasets PR:** probably not needed. The tests use 10x heads that are
      already on the `modules` branch. Whether to rely on them, and the check of their licence
      (10x's public datasets are generally CC BY 4.0), is the maintainer's decision. A new head would
      need its own PR, made before the modules PR.
- [ ] **nf-core/modules PR:** copy `nf-core/modules/tchem/` to `modules/nf-core/tchem/`, fill
      the containers, and run `nf-core modules lint tchem/detect tchem/fetchlists` and
      `nf-test test --profile docker --tag tchem`. Then open the pull request. The one lint
      failure left here, `bioconda_version`, clears once the package is on bioconda. If the
      lists' sources or the tool version change, regenerate the snapshots (they are
      committed). Tell the reviewers that the detect test's `setup` runs `TCHEM_FETCHLISTS`. For the two
      tests with real data, that downloads every list from its third-party source (GitHub raw,
      Box, an S3 bucket). The stub test stubs the fetch too. One unreachable source fails the
      test.

## How to repeat the local checks

The scripts used here are not committed: they are local scaffolding under `target/pkg/`. The
steps are:

1. Build an image `FROM condaforge/miniforge3` with `mamba create -n bu bioconda-utils=5.0.0`.
2. Make a scratch tree with bioconda's `config.yml`, an empty `build-fail-blacklist`, and
   `recipes/tussl-chemistry/` holding this recipe with `url: file:///src/v0.1.0.tar.gz` and the
   archive's sha256. Commit it as a git repository, because bioconda-utils needs one.
3. Run `bioconda-utils build recipes config.yml --packages tussl-chemistry` and
   `bioconda-utils lint recipes config.yml --packages tussl-chemistry`.
4. Install the package from the output channel into an environment, and copy that environment
   onto `quay.io/bioconda/base-glibc-busybox-bash` to stand in for BioContainers. Tag it `quay.io/biocontainers/tussl-chemistry:0.1.0--FILL_IN`, so the
   modules' container string resolves locally.
5. Copy the modules into an nf-core/modules clone, then run nf-core/tools and nf-test in a
   container. The container needs the host Docker socket and the work directory mounted at its
   host path.
