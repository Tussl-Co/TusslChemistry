#!/usr/bin/env bash
# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

# Is the coverage run behind data/validation.tsv still the code that detects? The run records its
# commit (bench/coverage/results/calls.tsv, first line). If anything that can change a call has
# changed since (the library outside its output writers, the reference data, the dependencies),
# the validation printed with every match describes an older binary: exit 1 and say so.
# Given a commit instead, it asks the same of that commit: bench/coverage/scripts/run.py --reuse
# carries a call over only from a run whose detection is the same as now.
# Usage: scripts/coverage_fresh.sh [commit]   (exit 0 fresh, 1 stale)
set -euo pipefail
cd "$(dirname "$0")/.."
calls=bench/coverage/results/calls.tsv
if [ $# -ge 1 ]; then
  commit=$1
  git rev-parse --verify --quiet "$commit^{commit}" > /dev/null \
    || { echo "coverage: $commit is not a commit"; exit 1; }
else
  commit=$(sed -n '1s/.*commit \([0-9a-f]\{7,\}\).*/\1/p' "$calls")
  [ -n "$commit" ] || { echo "coverage: no commit in the header of $calls"; exit 1; }
fi
# Paths that cannot change a call: the writers, --emit, the accession fetcher (coverage reads
# files on disk), the table itself, and the help text (cli/usage.rs holds only that string). The
# rest of the command line and main.rs stay in: they decide which file by position is R1, R2 or
# an index read, and the coverage run gives files by position.
if git diff --quiet "$commit" -- src data Cargo.toml Cargo.lock \
    ':!src/output' ':!src/emit' ':!src/remote.rs' ':!src/cli/usage.rs' \
    ':!src/layouts/validation.rs' \
    ':!data/validation.tsv' ':!data/validation_family.tsv'; then
  echo "coverage: fresh (run at $commit; nothing that detects has changed since)"
else
  echo "coverage: STALE, run at $commit; detection has changed since. Re-run bench/coverage:"
  git diff --stat "$commit" -- src data Cargo.toml Cargo.lock ':!src/output' ':!src/emit' ':!src/remote.rs' \
    ':!src/cli/usage.rs' \
    ':!src/layouts/validation.rs' ':!data/validation.tsv' ':!data/validation_family.tsv' | tail -1
  exit 1
fi
