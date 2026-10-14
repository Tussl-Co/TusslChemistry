#!/usr/bin/env bash
# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

# The gate every change passes before it is committed:
#   1. clippy with the crate's strict lint table: zero findings,
#   2. unit tests, including the allocation proofs,
#   3. the parity oracle against the release binary: every report byte-identical,
#   4. the same oracle against the debug binary: overflow checks, debug assertions and the
#      order-independence cross-check all armed, and the reports must still be identical,
#   5. the JSON oracle against the release binary, and every recorded JSON report against the
#      published schema (schema/report.v1.json; needs the jsonschema package),
#   6. data/validation.tsv and the README's recall table against the recorded coverage run and its
#      adjudications, and whether that run is still the code that detects (a warning here;
#      scripts/release_check.sh fails),
#   7. scripts/validate.sh: known samples and damaged copies, each verdict checked,
#   8. the chaos suite (bench/chaos/scripts/chaos.py: damaged FASTQs) and the chaos probes
#      (probes.py: hostile names, paths and bases), against the release binary: no abort, no
#      hang, no HIGH call wrong, and each probe's stated result. Their tables go to a scratch
#      directory; the committed ones are re-recorded by hand when a change is meant,
#   9. the documentation's Simplified Technical English rules (scripts/ste_check.py, terms in
#      docs/glossary.md).
# Usage: TUSSL_CHEMISTRY_DATA=<benchmark data root> scripts/gate.sh [oracle_dir]   (see bench/parity/scripts/oracle.py)
set -euo pipefail
cd "$(dirname "$0")/.."
T=${CARGO_TARGET_DIR:-target}
R=${TUSSL_CHEMISTRY_DATA:?set TUSSL_CHEMISTRY_DATA to the benchmark data root}
WL=${TUSSL_CHEMISTRY_WHITELISTS:-$R/whitelists}
ORACLE=${1:-$R/parity_oracle}
ORACLE_JSON=${ORACLE_JSON:-$R/parity_oracle_json}
echo "== clippy"
n=$(cargo clippy --release --all-targets --message-format=short 2>&1 | grep -cE '^src/' || true)
[ "$n" = 0 ] || { cargo clippy --release --all-targets --message-format=short 2>&1 | grep -E '^src/'; exit 1; }
echo "== tests"; cargo test --release 2>&1 | grep -E 'test result'
echo "== build"; cargo build --release -q && cargo build -q
echo "== oracle (release)"; python3 -I bench/parity/scripts/oracle.py compare "$T/release/tchem" "$WL" "$ORACLE"
echo "== oracle (debug)"; python3 -I bench/parity/scripts/oracle.py compare "$T/debug/tchem" "$WL" "$ORACLE"
echo "== oracle (json, release)"
ORACLE_ARGS="--format json" python3 -I bench/parity/scripts/oracle.py compare "$T/release/tchem" "$WL" "$ORACLE_JSON"
echo "== schema"; python3 -I bench/parity/scripts/schema_check.py schema/report.v1.json "$ORACLE_JSON"
echo "== validation"
python3 -I bench/coverage/scripts/tally.py bench/coverage/results/calls.tsv \
  bench/coverage/results/adjudication.tsv data/validation.tsv --check | tail -1
python3 -I bench/coverage/scripts/recall.py bench/coverage/results/calls.tsv \
  bench/coverage/results/adjudication.tsv bench/coverage/results/recall.tsv README.md --check
scripts/coverage_fresh.sh || echo "   (warning only: scripts/release_check.sh refuses a release on it)"
SCRATCH=$(mktemp -d)
trap 'rm -rf "$SCRATCH"' EXIT
# step <name> <command...>: the last line on success; on failure every line that did not pass.
step() {
  local name=$1; shift
  echo "== $name"
  "$@" > "$SCRATCH/$name.log" 2>&1 || { grep -v '^PASS' "$SCRATCH/$name.log"; exit 1; }
  tail -1 "$SCRATCH/$name.log"
}
step docs python3 -I scripts/ste_check.py
step validate zsh scripts/validate.sh
step probes python3 -I bench/chaos/scripts/probes.py "$T/release/tchem" "$WL" "$SCRATCH/p" \
  "$SCRATCH/probes.tsv"
step chaos python3 -I bench/chaos/scripts/chaos.py "$T/release/tchem" "$WL" "$SCRATCH/c" \
  "$SCRATCH/chaos.tsv"
echo "== gate passed"
