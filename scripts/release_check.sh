#!/usr/bin/env bash
# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

# What a release needs beyond the gate. Each check fails the release; none is a warning.
#   1. the tree is committed, so the release is a commit someone can name,
#   2. the changelog has a section for the crate's version,
#   3. both published schemas are valid JSON Schemas,
#   4. the coverage run behind data/validation.tsv is the code being released, so the
#      validation every match prints describes the binary that ships.
# Usage: scripts/release_check.sh   (run scripts/gate.sh first)
set -euo pipefail
cd "$(dirname "$0")/.."
git diff --quiet HEAD -- || { echo "release: uncommitted changes"; exit 1; }
version=$(grep -m1 '^version' Cargo.toml | cut -d'"' -f2)
grep -q "^## ${version}" CHANGELOG.md || { echo "release: CHANGELOG.md has no ## ${version}"; exit 1; }
python3 -I -c "import json, jsonschema
for f in ['schema/report.v1.json', 'schema/emit.v1.json']:
    jsonschema.Draft202012Validator.check_schema(json.load(open(f)))"
scripts/coverage_fresh.sh
echo "== release checks passed (tchem ${version})"
