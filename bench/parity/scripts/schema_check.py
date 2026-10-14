# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Checks every JSON report in a parity oracle against the published schema, so the contract in
schema/ and what the binary writes cannot drift apart. Reports that are refusals (no JSON on
stdout) are counted and skipped. Fails unless short and long reports were both checked.

Usage: python schema_check.py <schema.json> <json_oracle_dir>
Needs the jsonschema package (pip install jsonschema).
"""

import glob
import json
import os
import sys

import jsonschema

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

REPORTS_MAX = 10_000  # far above the oracle's 264 cases; a larger directory is the wrong one


STDERR_MARK = "\n--- stderr errors\n"  # oracle.py's separator between stdout and error lines


def body(path):
    """The report's stdout: oracle.py writes `rc=<code>`, stdout, then the stderr error lines."""
    with open(path, encoding="utf-8") as f:
        first, rest = f.readline(), f.read()
    invariant(first.startswith("rc="), f"{path}: not an oracle case file")
    stdout, mark, _ = rest.rpartition(STDERR_MARK)
    invariant(mark, f"{path}: no stderr section")
    return stdout


def main(schema_path, oracle_dir):
    with open(schema_path, encoding="utf-8") as f:
        schema = json.load(f)
    jsonschema.Draft202012Validator.check_schema(schema)
    validator = jsonschema.Draft202012Validator(schema)
    paths = sorted(glob.glob(os.path.join(oracle_dir, "*.out")))
    invariant(paths, f"no oracle cases in {oracle_dir}")
    invariant(len(paths) <= REPORTS_MAX, f"{len(paths)} cases: is {oracle_dir} the JSON oracle?")
    kinds, refusals, failures = {"short": 0, "long": 0}, 0, []
    for path in paths:
        text = body(path).strip()
        if not text.startswith("{"):
            refusals += 1
            continue
        report = json.loads(text)
        errors = sorted(validator.iter_errors(report), key=lambda e: list(e.path))
        if errors:
            failures.append((os.path.basename(path), errors[0]))
            continue
        kinds[report["read_kind"]] += 1
    for name, error in failures[:20]:
        where = "/".join(str(p) for p in error.path) or "(top)"
        print(f"FAIL  {name}: {where}: {error.message}")
    checked = sum(kinds.values())
    print(f"schema: {checked} reports valid ({kinds['short']} short, {kinds['long']} long), "
          f"{refusals} refusals, {len(failures)} invalid")
    invariant(checked + refusals + len(failures) == len(paths), "every case is counted once")
    if failures or kinds["short"] == 0 or kinds["long"] == 0:
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
