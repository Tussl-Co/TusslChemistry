# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Runs one binary over every held-out case (sources.py) and records each call, so the validation
counts describe the binary that ships, not the binary each benchmark first ran. The output's
header names the commit and the binary's checksum; tally.py reads only this table.

Re-running measures; it does not tune. A run whose failure is used to change the code moves to its
catalogue's design split before the next coverage run, never after.

With --reuse, a case's call is carried over from an earlier coverage table instead of being run
again, when nothing that detects has changed between that table's commit and now
(scripts/coverage_fresh.sh <commit>). A call depends only on the reads and on the code that
detects; the case's label (truth, referee) is always taken from the sources as they are now. Each
row names the commit it was called at (`called_at`). A case not in the earlier table, or whose row
there is not a clean result (ok or refused), is run now. Then VERIFY_SAMPLE carried cases whose
reads are still on disk are run again, and any call that differs stops the run. This lets the
read heads of scored runs be deleted to free disk; they are fetched again
(bench/heldout/scripts/fetch_heads.py) when detection changes.

Usage: python run.py <tchem> <whitelist_dir> <data_dir> <out.tsv> [--reuse <calls.tsv>]
(TCHEM_COVERAGE_WORKERS sets how many cases run at once; default 8.)
"""

import csv
import hashlib
import json
import os
import random
import re
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sources  # noqa: E402

FIELDS = ["source", "run", "study", "truth", "resolution", "referee", "read_kind", "status",
          "verdict", "confidence", "layout", "called_at"]
CALL_FIELDS = ("status", "verdict", "confidence", "layout")  # what the binary decides
CARRIED_STATUSES = ("ok", "refused")  # clean results; a missing head or a timeout is run again
VERIFY_SAMPLE = 200  # carried cases run again to check the carried call is what the binary says
VERIFY_SEED = 20261010
TRUTH_SEPARATOR = " | "  # layout names contain commas, slashes and quotes, never " | "
TIMEOUT_SECONDS = 600  # a short-read sample takes about 2 s; a stuck run is recorded, not waited on
# Cases called at once; each is one tchem process. The calls do not depend on it.
WORKERS = int(os.environ.get("TCHEM_COVERAGE_WORKERS", "8"))


def call_row(case):
    """The case's own columns, as the sources give them now."""
    return {"source": case.source, "run": case.run, "study": case.study,
            "truth": TRUTH_SEPARATOR.join(case.truth), "resolution": case.resolution,
            "referee": case.referee, "read_kind": case.read_kind}


def call(tool, whitelists, case):
    row = call_row(case)
    if not all(os.path.exists(f) for f in case.files):
        return {**row, "status": "no reads cached"}
    if len(case.files) == 2:
        inputs = ["--r1", case.files[0], "--r2", case.files[1]]
    else:
        inputs = list(case.files)
    if case.reads:
        inputs += ["--reads", str(case.reads)]
    try:
        run = subprocess.run([tool, "--whitelists", whitelists, "--format", "json", *inputs],
                             capture_output=True, text=True, timeout=TIMEOUT_SECONDS)
    except subprocess.TimeoutExpired:
        return {**row, "status": "timed out"}
    if run.returncode != 0:
        reason = next((l for l in run.stderr.splitlines() if l.startswith("Error")), "")
        # The data root as the recorded tables name it, so a refusal reads the same on any machine.
        data = os.environ.get("TUSSL_CHEMISTRY_DATA")
        if data:
            reason = reason.replace(data.rstrip("/"), "$TUSSL_CHEMISTRY_DATA")
        return {**row, "status": f"refused: {reason[:160]}"}
    report = json.loads(run.stdout)
    layout = report["layout"]["name"] if report["layout"] else ""
    return {**row, "status": "ok", "verdict": report["verdict"],
            "confidence": report["confidence"] or "", "layout": layout}


def provenance(tool):
    """The commit (a plain hash: a tag must not change how it reads), whether the tree has
    uncommitted changes, the binary's checksum and its version."""
    repo = os.path.dirname(os.path.abspath(__file__))
    git = lambda *a: subprocess.run(["git", "-C", repo, *a], capture_output=True, text=True)
    commit = git("rev-parse", "--short", "HEAD").stdout.strip()
    dirty = bool(git("status", "--porcelain", "--untracked-files=no").stdout.strip())
    invariant(re.fullmatch(r"[0-9a-f]{7,40}", commit), f"not a commit hash: {commit!r}")
    with open(tool, "rb") as f:
        digest = hashlib.sha256(f.read()).hexdigest()[:16]
    version = subprocess.run([tool, "--version"], capture_output=True, text=True).stdout.strip()
    return commit, dirty, digest, version


def detection_unchanged_since(commit):
    """Whether nothing that detects has changed between `commit` and the tree (coverage_fresh.sh)."""
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")
    script = os.path.join(root, "scripts", "coverage_fresh.sh")
    fresh = subprocess.run([script, commit], capture_output=True, text=True)
    invariant(fresh.returncode in (0, 1), f"coverage_fresh.sh failed: {fresh.stderr.strip()}")
    return fresh.returncode == 0


def carried_calls(reuse_path):
    """(source, run, read kind) -> the call columns and the commit it was made at, from an earlier
    coverage table whose detection is the same as now; only clean results are carried."""
    with open(reuse_path) as f:
        header = f.readline()
        rows = list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))
    found = re.search(r"commit ([0-9a-f]{7,40})", header)
    invariant(found, f"{reuse_path}: no commit in its header")
    run_commit = found.group(1)
    carried = {}
    for r in rows:
        if r["status"].split(":")[0] not in CARRIED_STATUSES:
            continue
        called_at = r.get("called_at") or run_commit
        carried[(r["source"], r["run"], r["read_kind"])] = (
            {k: r.get(k, "") for k in CALL_FIELDS}, called_at)
    # Each row was called at a commit whose detection matched its table's run; checking every
    # commit named here against the tree keeps that true however many runs a call is carried.
    for called_at in sorted({c for _, c in carried.values()}):
        invariant(detection_unchanged_since(called_at),
                  f"{reuse_path}: calls made at {called_at}, and detection has changed since")
    return carried


def verify_carried(tool, whitelists, cases, results):
    """Run a seeded sample of carried cases whose reads are on disk again; every call must match."""
    carried = [(c, r) for c, r in zip(cases, results)
               if r["called_at"] != "now" and all(os.path.exists(f) for f in c.files)]
    sample = random.Random(VERIFY_SEED).sample(carried, min(VERIFY_SAMPLE, len(carried)))
    with ThreadPoolExecutor(WORKERS) as pool:
        again = list(pool.map(lambda cr: call(tool, whitelists, cr[0]), sample))
    for (case, kept), fresh in zip(sample, again):
        invariant(all(kept.get(k, "") == fresh.get(k, "") for k in CALL_FIELDS),
                  f"{case.run}: carried call {[kept.get(k) for k in CALL_FIELDS]} but the binary "
                  f"now says {[fresh.get(k) for k in CALL_FIELDS]}")
    return len(sample)


def main(tool, whitelists, data, out_path, reuse_path=None):
    commit, dirty, digest, version = provenance(tool)
    invariant(not dirty, "run coverage from a committed tree, so it can be named")
    cases = sources.cases(data)
    invariant(cases, "no held-out cases found")
    carried = carried_calls(reuse_path) if reuse_path else {}

    def one(case):
        kept = carried.get((case.source, case.run, case.read_kind))
        if kept is None:
            return {**call(tool, whitelists, case), "called_at": "now"}
        calls, called_at = kept
        row = call_row(case)
        return {**row, **calls, "called_at": called_at}

    with ThreadPoolExecutor(WORKERS) as pool:
        results = list(pool.map(one, cases))
    invariant(len(results) == len(cases), "one row per case")
    verified = verify_carried(tool, whitelists, cases, results) if carried else 0
    for r in results:
        if r["called_at"] == "now":
            r["called_at"] = commit
    n_carried = sum(r["called_at"] != commit for r in results)
    with open(out_path, "w") as f:
        f.write(f"# Coverage run: {version}, commit {commit}, binary sha256 {digest}.\n")
        if reuse_path:
            f.write(f"# {n_carried} calls carried from earlier runs with the same detection "
                    f"(called_at); {len(results) - n_carried} run now; {verified} carried calls "
                    "run again and matched.\n")
        f.write("# One row per held-out run; truth lists the layouts it may rightly be called.\n")
        f.write("\t".join(FIELDS) + "\n")
        for r in results:
            f.write("\t".join(str(r.get(k, "")) for k in FIELDS) + "\n")
    statuses = {}
    for r in results:
        statuses[r["status"].split(":")[0]] = statuses.get(r["status"].split(":")[0], 0) + 1
    print(f"{len(results)} cases at {commit}: {statuses}")


if __name__ == "__main__":
    if len(sys.argv) == 5:
        main(*sys.argv[1:])
    elif len(sys.argv) == 7 and sys.argv[5] == "--reuse":
        main(*sys.argv[1:5], reuse_path=sys.argv[6])
    else:
        sys.exit(__doc__)
