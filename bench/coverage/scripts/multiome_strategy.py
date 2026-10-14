# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""What ENA records as the library strategy of each held-out run labelled 10x Multiome GEX, beside
what tchem said. A Multiome label comes from protocol text, which covers both libraries of a
Multiome experiment: gene expression and ATAC. tchem detects only the gene-expression library, so
a run that is the ATAC library has no right match. This table shows how many of the no-calls ENA
itself records as ATAC.

Reads ENA's metadata only (portal API, read_run: library_source, library_strategy); never reads.

Usage: python multiome_strategy.py <calls.tsv> <out.tsv>
"""

import collections
import os
import sys
import time
import urllib.parse
import urllib.request

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from invariant import invariant  # noqa: E402
import tally  # noqa: E402

LABEL = "10x Multiome GEX"
PORTAL = "https://www.ebi.ac.uk/ena/portal/api/search"
BATCH_RUNS = 40
ATTEMPTS = 4
RESPONSE_BYTES_MAX = 1 << 20


def labelled(calls):
    rows = [c for c in calls if c["truth"] == LABEL and tally.independent(c["referee"])]
    invariant(len(rows) > 0, f"held-out runs labelled {LABEL}")
    return rows


def strategies(runs):
    """run -> (library_source, library_strategy), for the runs ENA answers for."""
    found = {}
    for start in range(0, len(runs), BATCH_RUNS):
        batch = runs[start:start + BATCH_RUNS]
        query = " OR ".join(f"run_accession={run}" for run in batch)
        body = urllib.parse.urlencode({"result": "read_run", "query": query, "format": "tsv",
                                       "fields": "run_accession,library_source,library_strategy",
                                       "limit": 0}).encode()
        text = None
        for attempt in range(ATTEMPTS):
            try:
                with urllib.request.urlopen(PORTAL, data=body, timeout=120) as response:
                    text = response.read(RESPONSE_BYTES_MAX + 1).decode()
                break
            except OSError:
                time.sleep(5 * (attempt + 1))
        invariant(text is not None and len(text) <= RESPONSE_BYTES_MAX, f"ENA answered batch {start}")
        for line in text.splitlines()[1:]:
            run, source, strategy = line.split("\t")
            invariant(run in batch, f"{run}: ENA answers only for runs asked")
            found[run] = (source, strategy)
    invariant(len(found) <= len(runs), "no more answers than runs")
    return found


def outcome(call):
    if not call["status"].startswith("ok"):
        return "no_call"
    return "match" if call["verdict"] == "match" else "no_call" if call["verdict"] == "no_match" \
        else call["verdict"]


def main(calls_path, out_path):
    rows = labelled(tally.rows(calls_path))
    found = strategies([r["run"] for r in rows])
    out = [f"# Held-out runs labelled {LABEL}: ENA's library strategy beside tchem's outcome "
           "(bench/coverage/scripts/multiome_strategy.py).",
           f"# From: {tally.comment(calls_path).lstrip('# ')}",
           f"# ENA queried {time.strftime('%Y-%m-%d')}; '-' where ENA gave no record.",
           "run\tstudy\tlibrary_source\tlibrary_strategy\toutcome"]
    summary = collections.Counter()
    for r in sorted(rows, key=lambda r: r["run"]):
        source, strategy = found.get(r["run"], ("-", "-"))
        out.append("\t".join([r["run"], r["study"], source, strategy, outcome(r)]))
        summary[(strategy, outcome(r))] += 1
    invariant(len(out) == len(rows) + 4, "one row per labelled run")
    with open(out_path, "w") as f:
        f.write("\n".join(out) + "\n")
    for (strategy, result), n in sorted(summary.items()):
        print(f"{strategy:16} {result:10} {n}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
