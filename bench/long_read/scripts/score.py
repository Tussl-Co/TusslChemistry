# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Scores the long-read benchmark: each call against the catalogue's stated kit. Runs whose kit no
source states ("uncertain") are reported but not scored.

Usage: python score.py <results.tsv>
"""

import collections
import csv
import sys

TRUTH = {  # the catalogue's kit wording -> tchem's layout name
    "3prime v2": "10x 3' v2",
    "3prime v3": "10x 3' v3/v3.1",
    "3prime v3.1": "10x 3' v3/v3.1",
    "3prime v3 (v3 or v3.1)": "10x 3' v3/v3.1",
    "3prime v4 (GEM-X)": "10x 3' v4 (GEM-X)",
    "5prime v2": "10x 5' v1/v2",
    "5prime v3 (GEM-X)": "10x 5' v3 (GEM-X)",
}


def main(path):
    rows = list(csv.DictReader(open(path), delimiter="\t"))
    tally = collections.Counter()
    for r in rows:
        truth = TRUTH.get(r["kit"])
        if r["status"] != "ok":
            outcome = "not run: " + r["status"][:60]
        elif truth is None:
            outcome = f"unscored (kit {r['kit']!r}): called {r['layout'] or r['verdict']}"
        elif r["verdict"] == "match" and r["layout"] == truth:
            outcome = f"right ({r['confidence']})"
        elif r["verdict"] == "match":
            outcome = f"WRONG ({r['confidence']}): called {r['layout']}"
        else:
            outcome = f"no call ({r['verdict']})"
        tally[outcome.split(":")[0] if outcome.startswith(("not run", "unscored")) else outcome] += 1
        print(f"{r['run'][:46]:46} {r['platform'][:6]:6} {(truth or r['kit'])[:22]:22} {outcome}")
    print("\n" + "\n".join(f"{n:4}  {k}" for k, n in tally.most_common()))


if __name__ == "__main__":
    main(sys.argv[1])
