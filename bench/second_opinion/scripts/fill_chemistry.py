# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fills cellranger.tsv's chemistry column from each run's outs/web_summary.html ("Chemistry" row).

Usage: python fill_chemistry.py <cellranger.tsv> <workdir> <out.tsv>
"""

import csv
import os
import re
import sys


def main(path, work, out):
    rs = list(csv.DictReader(open(path), delimiter="\t"))
    for r in rs:
        p = os.path.join(work, r["run"], "outs", "web_summary.html")
        m = re.search(r'"Chemistry"\s*,\s*"([^"]+)"', open(p).read()) if os.path.exists(p) else None
        r["cellranger"] = m.group(1) if m else r["cellranger"]
    with open(out, "w") as f:
        w = csv.DictWriter(f, list(rs[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rs)


if __name__ == "__main__":
    main(*sys.argv[1:4])
