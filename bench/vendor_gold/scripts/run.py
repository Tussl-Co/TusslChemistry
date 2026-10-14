# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Runs tchem on every fetched head pair (JSON) and records the call beside the vendor's stated
chemistry.

Usage: python run.py <catalogue.tsv> <heads_dir> <tchem> <whitelist_dir> <out.tsv>
"""

import csv
import json
import os
import re
import subprocess
import sys

FIELDS = ["dataset", "vendor", "chemistry", "species", "sample_type", "status", "verdict",
          "confidence", "layout", "candidates", "swapped", "cdna_strand", "library", "notes", "hint"]


def run_one(row, heads, tool, whitelists):
    name = re.sub(r"[^A-Za-z0-9._-]", "_", row["dataset"])
    r1, r2 = (os.path.join(heads, f"{name}_{r}.fq.gz") for r in ("R1", "R2"))
    out = {k: row.get(k, "") for k in ("dataset", "vendor", "chemistry", "species", "sample_type")}
    if not (os.path.exists(r1) and os.path.exists(r2)):
        return {**out, "status": "no head"}
    p = subprocess.run([tool, "--whitelists", whitelists, "--r1", r1, "--r2", r2, "--format", "json"],
                       capture_output=True, text=True, timeout=900)
    if p.returncode != 0:
        return {**out, "status": "refused: " + p.stderr.strip().splitlines()[0][:150]}
    d = json.loads(p.stdout)
    return {**out, "status": "ok", "verdict": d["verdict"], "confidence": d["confidence"] or "",
            "layout": (d["layout"] or {}).get("name", ""), "candidates": ";".join(d["candidates"]),
            "swapped": d["swapped"], "cdna_strand": d["cdna_strand"] or "",
            "library": (d["library"] or {}).get("kind", ""), "notes": ";".join(d["notes"]),
            "hint": d["hint"] or ""}


def main(catalogue, heads, tool, whitelists, out_path):
    rows = list(csv.DictReader(open(catalogue), delimiter="\t"))
    with open(out_path, "w") as out:
        w = csv.DictWriter(out, FIELDS, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(run_one(row, heads, tool, whitelists))
            out.flush()


if __name__ == "__main__":
    main(*sys.argv[1:6])
