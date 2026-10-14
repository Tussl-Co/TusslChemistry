# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Runs tchem on every fetched head (2,000 reads each, JSON) and records the call beside the
catalogue's stated kit. Heads are truncated mid-file, so the read count stays well inside each.

Usage: python run.py <catalogue.tsv> <heads_dir> <tchem> <whitelist_dir> <out.tsv>
"""

import csv
import json
import os
import subprocess
import sys

sys.path.insert(0, os.path.dirname(__file__))
from fetch_heads import head_name, url_of  # noqa: E402

READS = 2000
FIELDS = ["run", "platform", "data_form", "kit", "kit_evidence", "status", "read_kind", "verdict",
          "confidence", "layout", "long_read_kit", "adapter_share", "list_shares", "poly_t_share",
          "tso_share", "umi_length", "read_length_mode", "hint"]


def run_one(row, heads, tool, whitelists):
    url = url_of(row)
    out = {k: row.get(k, "") for k in ("run", "platform", "data_form", "kit", "kit_evidence")}
    path = os.path.join(heads, head_name(row, url)) if url else None
    if path is None or not os.path.exists(path):
        return {**out, "status": "no head"}
    p = subprocess.run([tool, "--whitelists", whitelists, "--r1", path, "--reads", str(READS),
                        "--format", "json"], capture_output=True, text=True, timeout=600)
    if p.returncode != 0:
        return {**out, "status": "refused: " + p.stderr.strip().splitlines()[0][:150]}
    d = json.loads(p.stdout)
    e = d.get("long_read_evidence") or {}
    shares = ";".join(f"{x['list']}={x['share']:.3f}" for x in e.get("barcode_list_shares", []))
    return {**out, "status": "ok", "read_kind": d["read_kind"], "verdict": d["verdict"],
            "confidence": d["confidence"] or "", "layout": (d["layout"] or {}).get("name", ""),
            "long_read_kit": (d["layout"] or {}).get("long_read_kit") or "",
            "adapter_share": e.get("adapter_share", ""), "list_shares": shares,
            "poly_t_share": e.get("poly_t_share", ""), "tso_share": e.get("tso_share", ""),
            "umi_length": e.get("umi_length", ""), "read_length_mode": e.get("read_length_mode", ""),
            "hint": d.get("hint") or ""}


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
