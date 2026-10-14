# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fetch ENA read_run metadata for every input SRX, in chunks of 40, via the portal API search.
Usage: python -I fetch_ena.py <auditdir>"""
import csv, json, os, sys, time, urllib.parse, urllib.request

D = sys.argv[1]
FIELDS = "run_accession,experiment_accession,experiment_title,sample_title,library_name,library_construction_protocol,sample_alias,study_title,experiment_alias,study_accession"
rows = list(csv.DictReader(open(f"{D}/input_runs.tsv"), delimiter="\t"))
srxs = [r["srx"] for r in rows]
out = []
for i in range(0, len(srxs), 40):
    chunk = srxs[i:i + 40]
    q = " OR ".join(f"experiment_accession={x}" for x in chunk)
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(
        {"result": "read_run", "query": q, "fields": FIELDS, "format": "tsv", "limit": 0})
    body, err = None, ""
    for attempt in range(2):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                body = resp.read().decode("utf-8")
            break
        except Exception as e:
            err = repr(e); time.sleep(3)
    path = f"{D}/raw/ena/chunk_{i//40:03d}.tsv"
    if body is not None:
        open(path, "w").write(body)
    out.append({"chunk": i // 40, "url": url, "srxs": chunk, "path": path if body is not None else "", "error": err if body is None else ""})
    time.sleep(0.5)
json.dump(out, open(f"{D}/raw/ena/index.json", "w"), indent=1)
print(len(out), "chunks;", sum(1 for o in out if o["error"]), "failed")
