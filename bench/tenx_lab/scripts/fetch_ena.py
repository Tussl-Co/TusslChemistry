# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Copy of bench/head_to_head/scripts/lab_labels/fetch_ena.py; changes: input from work/prep.json,
fields + fastq_ftp,first_public, cached chunks are reused (resume). Usage: python3 -I -B fetch_ena.py <W>"""
import json, os, sys, time, urllib.parse, urllib.request

D = sys.argv[1]
FIELDS = "run_accession,experiment_accession,experiment_title,sample_title,library_name,library_construction_protocol,sample_alias,study_title,experiment_alias,study_accession,fastq_ftp,first_public"
srxs = json.load(open(f"{D}/work/prep.json"))["query"]
out = []
for i in range(0, len(srxs), 40):
    chunk = srxs[i:i + 40]
    q = " OR ".join(f"experiment_accession={x}" for x in chunk)
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + urllib.parse.urlencode(
        {"result": "read_run", "query": q, "fields": FIELDS, "format": "tsv", "limit": 0})
    path = f"{D}/raw/ena/chunk_{i//40:03d}.tsv"
    body, err = None, ""
    if os.path.exists(path):
        body = open(path, encoding="utf-8").read()
    else:
        for attempt in range(2):
            try:
                with urllib.request.urlopen(url, timeout=120) as resp:
                    body = resp.read().decode("utf-8")
                break
            except Exception as e:
                err = repr(e); time.sleep(3)
        if body is not None:
            open(path, "w", encoding="utf-8").write(body)
        time.sleep(0.5)
    out.append({"chunk": i // 40, "url": url, "srxs": chunk, "path": path if body is not None else "", "error": err if body is None else ""})
json.dump(out, open(f"{D}/raw/ena/index.json", "w"), indent=1)
print(len(out), "chunks;", sum(1 for o in out if o["error"]), "failed")
