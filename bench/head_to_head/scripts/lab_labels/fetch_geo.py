# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Fetch GEO brief text for every GSM linked from the input runs' ENA records (sample_alias,
experiment_alias, or experiment_title). At most 3 requests/s; one retry. Body must start with ^SAMPLE = <GSM>.
Usage: python -I fetch_geo.py <auditdir>"""
import csv, glob, json, re, sys, time, urllib.request

D = sys.argv[1]
inp = {r["srx"]: r for r in csv.DictReader(open(f"{D}/input_runs.tsv"), delimiter="\t")}
gsms = set()
for p in sorted(glob.glob(f"{D}/raw/ena/chunk_*.tsv")):
    for r in csv.DictReader(open(p), delimiter="\t"):
        x = r["experiment_accession"]
        if x in inp and r["run_accession"] == inp[x]["srr"]:
            gsms |= set(re.findall(r"GSM\d+", " ".join([r["sample_alias"], r["experiment_alias"], r["experiment_title"]])))
status = {}
last = 0.0
def get(url):
    global last
    wait = 0.34 - (time.time() - last)
    if wait > 0: time.sleep(wait)
    last = time.time()
    with urllib.request.urlopen(url, timeout=60) as resp:
        return resp.read().decode("utf-8", "replace")
for g in sorted(gsms):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=brief"
    err = ""
    for attempt in range(2):
        try:
            body = get(url)
            if re.match(rf"\^SAMPLE = {g}\s", body):
                open(f"{D}/raw/geo/{g}.txt", "w").write(body); err = ""; break
            err = "body does not start with ^SAMPLE = " + g + ": " + body[:80].replace("\n", " ")
        except Exception as e:
            err = repr(e)
        time.sleep(1)
    status[g] = {"url": url, "error": err}
json.dump(status, open(f"{D}/raw/geo/index.json", "w"), indent=1)
print(len(status), "GSM;", sum(1 for s in status.values() if s["error"]), "failed")
