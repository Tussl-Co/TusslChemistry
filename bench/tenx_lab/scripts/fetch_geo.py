# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Copy of bench/head_to_head/scripts/lab_labels/fetch_geo.py; changes: GSMs from every run (not one srr) of
the input runs (work/input_runs.tsv, after exclusions), cached GSM texts reused (resume). URL, ^SAMPLE check, 0.34 s throttle
(<3 requests/s), one retry: unchanged. Usage: python3 -I -B fetch_geo.py <W>"""
import csv, glob, json, os, re, sys, time, urllib.request

D = sys.argv[1]
csv.field_size_limit(10**8)
rem = {r["srr"] for r in csv.DictReader(open(f"{D}/work/input_runs.tsv", encoding="utf-8"), delimiter="\t")}
gsms = set()
for p in sorted(glob.glob(f"{D}/raw/ena/chunk_*.tsv")):
    for r in csv.DictReader(open(p, encoding="utf-8"), delimiter="\t"):
        if r["run_accession"] in rem:
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
for n, g in enumerate(sorted(gsms)):
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=brief"
    if os.path.exists(f"{D}/raw/geo/{g}.txt"):
        status[g] = {"url": url, "error": ""}; continue
    err = ""
    for attempt in range(2):
        try:
            body = get(url)
            if re.match(rf"\^SAMPLE = {g}\s", body):
                open(f"{D}/raw/geo/{g}.txt", "w", encoding="utf-8").write(body); err = ""; break
            err = "body does not start with ^SAMPLE = " + g + ": " + body[:80].replace("\n", " ")
        except Exception as e:
            err = repr(e)
        time.sleep(1)
    status[g] = {"url": url, "error": err}
    if n % 200 == 0: print(n, flush=True)
json.dump(status, open(f"{D}/raw/geo/index.json", "w"), indent=1)
print(len(status), "GSM;", sum(1 for s in status.values() if s["error"]), "failed")
