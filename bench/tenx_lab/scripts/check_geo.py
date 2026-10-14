# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Refetch every GSM text written before the rerun (mtime < cutoff) and compare byte for byte. 0.34 s throttle,
one retry. Usage: python3 -I -B check_geo.py <W> <cutoff_epoch>"""
import glob, os, sys, time, urllib.request
D, cut = sys.argv[1], float(sys.argv[2])
files = sorted(p for p in glob.glob(f"{D}/raw/geo/GSM*.txt") if os.path.getmtime(p) < cut)
last = 0.0; diff = []; fail = []
for p in files:
    g = os.path.basename(p)[:-4]
    url = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={g}&targ=self&form=text&view=brief"
    body = None
    for _ in range(2):
        try:
            w = 0.34 - (time.time() - last)
            if w > 0: time.sleep(w)
            last = time.time()
            with urllib.request.urlopen(url, timeout=60) as r:
                body = r.read().decode("utf-8", "replace"); break
        except Exception as e:
            err = repr(e); time.sleep(1)
    if body is None:
        fail.append((g, err)); continue
    with open(f"{D}/raw/geo_check/{g}.txt", "w", encoding="utf-8") as f:
        f.write(body)
    if open(p, encoding="utf-8").read() != body:
        diff.append(g)
print(len(files), "checked;", len(diff), "differ:", diff, ";", len(fail), "failed:", fail)
