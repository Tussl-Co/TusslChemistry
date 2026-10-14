# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Scorecard for the Arc benchmark: merged ENA + SRA-recovered results, with hand adjudications.

Usage: python score.py <results.tsv> <results_sra.tsv> <adjudication.tsv>
"""
import csv
import sys
from collections import Counter

main = {r["srx"]: r for r in csv.DictReader(open(sys.argv[1]), delimiter="\t")}
sra = {r["srx"]: r for r in csv.DictReader(open(sys.argv[2]), delimiter="\t")}
adj = {r["srx"]: r for r in csv.DictReader((l for l in open(sys.argv[3]) if not l.startswith("#")), delimiter="\t")}
rows = [sra.get(s, r) for s, r in main.items()]
usable = [r for r in rows if r.get("verdict")]
print(f"samples {len(rows)}; usable {len(usable)} (ENA as exported {sum(r['srx'] not in sra for r in usable)}, "
      f"recovered from SRA {sum(r['srx'] in sra for r in usable)}); unusable {len(rows) - len(usable)}")
print(f"ENA export dropped the barcode read: {sum(1 for r in main.values() if r['condition'] == 'single read file')} of {len(main)}")

gex = [r for r in usable if r["tech_10x"] in ("3_prime_gex", "5_prime_gex") and r["arc_whitelist"] != "unassigned"]
wl = sum(r["our_whitelist"] == r["arc_whitelist"] for r in gex)
amb = [r for r in gex if r["our_end"] == "ambiguous"]
called = [r for r in gex if r["our_end"] in ("3'", "5'")]
end_agree = sum(r["end_agree"] == "True" for r in called)
print(f"\nIn-scope GEX ({len(gex)}): whitelist matches Arc's read-based choice {wl}/{len(gex)}")
print(f"  3'/5': called {len(called)}, agree with Arc's metadata label {end_agree}; ambiguous {len(amb)}")
for r in called:
    if r["end_agree"] != "True":
        a = adj.get(r["srx"], {})
        print(f"    disagree {r['srx']}: Arc {r['arc_end']}, us {r['our_end']} -> {a.get('category', 'UNADJUDICATED')}")
conf = Counter(r["confidence"] for r in called)
print(f"  confidence of calls: {dict(conf)}")
oos = [r for r in usable if r not in gex]
high = [r for r in oos if r["confidence"] == "HIGH"]
print(f"\nOut of scope by Arc's label ({len(oos)}): HIGH calls {len(high)}")
for r in high:
    a = adj.get(r["srx"], {})
    print(f"    {r['srx']} Arc {r['tech_10x']}, us {r['verdict']} -> {a.get('category', 'UNADJUDICATED')}: {a.get('true_label', '')}")
