# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Scores the head-to-head per the frozen plan in selection.tsv.

Usage: python score.py <results.tsv> <adjudication.tsv>
"""

import csv
import os
import statistics
import sys
from collections import Counter

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "arc", "scripts"))
from run_bench import VERDICT_MAP  # noqa: E402

LISTS = {"737K-august-2016", "3M-february-2018", "3M-3pgex-may-2023", "3M-5pgex-jan-2023"}
SANGER = {
    "10x_3p_v2_or_5p_v1_v2": ("737K-august-2016", "unresolved"),
    "10x_3p_v3_family": ("3M-february-2018", "3'"),
    "10x_3p_v4_gemx": ("3M-3pgex-may-2023", "3'"),
    "10x_5p_v3_gemx": ("3M-5pgex-jan-2023", "5'"),
}


def ours(r):
    v = r.get("ours_verdict") or ""
    if v.startswith("AMBIGUOUS between 10x"):
        return "737K-august-2016", "unresolved"
    return VERDICT_MAP.get(v, ("", ""))


def sanger(r):
    return SANGER.get(r.get("sanger_id") or "", ("", ""))


def pct(n, d):
    return f"{n}/{d} ({100 * n / d:.1f}%)" if d else "0/0"


def main(results, adjudication):
    rows = {}
    for r in csv.DictReader(open(results), delimiter="\t"):
        rows[r["srx"]] = r  # last attempt wins
    adj = {r["srx"]: r for r in csv.DictReader(open(adjudication), delimiter="\t")}
    usable = [r for r in rows.values() if r.get("condition") == "ok"]
    print(f"samples {len(rows)}, usable {len(usable)}; source {dict(Counter(r['source'] for r in usable))}")
    print(f"unusable: {dict(Counter(r['condition'][:40] for r in rows.values() if r.get('condition') != 'ok'))}")

    def in_scope(r):
        a = adj.get(r["srx"], {})
        if r["arc_whitelist"] not in LISTS:
            return False
        return r["tech_10x"] in ("3_prime_gex", "5_prime_gex") or a.get("lib_type") == "gex"

    gex = [r for r in usable if in_scope(r)]
    print(f"\nIN-SCOPE 10x GEX: {len(gex)} (incl. {sum(r['tech_10x'] not in ('3_prime_gex', '5_prime_gex') for r in gex)} "
          "Arc labels out of scope but ENA titles say gene expression)")
    for name, f in (("tussl-chemistry", ours), ("Sanger infer_10x", sanger)):
        wl_ok = sum(f(r)[0] == r["arc_whitelist"] for r in gex)
        resolved = [r for r in gex if f(r)[1] in ("3'", "5'")]
        determined = [r for r in resolved if adj.get(r["srx"], {}).get("truth_end", "undetermined") != "undetermined"]
        right = [r for r in determined if f(r)[1] == adj[r["srx"]]["truth_end"]]
        no_call = sum(f(r)[0] == "" for r in gex)
        print(f"  {name:17} whitelist = Arc {pct(wl_ok, len(gex))}; 3'/5' resolved {pct(len(resolved), len(gex))}; "
              f"correct vs adjudicated truth {pct(len(right), len(determined))}; no call {pct(no_call, len(gex))}")
        for r in determined:
            if f(r)[1] != adj[r["srx"]]["truth_end"]:
                a = adj[r["srx"]]
                print(f"      wrong end {r['srx']}: called {f(r)[1]}, truth {a['truth_end']} "
                      f"(title {a['title_end'] or '-'}, nuclear {a['nuc_sense']}/{a['nuc_anti']}) {a['title'][:70]}")
        for r in gex:
            if f(r)[0] and f(r)[0] != r["arc_whitelist"]:
                print(f"      wrong list {r['srx']}: called {f(r)[0]}, Arc {r['arc_whitelist']}")

    arc_vs_truth = [r for r in gex if r["tech_10x"] in ("3_prime_gex", "5_prime_gex")
                    and adj.get(r["srx"], {}).get("truth_end", "undetermined") != "undetermined"]
    arc_right = sum(("3'" if r["tech_10x"] == "3_prime_gex" else "5'") == adj[r["srx"]]["truth_end"] for r in arc_vs_truth)
    print(f"  {'Arc metadata label':17} 3'/5' correct vs adjudicated truth {pct(arc_right, len(arc_vs_truth))}")
    print(f"  truth determinable for {pct(sum(adj.get(r['srx'], {}).get('truth_end', 'undetermined') != 'undetermined' for r in gex), len(gex))}")

    high = [r for r in gex if r.get("ours_conf") == "HIGH"]
    high_wrong = [r for r in high if ours(r)[0] != r["arc_whitelist"] or (
        adj.get(r["srx"], {}).get("truth_end", "undetermined") not in ("undetermined", ours(r)[1]))]
    print(f"  tussl-chemistry HIGH calls {len(high)}, HIGH and wrong {len(high_wrong)}: {[r['srx'] for r in high_wrong]}")

    oos = [r for r in usable if r not in gex]
    print(f"\nOUT OF SCOPE (not 10x GEX by Arc's list or by title): {len(oos)}")
    print(f"  tussl-chemistry HIGH 10x GEX calls: {sum(r.get('ours_conf') == 'HIGH' and ours(r)[0] in LISTS for r in oos)}; "
          f"Sanger 10x calls: {sum(sanger(r)[0] in LISTS for r in oos)}")

    def cpu(key):
        v = [float(r[key]) for r in usable if r.get(key)]
        return f"median {statistics.median(v):.2f} s, mean {statistics.mean(v):.2f} s, max {max(v):.1f} s"
    print(f"\nCPU per sample: tussl-chemistry {cpu('ours_cpu')}; Sanger {cpu('sanger_cpu')}")


if __name__ == "__main__":
    main(*sys.argv[1:3])
