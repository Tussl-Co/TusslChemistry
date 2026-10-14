# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""The evidence for adjudicating the HIGH calls a tally marks wrong (bench/VALIDATION.md, step 5):
for each, the label and where it comes from, and the reads' own barcode-list evidence from the
frozen binary's report (`--all`): the hit rate of the list the call names, and the next two
layouts' best list rates. Adjudication is then decided by the fixed rule, once, by a person.

Usage: python evidence.py <tally output> <tchem> <whitelist_dir> <data_dir> <out.tsv>
  <tally output> is what tally.py printed (its WRONG lines).
"""

import collections
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", ".."))
import sources  # noqa: E402
from invariant import invariant  # noqa: E402

BENCH = sources.BENCH
WRONG = re.compile(r"WRONG  (.+?) on (\S+) (\S+): truth (.*)")
CANDIDATE = re.compile(r"  (\S.*?)\s+on R[12]\s+support")
LIST_RATE = re.compile(r"(\S+)@\d+ (\d+)%")
CATALOGUE_BENCHES = ("flex", "singleron", "bd", "mgi", "pipseq", "parse", "indrops", "seekgene", "tenx_lab", "tenx_ena")


def label_quotes():
    """run -> (where the label comes from, its quoted evidence)."""
    quotes = {}
    for name in CATALOGUE_BENCHES:
        path = os.path.join(BENCH, name, "results", "catalogue.tsv")
        for r in sources.rows(path):
            quotes[r["run"]] = (f"catalogue {name}", r.get("evidence", ""))
    for r in sources.rows(os.path.join(BENCH, "bd", "results", "exact_labels.tsv")):
        quotes[r["run"]] = ("BD part", r["evidence"])
    statements = collections.defaultdict(list)
    for r in sources.rows(os.path.join(BENCH, "head_to_head", "results",
                                       "lab_label_statements.tsv")):
        if r["mapped"]:
            statements[r["srx"]].append(f"[{r['field']}|{r['mapped']}] {r['quoted_sentence'][:160]}")
    for srx, said in statements.items():
        quotes.setdefault(srx, ("lab record", " || ".join(said)))
    invariant(quotes, "labels have quotes")
    return quotes


def list_rates(tool, whitelists, case):
    """layout -> its best barcode-list hit rate in the binary's --all report."""
    inputs = ["--r1", case.files[0], "--r2", case.files[1]] if len(case.files) == 2 else list(case.files)
    report = subprocess.run([tool, "--whitelists", whitelists, "--all", *inputs],
                            capture_output=True, text=True, timeout=600).stdout.split("\n")
    rates = {}
    for i, line in enumerate(report):
        m = CANDIDATE.match(line)
        if m and i + 1 < len(report):
            hits = [int(p) for _, p in LIST_RATE.findall(report[i + 1])]
            rates[m.group(1)] = max(hits) if hits else None
    return rates


def main(tally_output, tool, whitelists, data, out_path):
    wrong = WRONG.findall(open(tally_output).read())
    cases = {c.run: c for c in sources.cases(data)}
    quotes = label_quotes()
    with open(out_path, "w") as out:
        out.write("source\trun\tstudy\tcalled\ttruth\tcalled_list_rate\tnext_best\tlabel_evidence\n")
        for called, source, run, truth in wrong:
            invariant(run in cases, f"{run}: a wrong call on a run coverage does not hold")
            case = cases[run]
            rates = list_rates(tool, whitelists, case)
            others = sorted(((v, k) for k, v in rates.items() if k != called and v is not None),
                            reverse=True)[:2]
            where, quote = quotes.get(run, ("(no quote)", ""))
            out.write("\t".join([source, run, case.study, called, truth.strip("'\""),
                                 str(rates.get(called)), "; ".join(f"{k} {v}%" for v, k in others),
                                 f"{where}: {quote}"[:600].replace("\t", " ")]) + "\n")
    print(f"{len(wrong)} wrong calls with evidence -> {out_path}")


if __name__ == "__main__":
    if len(sys.argv) != 6:
        sys.exit(__doc__)
    main(*sys.argv[1:])
