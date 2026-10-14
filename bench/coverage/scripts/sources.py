# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

"""Every held-out run the benchmarks have, as one list of cases: where its reads are cached, which
layouts its truth admits, at what resolution, and who says so. The coverage run feeds these to the
current binary; the tally turns the calls into each layout's validation.

A case's truth is the set of layout names it may rightly be called. Resolution is "exact" when the
truth names one layout and "family" when it admits several (BD archive runs rarely state the bead
version; Parse's kit version does not fix its layout; a 737K list without a known 3'/5' end). An
empty truth means no layout is right (a kit the tool does not claim): any confident call is wrong.

Design runs (those that settled a layout) are left out, from each catalogue's `split` column or,
for the vendor gold set and the older layouts, from DESIGN below, which cites where each was named.
"""

import collections
import csv
import glob
import importlib.util
import os
import re
import sys
from dataclasses import dataclass

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from invariant import invariant  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.dirname(os.path.dirname(HERE))


def _load(name, path):
    """A bench script as a module, by path: several benches name their runner run.py."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


# The held-out catalogues' kit -> layouts, and the vendor gold wording -> layout: one source each.
EXPECTED = _load("heldout_run", os.path.join(BENCH, "heldout", "scripts", "run.py")).EXPECTED
_score = _load("vendor_score", os.path.join(BENCH, "vendor_gold", "scripts", "score.py"))
UNSUPPORTED, vendor_expected = _score.UNSUPPORTED, _score.expected

V2, V5 = "10x 3' v2", "10x 5' v1/v2"
# Catalogue kits whose metadata names a kit the tool does not claim: no layout is right for them.
NOT_CLAIMED = {"10x_gemx_flex_v2", "indrops_v3", "seekgene_mm"}
# Runs that settled a layout, outside any `split` column, with where each was named.
DESIGN = {
    "pbmc_granulocyte_sorted_10k": "Multiome GEX design run (81da1b8)",
    "10k_K562_singleplex_Multiplex": "Flex v1 design run (b81d304)",
    "ScaleBio_test_202506_quantumV1": "QuantumScale design run (c075d12)",
    "BD_WTANext_3070": "BD Enhanced V2/V3: WTA Next demos read before the layout (b03a942)",
    "BD_WTANext_AbSeq_3070": "BD Enhanced V2/V3: WTA Next demos read before the layout (b03a942)",
    "BD_WTANext_TCRBCR_D40_3074": "BD Enhanced V2/V3: WTA Next demos read before the layout (b03a942)",
}
# Long-read heads end mid-gzip (fetched by byte range); the long-read bench reads 2,000 reads,
# well inside each, and so does the coverage run.
LONG_READS = 2000
# The 10x lists Arc's pipeline recovered, as layouts (head-to-head truth). A 737K list is 3' v2 or
# 5' v1/v2: the end decides, from the run's title when it names one.
ARC_LIST = {
    "3M-february-2018": ("10x 3' v3/v3.1",),
    "3M-3pgex-may-2023": ("10x 3' v4 (GEM-X)",),
    "3M-5pgex-jan-2023": ("10x 5' v3 (GEM-X)",),
    "737K-arc-v1": ("10x Multiome GEX",),
}
LONG_KIT = [  # long-read catalogue wording -> layout, first match wins
    (r"^3prime v2", "10x 3' v2"), (r"^3prime v3", "10x 3' v3/v3.1"), (r"^3prime v4", "10x 3' v4 (GEM-X)"),
    (r"^5prime v(1|2)", V5), (r"^5prime v3", "10x 5' v3 (GEM-X)"),
]


@dataclass(frozen=True)
class Case:
    source: str
    run: str  # run accession, or the vendor's dataset name
    study: str  # empty when unknown: then the run counts as its own study
    truth: tuple  # layout names this run may rightly be called
    referee: str  # who says what the truth is
    read_kind: str  # short or long
    files: tuple  # cached read files, in the order to give tchem
    also: tuple = ()  # other accessions of the same run (an SRX's SRR), for de-duplication
    reads: int = 0  # --reads to pass; 0 for the tool's default

    @property
    def resolution(self):
        return {0: "unclaimed", 1: "exact"}.get(len(self.truth), "family")


def rows(path):
    with open(path, newline="") as f:
        yield from csv.DictReader((line for line in f if not line.startswith("#")), delimiter="\t")


# A catalogue whose labels are not the run's protocol alone, with its referee.
REFEREE = {"tenx_lab": "the lab's own record (bench/tenx_lab/results/statements_drawn.tsv.gz)"}


def heldout_cases(data, name, catalogue):
    """A catalogue in the held-out format (bench/heldout): kit, study, run, split. A run listed in
    `exact_labels.tsv` beside the catalogue takes that file's layout and referee instead of its
    kit's (BD runs whose protocol cites a bead-specific part, bench/bd/BEAD_PARTS.md)."""
    heads = os.path.join(data, f"{name}_heldout", "heads")
    exact_path = os.path.join(os.path.dirname(catalogue), "exact_labels.tsv")
    exact = {r["run"]: r for r in rows(exact_path)} if os.path.exists(exact_path) else {}
    for r in rows(catalogue):
        # Only rows marked held out count: design rows shaped the layout, and rows of studies whose
        # reads were examined while the kit was researched are not independent of it.
        if r.get("split", "held_out") != "held_out":
            continue
        truth = () if r["kit"] in NOT_CLAIMED else EXPECTED[r["kit"]]
        referee = REFEREE.get(name, "the run's own protocol (catalogue)")
        if r["run"] in exact:
            invariant(exact[r["run"]]["layout"] in truth, f"{r['run']}: exact label outside its kit")
            truth, referee = (exact[r["run"]]["layout"],), exact[r["run"]]["referee"]
        files = tuple(os.path.join(heads, f"{r['run']}_R{i}.fq.gz") for i in (1, 2))
        yield Case(name, r["run"], r["study"], truth, referee, "short", files)


def vendor_cases(data, results):
    heads = os.path.join(data, "vendor_gold", "heads")
    for r in rows(results):
        if r["dataset"] in DESIGN or r["status"] != "ok":
            continue
        want = vendor_expected(r["vendor"], r["chemistry"], r["dataset"])
        if want is None:
            continue
        truth = () if want.startswith(UNSUPPORTED) else (want,)
        files = tuple(os.path.join(heads, f"{r['dataset']}_R{i}.fq.gz") for i in (1, 2))
        # One vendor's demos share one source and one processing: they count as one study.
        yield Case("vendor_gold", r["dataset"], f"vendor:{r['vendor']}", truth,
                   "the vendor's stated chemistry", "short", files)


# A label the lab's own record states (bench/head_to_head/LAB_LABELS.md): independent of the reads.
LAB_REFEREE = "the lab's own record (bench/head_to_head/results/lab_labels.tsv)"


def lab_labels(bench):
    """srx -> the one layout its own records state, for runs whose records state exactly one
    layout tchem has; partial, conflicting and unstated labels are left out."""
    out = {}
    for r in rows(f"{bench}/head_to_head/results/lab_labels.tsv"):
        if r["label"].startswith("10x "):
            out[r["srx"]] = r["label"]
    invariant(out, "lab labels are recorded")
    return out


def head_to_head_cases(data, bench):
    studies = {r["srx"]: r["study"] for r in rows(f"{bench}/head_to_head/results/studies.tsv")}
    stated = lab_labels(bench)
    ends = {r["srx"]: r for r in rows(f"{bench}/head_to_head/results/adjudication_b.tsv")}
    for r in rows(f"{bench}/head_to_head/results/results.tsv"):
        pair = [p for p in r["ours_pair"].split(",") if p]
        if r["condition"] != "ok" or not pair or r["arc_whitelist"] not in (*ARC_LIST, "737K-august-2016"):
            continue
        cache = os.path.join(data, "h2h", "cache")
        if r["source"] == "SRA":
            cache = os.path.join(cache, f"sra_{r['srr']}")
        title_end = ends.get(r["srx"], {}).get("title_end", "")
        if r["srx"] in stated:
            truth, referee = (stated[r["srx"]],), LAB_REFEREE
        elif r["arc_whitelist"] in ARC_LIST:
            truth, referee = ARC_LIST[r["arc_whitelist"]], "Arc's recovered list"
        elif title_end in ("3'", "5'"):
            truth, referee = ((V2,) if title_end == "3'" else (V5,)), "Arc's list; end from title"
        else:
            truth, referee = (V2, V5), "Arc's list; end unknown"
        yield Case("head_to_head", r["srx"], studies.get(r["srx"], ""), truth, referee, "short",
                   tuple(os.path.join(cache, p) for p in pair[:2]), also=(r["srr"],))


def long_read_cases(data, bench):
    heads = os.path.join(data, "longread", "heads_heldout")
    for r in rows(f"{bench}/long_read/results/catalogue_heldout.tsv"):
        kit = r["kit"]
        if kit.startswith("uncertain"):
            continue
        truth = () if kit.startswith("negative") else tuple(
            layout for pattern, layout in LONG_KIT if re.match(pattern, kit))[:1]
        invariant(kit.startswith("negative") or truth, f"{r['run']}: unmapped long-read kit {kit!r}")
        found = [f for f in (f"{r['run']}.fq.gz", f"{r['run']}.fq", f"{r['run']}.bam")
                 if os.path.exists(os.path.join(heads, f))]
        files = (os.path.join(heads, found[0]),) if found else (os.path.join(heads, r["run"]),)
        study = re.split(r"[ :/]", r["study"], maxsplit=1)[0]
        yield Case("long_read", r["run"], study, truth, "the run's own protocol (catalogue)",
                   "long", files, reads=LONG_READS)


def small_cases(data):
    """Runs held out by name in a commit, with no catalogue of their own. (Singleron's three, held
    out by their reads' linker positions in 2d0140e, are now rows of bench/singleron's catalogue,
    labelled from their protocols, 9 Oct 2026.)"""
    arc = os.path.join(data, "arc_bench")
    for srx, study, files in [
        ("SRX15432977", "PRJNA841957", ("cache/SRR19374511_1.fastq.gz", "cache/SRR19374511_2.fastq.gz")),
        ("SRX18306640", "PRJNA903152", ("cache/SRR22333821_1.fastq.gz", "cache/SRR22333821_2.fastq.gz")),
        ("SRX22021743", "PRJNA1024986", ("cache/SRR26313439_1.fastq.gz", "cache/SRR26313439_2.fastq.gz")),
        ("SRX27291582", "PRJNA1208088", ("sra/SRR31934508/SRR31934508_2.fastq.gz",
                                         "sra/SRR31934508/SRR31934508_3.fastq.gz")),
        ("SRX27325568", "PRJNA1194194", ("sra/SRR31970277/SRR31970277_3.fastq.gz",
                                         "sra/SRR31970277/SRR31970277_4.fastq.gz")),
    ]:
        yield Case("multiome_archive", srx, study, ("10x Multiome GEX",),
                   "archive protocol names the kit (81da1b8)", "short",
                   tuple(os.path.join(arc, f) for f in files))


def design_studies(bench=BENCH):
    """(study, layout, read kind) -> the design run that shaped that layout from that study.

    Rule 2 (bench/VALIDATION.md, adopted 9 Oct 2026): a study whose reads shaped a layout tests
    that layout no more. Runs of one study share a lab, a protocol and often a library, so a
    design run's sister runs are not independent of it. Sources: every catalogue's `split=design`
    rows, the long-read design catalogue, and the vendor demos in DESIGN (all of one vendor's
    demos are one study)."""
    out = {}
    for path in sorted(glob.glob(f"{bench}/*/results/catalogue.tsv")):
        name = path.split(os.sep)[-3]
        if name == "long_read":
            continue  # its design catalogue has no split column: every row is design (below)
        for r in rows(path):
            if r.get("split") == "design":
                for layout in EXPECTED.get(r["kit"], ()):
                    out.setdefault((r["study"], layout, "short"), f"{name} design run {r['run']}")
    for r in rows(f"{bench}/long_read/results/catalogue.tsv"):
        study = re.split(r"[ :/]", r["study"], maxsplit=1)[0]
        for pattern, layout in LONG_KIT:
            if re.match(pattern, r["kit"]):
                out.setdefault((study, layout, "long"), f"long_read design run {r['run']}")
                break
    vendor = {r["dataset"]: r for r in rows(f"{bench}/vendor_gold/results/results_quantumscale.tsv")}
    for dataset in DESIGN:
        r = vendor.get(dataset)
        invariant(r is not None, f"DESIGN names {dataset}, not in the vendor gold results")
        want = vendor_expected(r["vendor"], r["chemistry"], dataset)
        invariant(want, f"DESIGN run {dataset} has no expected layout")
        out.setdefault((f"vendor:{r['vendor']}", want, "short"), f"vendor design demo {dataset}")
    invariant(out, "design runs are recorded")
    return out


def shaped_by_design(found, designs):
    """The cases `found` whose study shaped one of the layouts they may rightly be called, each with
    the design run that did (rule 2)."""
    out = []
    for case in found:
        why = next((designs[(case.study, layout, case.read_kind)] for layout in case.truth
                    if (case.study, layout, case.read_kind) in designs), None)
        if why is not None:
            out.append((case, why))
    invariant(all(c.study for c, _ in out), "a case without a study is its own study, never shared")
    return out


def drawn(data, bench=BENCH):
    """Every case the sources hold, each run once: a run in two sources is kept from the first
    listed. Rule 2 is not yet applied (`cases`)."""
    out, seen = [], set()
    generators = [
        *(heldout_cases(data, n, f"{bench}/{n}/results/catalogue.tsv")
          for n in ("flex", "bd", "mgi", "parse", "pipseq", "singleron", "indrops", "seekgene", "tenx_lab", "tenx_ena")),
        small_cases(data),
        long_read_cases(data, bench),
        head_to_head_cases(data, bench),
        vendor_cases(data, f"{bench}/vendor_gold/results/results_quantumscale.tsv"),
    ]
    for generator in generators:
        for case in generator:
            names = {case.run, *case.also}
            if names & seen:
                continue
            seen |= names
            out.append(case)
    invariant(len({c.run for c in out}) == len(out), "each run appears once")
    return out


def rule7_studies(bench=BENCH):
    """Studies of held-out runs that were used to change the code (rule 7, results/rule7.tsv): design
    data from then on, and with rule 2 their whole study leaves the held-out set."""
    path = f"{bench}/coverage/results/rule7.tsv"
    found = [r for r in rows(path)] if os.path.exists(path) else []
    invariant(all(r["run"] and r["study"] and r["why"] for r in found), "every rule 7 run says why")
    return {r["study"] for r in found}


def cases(data, bench=BENCH):
    """Every held-out case: those drawn, without those whose study shaped a layout they test
    (rule 2, design_studies) or holds a run used to change the code (rule 7, rule7_studies)."""
    out = drawn(data, bench)
    shaped = {c.run for c, _ in shaped_by_design(out, design_studies(bench))}
    shaped |= {c.run for c in out if c.study in rule7_studies(bench)}
    kept = [c for c in out if c.run not in shaped]
    invariant(len(kept) + len(shaped) == len(out), "each case kept or set aside once")
    return kept
