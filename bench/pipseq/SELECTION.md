# Revised PIPseq draw (9 Oct 2026)

`select.py` is bench/pipseq/scripts/select.py with only the version token revised. The
docstring records the revision, its reason and the date. The draw uses today's
bench/heldout/scripts/catalogue.py, unchanged. Nothing here is applied to the repository.

```sh
R=<repository root>
cd <scratch>/exec2/pipseq/draw
python3 -I -B run_draw.py $R .
```

`run_draw.py` calls `catalogue.draw(QUERY, classify, shape)` from the revised select.py. It wraps
`ena_rows`, `classify`, `shape` and `read_length` only to observe them, and their results are
unchanged. It then compares the draw's design rows with those in
bench/pipseq/results/catalogue.tsv. If they differ, it stops before it writes `catalogue.tsv`.
They agree: v3 PRJEB74471/ERR12826674, v4 PRJNA1014494/SRR25991440 and
v5 PRJNA1251639/SRR33161101.

Outputs:
- `catalogue.tsv`: catalogue.OUT_FIELDS, as `catalogue.write` writes it. 76 runs from 39 studies.
- `refused.tsv`: each study that the query returns but the draw does not include. Columns are
  study, reason, detail and protocol. The reason is the first rule that the study fails:
  `no_paired_run`, `other_vendor`, `no_single_version` or `shape_refused`.

Reads: as the shared rule requires, the draw reads modal read lengths only (the first 64 KB of
each file) on up to 6 runs per classified study. It reads no other sequence.

# Second extension, 9 Oct 2026: PIPseq V by name

`scripts/select.py` was revised a second time, before any added study was run. The revision and
its reason are in the docstring. A study that states no version token but names PIPseq V, or
Illumina's Single Cell 3' RNA Prep, is drawn as `pipseq_v5`. The redraw left every row of a study
already in the catalogue unchanged. Ten held-out rows from five new v5 studies were appended:
PRJNA1284893, PRJNA1377937, PRJNA1445686, PRJNA1477052 and PRJNA1480533.

The PIPseq V research (`research_v5/studies_table.md`) also lists four single-kit studies that this
draw does not reach: PRJNA1268916, PRJNA1269311, PRJNA1438756 and PRJNA1481985. They were found
through GEO, and ENA's query or the shape rule does not return them. They are left out, not
chased one by one.

PRJNA1216192 and PRJNA1373237 were drawn on 8 Oct as `pipseq_v4`, but their protocols also name
PIPseq V or the Illumina kit. They stay as drawn, because the catalogue is append-only. A
disagreement on them is adjudicated with this on record.
