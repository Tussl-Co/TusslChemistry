# Validation

This file tells how tussl-chemistry earns trust in a call. It also tells how tchem shows users the
limits of that trust. The aim is to be scientific and transparent. We say what the evidence for
each kit is. We say where the tool can be wrong. We print the record as it is, not as we want it
to be.

[docs/glossary.md](../docs/glossary.md) defines each term in this file. Each term has one meaning.

Status: adopted 9 Oct 2026, with the design-study rule and the independent-label rule. The README
table ("How well it works") shows the level of each kit. `data/validation.tsv` is its source. This
file defines what that table means.

Trust has three layers. Each layer finds errors that the layer before it cannot find.

1. **The mechanism.** A HIGH call needs evidence that the reads cannot make by chance.
   When the reads cannot separate two kits, tchem says so. It does not guess.
2. **The tests.** The gate runs every layout on every commit, complete and damaged.
   A wrong HIGH call stops the commit.
   We record each pair of kits that the reads cannot separate.
   For each pair, a test fixes what tchem does.
3. **The held-out record.** This is a set of public runs that did not help to make a layout.
   We compare each call with what an independent source says the kit is.
   Every match prints the record, and the README shows it, at its current level.

## 1. The mechanism

tchem scores every known layout against a sample of reads (README, "How tchem decides"). The support
of a layout comes from three types of evidence:

- barcodes on the vendor's barcode list;
- fixed linkers at the positions that the layout gives;
- for kits with no public barcode list, how much the barcodes concentrate on few sequences.

| Verdict | Rule (`src/detect/decide.rs`, `Confidence`, `decide`) |
|---|---|
| HIGH | At least 70% support. Nothing in the reads contradicts the layout. No other layout has half as much support. The evidence includes a barcode list or a linker. Concentration alone never gives HIGH. A layout with no barcode list also needs each linker in at least 70% of reads (`UNLISTED_LINKER_HIGH`). |
| MEDIUM | A clear lead, but less than HIGH needs. |
| LOW | The barcode list agrees with the layout (80% of barcodes listed), but another part of the layout does not fit. The report shows the contradictions. |
| Ambiguous | Two or more layouts share their barcode lists and linkers and fit equally well. tchem tests the tail after the UMI, then the strand of the cDNA read. If neither decides, the verdict is ambiguous, and the report says why. |
| No match | No layout has the minimum support. |

Every rate is a share of all reads in the sample. A read that is too short to reach a segment
counts as a read without that segment. Thus short reads of a different kit cannot leave the count
(`segment_evidence_judge`).

Most of the safety comes from this mechanism. A barcode list match is almost all or nothing. A
10x 3' v3 run matches the v3 list at about 97% and the other 10x lists at about 1%. Errors come
less from chance than from cases that the mechanism was not made for. Layer 2 tests those cases.

## 2. The tests

`scripts/gate.sh` runs two suites before every commit:

- `scripts/validate.sh`, with one row for each kit in `scripts/validate_kits.tsv`. It runs tchem
  on regression inputs and compares each verdict word for word.
- `bench/chaos/scripts/chaos.py`. It damages one chaos base for each layout in 17 ways. Examples
  are a trimmed or missing barcode read, swapped files, a different kit mixed in, N bases,
  substitutions and damaged gzip. For each chaos case, tchem must not abort or give a wrong HIGH
  call, and it must give the same report two times. A test (`tests/tidy.rs`,
  `chaos_covers_every_layout`) fails if a layout has no chaos base.

The chaos bases and the tests use design runs, vendor demonstration runs, or runs outside every
held-out set. We tune the code against each test that the gate runs on every commit. Thus such a
test uses up runs from the held-out pool (section 7, open item 1).

### What the reads cannot separate

| Kits | Why | What tchem does | Fixed by |
|---|---|---|---|
| 10x 3' v2 and 10x 5' v1/v2 | one barcode list (737K) and the same 16+10 geometry | the tail after the UMI, then the strand of the cDNA read; if neither decides, ambiguous with the reason (`decide.rs`, `resolve_tie`) | `validate.sh`, "pbmc_1k_v2, no R2" |
| 10x 3' v3, v3.1, LT and HT | one barcode list and one geometry | one layout, "10x 3' v3/v3.1" | the name of the layout |
| Parse Evercode WT versions | the version that a lab names does not set the layout: early v3 kits read like v2 | two layouts by structure; held-out runs scored by family (`bench/heldout/scripts/run.py`, `EXPECTED`) | `validate_kits.tsv`, both Parse rows |
| Parse kit sizes (Mini, WT, Mega) | the same reads | not reported | none needed |
| Singleron GEXSCOPE MicroBead and Drop-seq | the same structure | called Drop-seq (2d0140e) | none yet |
| Drop-seq and all bead kits with unlisted barcodes | no public barcode list, so only barcode concentration | never HIGH, because concentration alone cannot give HIGH | `validate.sh`, "Drop-seq SRR1853178" (MEDIUM) |

## 3. The held-out record

### The question it answers

When tchem says "this is kit X" at HIGH confidence, how frequently is it wrong?

The record measures the precision of a HIGH call for one layout and one read type (short or long).
It does not measure recall, that is, how frequently tchem calls kit X when the reads are kit X. It
gives no promise for MEDIUM or LOW calls.

### Terms

[docs/glossary.md](../docs/glossary.md) defines run, study, label, exact label, family label, design
run, examined study, held-out run, right, wrong and disputed. Three details apply to this file:

- All demonstration runs of one vendor are one study, `vendor:<name>` (`sources.vendor_cases`).
- `examined.tsv` lists the examined studies of each bench.
- Adjudication can make a call right or disputed (section 5, step 5). Without adjudication, a HIGH
  call that the label does not permit is wrong.

### Levels

Each layout has one level. Every match prints the level (`layout.validation`). The README shows the
level of each layout with its counts. The levels are a report, not targets.

| Level | Rule | What it claims |
|---|---|---|
| `unvalidated` | no held-out study with an exact label | nothing yet |
| `checked` | 1 to 29 exact studies, none wrong | the counts, and no rate |
| `validated` | 30 or more exact studies, none wrong | a HIGH call's error rate is under 10% (95% confidence) |
| `error_under_3pct` | 100 or more exact studies, none wrong | a HIGH call's error rate is under 3% (95% confidence) |
| `contradicted` | any held-out HIGH call wrong | the layout's HIGH calls have failed |

Source: `src/layouts/validation.rs`, `Validation::level`, `STUDIES_VALIDATED`,
`STUDIES_ERROR_UNDER_3PCT`. A wrong call has priority over every count.

### The family record (adopted 9 Oct 2026)

Four kit families have more than one layout. These are BD Rhapsody (three bead layouts), Parse
Evercode WT (two layouts), SeekGene SeekOne DD (3' and 5') and Singleron GEXSCOPE (three versions).
tchem separates these layouts from the reads. But labs rarely record the bead type or the version.
Thus most runs have a family label, and the exact record above cannot count them.

The family record asks a question that a family label can answer. Of the held-out HIGH calls that
name a layout of the family, how many runs are of that family? The rules are those of the exact
record: independent labels only, studies counted, the same levels. There is one difference. A call
that names the wrong layout in the correct family is right in the family record. The exact record
keeps it as disputed or wrong.

Two lists name the families: `src/layouts/validation.rs` `FAMILIES` and
`bench/coverage/scripts/tally.py` `FAMILIES`. A test makes sure that they agree. The record is
`data/validation_family.tsv`. Every match of a family layout prints it (`layout.family_validation`,
schema 1.1), and the README shows it.

**Why a bound of 3 / N.** Zero errors in N trials does not show an error rate of zero. By the rule
of three, the true rate is less than 3 / N at 95% confidence. If the rate were 3 / N, N trials
without an error would occur only about 5% of the time. At 30 and 100 studies, that bound goes
below 10% and 3% (`Validation::error_bound`).

**Why studies, not runs.** The rule of three needs independent trials. Runs from one study share a
lab, a protocol, a kit lot, a sequencer and frequently one library. If tchem calls one of them
right, it almost certainly calls the others right. Thus the level counts studies, and the record
shows runs next to them.

**Why held-out runs.** We make a layout from the reads of real runs. A score on those runs is like
a test of a student with the exam that the student studied.

## 4. Labels

A label is only as good as its source (`referee` in `bench/coverage/results/calls.tsv`):

| Referee | Independent of the reads? |
|---|---|
| The run's library protocol text in the archive (catalogue `evidence`) | yes |
| The lab's own record: the run's ENA fields and GEO sample, mapped by `label10x.py` (`bench/head_to_head/LAB_LABELS.md`) | yes |
| The vendor's statement on a demo dataset page | yes |
| A BD part number tied to one bead type by BD's own documents (`bench/bd/BEAD_PARTS.md`) | yes |
| A GEO sample record or the run's title (adjudication) | yes |
| Arc Institute's scBaseCount list, recovered by matching reads to 10x lists | **no**: the same kind of signal tchem uses |
| Linker positions in the run's own reads (Singleron, 2d0140e; retired 9 Oct 2026, those runs are now catalogue rows labelled by protocol) | **no**: read off the reads |
| Arc's metadata label (`tech_10x`, an LLM reading metadata; one adjudication) | **no**: treated as not independent |

A right call against Arc's list shows that two list matchers agree. It does not show that tchem
agrees with the lab's record. Before the independent-label rule, 10x 3' v3 was `validated` on 97
studies. Arc's list was the only label for 96 of them (coverage at 0f33fdb). With the
labs' own records, 10x 3' v3 was `validated` on 54 independent studies (coverage at d6d44ef).

**The independent-label rule** (adopted 9 Oct 2026): the level counts independent labels only. A
HIGH call that is right only against a read-derived label goes in its own column, "read-derived
label only". The level does not count it (`tally.READ_DERIVED`). A wrong call counts, whatever its
referee.

## 5. The procedure for a held-out record

1. **Write the selection rules and the labelling rules first.** These rules include:
   - the ENA query;
   - the classifier that names the kit of a run from its protocol text;
   - the read-length filter;
   - whether the record scores the kit by exact label or by family (`select.py`;
     `bench/heldout/scripts/run.py`, `EXPECTED`).

   Commit the rules before any draw. You can revise a rule after a draw and before any run.
   Record the revision and its reason (Flex and Singleron, 9 Oct 2026). After tchem calls runs,
   a revised rule applies only to studies not yet seen. It never relabels the runs whose calls
   showed the fault.
2. **Exclude what made the layout.** Design runs and examined studies are never held-out.
   **The design-study rule** (adopted 9 Oct 2026) adds this:
   - A study whose reads made a layout no longer tests that layout, the same as an examined study.
   - Its other runs share its lab, protocol and frequently its library, so they are not
     independent of it.
   - Coverage applies the rule when it selects its cases (`sources.design_studies`,
     `sources.shaped_by_design`), for each layout and read type.
   - At its adoption, the rule removed 17 runs.
3. **Draw and commit the catalogue before any run.** Use metadata and read lengths only.
   Never use read content (`catalogue.draw`, `catalogue.read_length`).
   Every row has a `split` and an `evidence` value.
4. **Run a frozen binary one time on every held-out case.** Use `bench/coverage/scripts/run.py`.
   - The coverage commit must be in the released history (`scripts/coverage_fresh.sh`).
   - A case that a binary with the same detection code called before keeps that call (`--reuse`).
     The call depends only on the reads and the detection code. Each row names its commit.
   - Only the new cases run. A seeded random selection of 200 carried cases runs again and must agree.
   - Thus we can delete the read heads of scored runs on `error_under_3pct` layouts to free disk.
   - When the detection code changes, no call is carried, and we download the heads again.
5. **Adjudicate disagreements by the rule fixed before the run.** Record each decision in
   `bench/coverage/results/adjudication.tsv`:
   - `right` when an independent source names the called kit;
   - `disputed` when a source names a different kit and only the reads' evidence supports the
     call. The evidence is the barcode list or, where two kits share a list, the strand;
   - otherwise `wrong`.

   Quote the evidence.
6. **Count studies and render the README table.** Use `bench/coverage/scripts/tally.py`.
   Render recall from the same coverage pass (`bench/coverage/scripts/recall.py`).
   Recall shows what the binary said on the runs with each kit label. No level counts recall.
7. **The design-data rule.** A held-out run that we use to change the code becomes a design run.
   This occurs before the next coverage pass.
   - A `contradicted` layout gets a new level only from studies not yet seen.
   - `bench/coverage/results/rule7.tsv` lists such runs, with the reason.
   - Coverage removes their full study, as the design-study rule does for a design run.

Each bench keeps these records:

- `scripts/select.py`;
- `results/catalogue.tsv`;
- `results/examined.tsv` (if it is absent, nothing is recorded; this does not mean that nothing
  was read);
- `RESEARCH.md`, which names every examined run and design run.

`bench/coverage/scripts/inventory.py` compares these records. It reports a status with no record
behind it as unknown.

## 6. What the record does not claim

- That the studies are typical. Catalogues come from archive studies whose protocol names the kit.
- That a label is correct. Adjudication finds some wrong labels, not all.
- Anything about the exact version from a call with a family label.
- That every kit can get to every level. A level needs public studies with exact labels. For some
  kits, the archives do not have sufficient studies. The README then shows the level that the
  archives support.

## 7. Decided, and open

Decided 9 Oct 2026:

- the design-study rule;
- the independent-label rule;
- the aim for each kit. Every layout has a chaos base and a test. Every layout has a
  held-out record from independent labels, at the highest level that the public archives permit.
  The README shows the level and the label sources.

Open:

1. **Held-out runs in the regression tests.** The header of `scripts/validate_kits.tsv` sets a rule.
   - The rule: its tests use design runs and vendor data, never a held-out set.
   - Its three Singleron rows and the three Singleron chaos bases are held-out rows of the
     `bench/singleron` catalogue.
   - The chaos bases for 10x 3' v4 and 5' v3 (`10k_hgmm_3p_gemx_Multiplex`,
     `10k_hgmm_5p_gemx_Multiplex`) are held-out 10x demonstration runs. Neither layout has a
     recorded design run.
   - We tune the code against each test that the gate runs on every commit. Thus those runs must
     become design runs, or the rule in the header must change.
   - Its 10x Multiome and Flex rows were held-out runs too. The design-study rule removed them.
2. **Labelling rules that the first independent coverage pass found at fault.**
   The source is `adjudication.tsv`, 9 Oct 2026.
   We fix them for the next draw and apply them only to studies not yet seen:
   - "GEM-X Flex" read as Flex v2. 10x also uses the name for the v1 kit (34 runs disputed).
   - A study that names two assays, catalogued as one assay (14 runs).
   - "ver.3" not read as a version by `label10x.py` (6 runs, adjudicated right).
3. **The Singleron kit-version assumption.** A kit "V2" next to GEXSCOPE is taken as the V2 layout.
   No held-out call contradicts it yet (8 studies V2, 2 V3).
