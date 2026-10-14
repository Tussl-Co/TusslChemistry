# inDrops held-out catalogue (draft, 9 Oct 2026)

Drawn from ENA metadata and read lengths only (catalogue.draw's 64 KB length check). No read
content was examined, tchem was not run, and no results file was opened.

## Commands

    python3 -I select.py . ../examined_wild.tsv

`ena_probe.json` is a metadata-only fetch of the same query, made the same day before the draw.
It was used to test the classifier and to analyse the misfires below. All 141 drawn runs are in it.

## Files

- `select.py`: the query, classifier, shape rule and draw, with the rules in its docstring.
- `catalogue.tsv`: catalogue.OUT_FIELDS. Evidence is a verbatim 120-character window.
- `examined.tsv`: rows of examined_wild.tsv whose study was drawn (run, study, what).
- `refused.tsv`: study, run, reason, lengths, protocol_window. A row with a blank run refuses the
  whole study. A row with a run is a run that the shape rule refused; its `lengths` are the modal
  read lengths of its files in name order. A length of 0 would mean the file could not be read.
  There are none.
- `versions.tsv`: every inDrop version token seen in a drawn study's protocol.
- `*_v1.tsv`: the first draw, before revisions 1-3. `*_v2.tsv`: the second draw, before revision 2b.

## Counts (third draw, after revision 2b, 9 Oct 2026)

ENA returned 186 studies.

| kit | split | runs | studies |
|---|---|---|---|
| indrops | held_out | 128 | 65 |
| indrops | examined | 2 | 1 (PRJEB55757) |
| indrops_v3 | held_out | 1 | 1 (PRJNA736741) |
| total | | 131 | 67 |

No row is design (`design=False`). Nothing is labelled TruDrop separately: TruDrop keeps the V2 barcode read (see below).

Refusals, by study (blank run): not paired 83, other vendor 15, shape 9, other method named 4,
other method named (unclear) 3, no inDrop statement 2, no fastq 2, mixed versions 1. Total 119.
Refusals, by run: shape 45 runs. No run was unreadable (length 0).

## Revisions (9 Oct 2026, after the first draw and before any read was run)

The first draw is kept as `catalogue_v1.tsv` (141 runs, 72 studies), `refused_v1.tsv` and
`versions_v1.tsv`.

1. The inDrop word rule became `\bin-?drops?(?:tm|\d+)?(?![a-z])|1cellbio`. The old `\b` failed on
   "inDropTM barcoded beads" (PRJNA996596) and "an inDrops21 setup" (PRJNA679887). The version
   rule's inDrop mention uses the same pattern.
2. A protocol that names another method used beside inDrop is refused as "other method named".
   The pattern is `drop-?seq|cel-?seq2|\b10x\b(?! with)`. A protocol that names Smart-seq is refused as
   "other method named (unclear)". A bare "CEL-Seq" is inDrop chemistry and is not refused.
   Examples: PRJNA438523 ("followed inDrop, Drop-seq and 10X standart protocol respectively"),
   PRJNA473915 ("For CEL-Seq2 sample"), PRJNA493298 ("10x cells were prepared").
3. A protocol that states both v2 and v3 for inDrop is refused as "mixed versions". The only
   case is PRJNA471420 ("Most samples were run using 'inDrops v2' ... one experiment used
   'inDrops v3' beads").

Effect against v1: 12 rows from 6 studies are removed (PRJNA438523, 471420, 473915, 493298,
698490, 860187). No row is added, and no kept row changed kit or split. Revision 1 added nothing.
PRJNA679887 has no paired run. PRJNA996596 is now refused by revision 2 because of "connected to
the Drop-seq apparatus" (see the open item below).

2b. (9 Oct 2026, after the second draw and before any read was run) Revision 2 fired where no
   other method was used:
   - "10X Danieau Buffer" (PRJNA445487)
   - "InDrop library protocol, modified CEL-Seq2" (PRJNA497610)
   - "inDrop(TM) barcoded beads in a home-made microfluidic device connected to the Drop-seq
     apparatus" (PRJNA996596)

   The pattern became `drop-?seq(?! apparatus)|(?<!modified )cel-?seq2|\b10x\b(?!\s+with)(?!\s+(?:\w+\s+)?buffer)`.
   The second draw is kept as `catalogue_v2.tsv`, `refused_v2.tsv` and `versions_v2.tsv`.
   Effect against v2: 2 rows are added (PRJNA996596: SRR25865112, SRR25865113; indrops,
   held_out). None is removed or changed. PRJNA445487 and PRJNA497610 are now refused as "not
   paired". `examined.tsv` is byte-identical.

## TruDrop

Southard-Smith et al., BMC Genomics 2020, https://pmc.ncbi.nlm.nih.gov/articles/PMC7331155/
(doi 10.1186/s12864-020-06843-0). The Methods say: "The sequence for the final TruDrop library is
as follows: Cell Barcode + UMI (P5): ...[cell barcode 1] GAGTGATTGCTTGTGACGCCTT [cell barcode 2]
[UMI]..." and "The new TruDrop library structure utilizes the standard Illumina TruSeq sequencing
primers." TruDrop also "builds upon the widely-used, commercially available V2 gel beads". The
barcode read keeps V2's BC1 + W1 + BC2 + UMI + polyT. Only the adapters and the side change (P5,
TruSeq Read 1). TruDrop studies therefore stay `indrops`.

## Rule notes

- "Near inDrop" for v3 has two parts. The first is a version token within 60 characters after an
  inDrop or 1CellBio mention, in the same clause. The second is "version N" or "vN
  beads/hydrogel" within 200 characters of a mention, across sentences. The second part was
  added before the draw, after PRJNA414337 ("Hydrogel beads with version 3 oligonucleotide
  design"). It also catches PRJNA516696 ("library version v3 where R1 is the biological read").
  Neither study is in the catalogue (both are not paired).
- See the final report for the misfires found in this draw and the proposed revisions.
