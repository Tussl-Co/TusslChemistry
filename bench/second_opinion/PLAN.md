# Second opinions on the 999 disputed HIGH calls: plan, fixed 10 Oct 2026 before any tool ran

Scope: every disputed HIGH call in bench/coverage/results/adjudication.tsv whose reads are on disk
(pruned heads re-fetched with fetch_heads.py). Inputs: the same head files coverage used.
None of these results changes a level: all are read-derived (A4); they are reported beside the record.

1. Sanger cellgeni infer_10x_chemistry.py @ a8611621 (nf-reprocessing-public-10x), --dry-run
   --allow-missing-whitelists, given the 10x lists it knows (737K-august-2016, 3M-february-2018,
   3M-3pgex-may-2023, 3M-5pgex-jan-2023, gex_737K-arc-v1), the same lists tchem uses.
   Scored on 10x-vs-10x disputes: Sanger agrees with tchem / agrees with the label / neither / no call.
   Its ids map: 10x_3p_v2_or_5p_v1_v2 -> {3' v2, 5' v1/v2}; 10x_3p_v3_family -> 3' v3;
   10x_3p_v4_gemx -> 3' v4; 10x_5p_v3_gemx -> 5' v3; 10x_multiome_arc_gex -> Multiome GEX.
   A 3'/5' split it does not resolve counts as matching both.
2. Strand referee (3' vs 5'), for disputes where label and call differ in end (3' vs 5'):
   N = strand of the cDNA read on spliced transcripts of 40 human genes disjoint from tchem's
   housekeeping set (any k-mer shared with data/reference/housekeeping_human.fa dropped), via
   bench/head_to_head/scripts/adjudicate.py strand() with its default bar (share >= 0.2%, ratio >= 3).
   Outcome: supports the call / supports the label / undetermined (includes non-human runs).
3. Cell Ranger 10.1.0 count --chemistry auto, human reference refdata-gex-GRCh38-2024-A, on a
   seeded sample (random.Random(20261010)) of human disputed 10x runs, one at a time.
   Scored as in 1 against its detected chemistry.
Set B genes: RPL4,RPL5,RPL7,RPL8,RPL11,RPL12,RPL18,RPL19,RPL23,RPL27A,RPL30,RPL32,RPL34,RPL35A,RPL37A,RPS3,RPS4X,RPS6,RPS8,RPS11,RPS13,RPS14,RPS15A,RPS19,RPS23,RPS24,RPS25,RPS27,RPS28,RPS29,EEF1G,HSP90AB1,PPIA,UBB,UBC,ATP5F1E,COX4I1,SERF2,TOMM7,NPM1

Amended 10 Oct 2026, before any tool ran: the disputes are 518 mouse, 273 human, 208 other or
unstated (ENA scientific_name, compare/species.tsv). The strand referee uses the human set B for
Homo sapiens runs and the same 40 genes from GENCODE M35 on GRCm39 for Mus musculus runs
(Rpl4,Rpl5,Rpl7,Rpl8,Rpl11,Rpl12,Rpl18,Rpl19,Rpl23,Rpl27a,Rpl30,Rpl32,Rpl34,Rpl35a,Rpl37a,Rps3,Rps4x,Rps6,Rps8,Rps11,Rps13,Rps14,Rps15a,Rps19,Rps23,Rps24,Rps25,Rps27,Rps28,Rps29,Eef1g,Hsp90ab1,Ppia,Ubb,Ubc,Atp5f1e,Cox4i1,Serf2,Tomm7,Npm1), each with k-mers shared with tchem's housekeeping_human.fa dropped. Other species:
undetermined. Cell Ranger: human runs only (one reference on disk).
