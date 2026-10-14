# BD Rhapsody part numbers and bead types (2026-10-09)

This directory has a reference table of the BD part numbers and the kit and bead names cited in
the 511 held-out BD runs. It says which ones identify a bead type.

Inputs were the audit files `audit/parse_bd/tokens_classified.tsv`, `runs.tsv` and `fields.tsv`
(the last one only to read the sentence around a number), plus public web pages.
- No sequence data was opened.
- `tchem` was not run.
- No bench results file was opened.

## Files
- `parts.tsv`: one row per six-digit part number (20) and per normalised kit or bead name (24).
  - Names are normalised: lower case, with ™, ®, "TM" and "?" removed, and a leading "BD " or "Rhapsody " removed.
  - The `quote` column holds one or more verbatim fragments, joined by ` … `.
  - `scripts/build.py` checks that every fragment is a substring of the saved source text after whitespace normalisation, and that each quote is no longer than 300 characters.
  - When a row has more than one source, the URLs are joined by ` ; `.
- `runs_beads.tsv`: one row per BD run (511).
  - `parts_cited` lists the part numbers and normalised names found in the tokens of classes B_CATALOGUE, B_KIT_NAME, B_BEAD and B_KIT_VERSION.
  - A label is implied only by a row whose class is a bead class.
- `raw/`: the fetched pages and PDFs, with their text extractions. `fetched_urls.tsv` lists the URL of each fetch.
- `scripts/`: `get.sh` (curl), `extract_pdp.py`, `totext.py` and `build.py`.

## Classification rules that were applied
- **633731 (Cartridge Reagent Kit) is `original_v1_beads`.**
  - Its product insert (Doc 54141) lists `650000089 Cell Capture Beads`.
  - BD protocol 210966 Rev 3 sets 650000089 apart from `Enhanced Cell Capture Beads (Cat. No. 700027881)`.
  - The BD flyer BD-52956 says "original BD Rhapsody™ Cell Capture Beads (V1)".
  - The product page alone could not decide this. Its text is word for word the same as the 664887 page.
- **664887 (Enhanced Cartridge Reagent Kit) is `not_bead_specific`, but only within the enhanced family.**
  - BD's 2021–2022 documents (instrument guide 214062 Rev 5, protocols 210966 and 210967 Rev 3) pair 664887 with "Enhanced Cell Capture Beads (Cat. No. 700027881)".
  - The current product insert, 23-24070(02) dated 2023-01, lists `700034960 Enhanced Cell Capture Beads v2`.
  - BD's pipeline guide and its bioinformatics docs treat "Enhanced beads" and "Enhanced V2/V3 beads" as separate layouts.
  - So one catalogue number has carried two bead parts in two layouts. Original V1 beads are excluded, but the number alone cannot tell Enhanced from Enhanced V2.
  - The same class goes to every item that contains or requires 664887: 665915, 665913, 666622 and 665345. It also goes to the names "enhanced cartridge reagent kit" and "enhanced cartridge reagent".
- **"enhanced cell capture beads" and "enhanced beads" (bare names) are `not_bead_specific`.**
  - BD's 2022 protocols use the bare name for 700027881.
  - BD's 2023 HT instrument guide uses the same bare name for the beads of the 664887 workflow, and that kit's insert lists v2.
- **"cartridge reagent kit", "cell capture beads", "beads" and "bead" (bare names) are `not_bead_specific`.**
  - BD's 01/2026 Protocol Planner prints "BD Rhapsody Cartridge Reagent Kit (Cat. No. 664887)".
  - The 664887 product page says "capture them on BD Rhapsody™ Cell Capture Beads".
- **"enhanced cartridge reagent kit v3" is `enhanced_v2v3_beads`.**
  - Insert 23-24738 for 667052 lists `91-1294 BD Rhapsody™ Enhanced Cell Capture Beads v3`.
  - BD's bioinformatics docs group "Enhanced V2/V3 beads".
- **633802 (WTA Reagent Kit, 4 Pack) is `not_bead_specific`.**
  - The 2019 workflow sheet (Doc 220524 Rev 4.0) says it includes Kit 3 = 633731.
  - BD instrument guide 214062 Rev 5.0 (2021-11) lists "Cartridge Reagent Kit or the ... Enhanced Cartridge Reagent Kit" as "included in each of the following 4-pack products", and 633802 is in that list.
  - No BD contents list shuts out the enhanced kit, so 633802 was not called V1.
- **Instruments, cartridges, cDNA, amplification and multiplexing kits, and primer panels are `not_bead_specific`.**
  - This covers 666625 (HT Xpress package), 666262 (8-Lane Cartridge) and 633733 (Cartridge Kit).
  - None of these contains beads, and BD's text ties none of them to one bead.
- **Rows set to `unknown`:**
  - 633663 and 633664: no BD page, document or search hit. The audit sentence calls them "BD Rhapsody cDNA kit" and "BD Target Amplification kit". These were not treated as typos.
  - "enhanced cell capture beads v1" and "standard bd rhapsody beads": not BD names.
  - "v1 bd rhapsody": a protocol-version phrase. Its sentences say the protocol was "applied to enhanced beads".
- No STARsolo adapter or linker strings in the run text were used to assign labels.

## Results
| Bead label implied | Runs | Studies |
|---|---:|---:|
| BD Rhapsody (original V1 beads) | 17 | 9 |
| BD Rhapsody (Enhanced V2/V3 beads) | 2 | 1 |
| none | 492 | 268 |

- No run has a conflict, meaning more than one bead label.
- 307 runs cite no part and no name at all.
- 40 runs (22 studies) cite an item from the enhanced family. That shows the beads were not original V1, but it does not say which enhanced layout.
- **Tension:** PRJNA1183407 (SRR31273785 and SRR31273786) cites 633731, a V1 kit, together with 666622. The 666622 component table lists 664887, an enhanced-family kit. The label comes only from 633731, so the run is not marked as a conflict.

## Sources consulted (all BD-owned)
- **bdbiosciences.com product pages** (en-us, en-in or en-de where en-us gave 404) for 633773, 633801, 633774, 633781, 633793, 633750, 633753, 633731, 633733, 664887, 667052, 665345, 665913, 665915, 666262, 666622, 666625. Also the Accessory Kits page and the Single-Cell Multiplexing Kit page (626545).
- **Product inserts:** 54141 (633731), 54138 (633733), 23-24070(02) (664887), 23-24738(01) (667052), 23-21714(03) (633801), 23-21337(01) (633781), 23-24069(01) (665345), and the 633772 Training Kit insert 219038.
- **Protocols and guides:**
  - 210966 Rev 3 / 23-22951(01) and 210967 Rev 3 / 23-22952(01): cell capture and cDNA synthesis.
  - Instrument guide 214062 Rev 5 / 23-21336(03).
  - HT instrument guide 23-24257(01).
  - WTA protocols 23-24117(01) and 23-24118(01), and the VDJ+WTA protocol 23-24017(01).
  - Sequence Analysis Pipeline guide 23-24580(01).
  - Workflow sheet 220524 Rev 4.0.
  - Older cell capture guides.
- **BD flyers:** BD-52956 (Enhanced beads) and BD-85136 (HT).
- **scomix.bd.com (BD support):** Protocol Planner 01/2026 and Reagent Menu v2 Dec-2025.
- **bd-rhapsody-bioinfo-docs.genomics.bd.com:** the cell-label step page and the release notes.
- No third-party reseller page was used as a source.

## Unreachable or not found
- 633663 and 633664: no product page in en-us, en-in, ja-jp, en-gb or en-de, and no search hit.
- 633802 and 626545 have no product page of their own (404 in all the locales tried). They were resolved from BD guides and from the Single-Cell Multiplexing Kit page.
- The en-us product pages for 665345, 665915 and 666622 return 404. The en-de and en-in pages were used.
- Product inserts 23-22961 and 23-24069 at the guessed DAM path returned 404. The 665345 insert was found at `665345.pdf`.
- No earlier revision of insert 23-24070, which might list 700027881 as a 664887 component, was found online. The Wayback Machine has no copy.
- bdbiosciences.com returns 403 to curl without browser-like headers. With those headers, every fetch worked.
- Product-page component tables (for example 667052) show "N/A" for the component names. Component names were taken from the PDF inserts.

## BD inconsistencies seen (recorded, not resolved)
- The Protocol Planner 01/2026 prints "BD OMICS-One WTA Next Amplification Kit (Cat No. 633801)" and "Targeted mRNA and AbSeq Amplification Kit (Cat. No. 572620)".
  - The Reagent Menu v2 (Dec 2025) gives 572620 as the WTA Next Amplification Kit.
  - The 633801 product page gives 633801 as the WTA Amplification Kit.
  - The Planner was not used as a source for 633801 or 633774.
- The 2019 workflow sheet gives the cDNA kit as 633732. Everywhere else it is 633773.
