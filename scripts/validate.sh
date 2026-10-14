#!/usr/bin/env zsh
# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

# Runs the detector over samples with known layouts and damaged copies, and checks each verdict.
# Inputs live under $TUSSL_CHEMISTRY_DATA (the benchmark data root); whitelists default to its whitelists/.
set -uo pipefail
B=${CARGO_TARGET_DIR:-target}/release/tchem
R=${TUSSL_CHEMISTRY_DATA:?set TUSSL_CHEMISTRY_DATA to the benchmark data root}
export TUSSL_CHEMISTRY_WHITELISTS=${TUSSL_CHEMISTRY_WHITELISTS:-$R/whitelists}
S=$R/samples; M=$S/preflight_mangled; P=$R/pbmc_1k_v3_fastqs
G=$R/giab-chr20-35x
fail=0
# check <name> <wanted substrings separated by '|'> <detector args...>
check() {
  local name=$1 wants=$2; shift 2
  local got; got=$("$B" "$@" 2>/dev/null | tr '\n' ' ')
  local want; local -a parts=("${(@s:|:)wants}")
  for want in "${parts[@]}"; do
    if [[ "$got" != *"$want"* ]]; then
      echo "FAIL  $name: wanted '$want', got: $(echo "$got" | grep -o 'VERDICT.*')"; fail=1; return
    fi
  done
  echo "PASS  $name"
}
check "pbmc_1k_v3"          "10x 3' v3/v3.1  (confidence HIGH)" --r1 $P/pbmc_1k_v3_S1_L001_R1_001.fastq.gz --r2 $P/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
check "5k_pbmc_v3"          "10x 3' v3/v3.1  (confidence HIGH)" --r1 $S/5k_pbmc_v3/R1.fastq.gz --r2 $S/5k_pbmc_v3/R2.fastq.gz
check "pbmc_10k_v3"         "10x 3' v3/v3.1  (confidence HIGH)" --r1 $S/pbmc_10k_v3/R1.fastq.gz --r2 $S/pbmc_10k_v3/R2.fastq.gz
check "10k nextgem v3.1"    "10x 3' v3/v3.1  (confidence HIGH)" --r1 $S/10k_PBMC_3p_nextgem_Chromium_Controller/R1.fastq.gz --r2 $S/10k_PBMC_3p_nextgem_Chromium_Controller/R2.fastq.gz
check "pbmc_1k_v2 (3' v2)"  "10x 3' v2  (confidence HIGH)" --r1 $S/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R1_001.fastq.gz --r2 $S/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R2_001.fastq.gz
check "pbmc_1k_v2, no R2"   "AMBIGUOUS between 10x 3' v2 and 10x 5' v1/v2" --r1 $S/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R1_001.fastq.gz
check "sc5p_v2 (5' v2)"     "10x 5' v1/v2  (confidence HIGH)" --r1 $S/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001_R1_001.fastq.gz --r2 $S/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001_R2_001.fastq.gz
check "v3 swapped"          "10x 3' v3/v3.1  (confidence HIGH)" --r1 $P/pbmc_1k_v3_S1_L001_R2_001.fastq.gz --r2 $P/pbmc_1k_v3_S1_L001_R1_001.fastq.gz
check "v3 R1 trimmed to 26" "but: R1 is 26 bases, too short" --r1 $M/v3_trim26_R1.fastq.gz --r2 $M/v3_R2_head.fastq.gz
check "bulk WGS (negative)" "no known layout" --r1 $G/HG002_chr20_R1.fq --r2 $G/HG002_chr20_R2.fq
E=$S/ena_
check "BD Rhapsody V1 beads"    "BD Rhapsody (original V1 beads)  (confidence HIGH)" --r1 $S/bd_v1_smk_head/R1.fastq.gz --r2 $S/bd_v1_smk_head/R2.fastq.gz
check "BD Rhapsody Enhanced"    "BD Rhapsody (Enhanced Beads)  (confidence HIGH)" --r1 $S/bd_eb_wta_head/R1.fastq.gz --r2 $S/bd_eb_wta_head/R2.fastq.gz
check "SPLiT-seq 2018"          "SPLiT-seq (Rosenberg 2018)  (confidence HIGH)" --r1 ${E}SRR6750044/SRR6750044_1.fastq.gz --r2 ${E}SRR6750044/SRR6750044_2.fastq.gz
check "Parse SRR21740391"       "Parse Evercode WT  (confidence HIGH)" --r1 ${E}SRR21740391/SRR21740391_1.fastq.gz --r2 ${E}SRR21740391/SRR21740391_2.fastq.gz
check "Parse SRR23972326"       "Parse Evercode WT  (confidence HIGH)" --r1 ${E}SRR23972326/SRR23972326_1.fastq.gz --r2 ${E}SRR23972326/SRR23972326_2.fastq.gz
check "Parse ERR14667668"       "Parse Evercode WT  (confidence HIGH)" --r1 ${E}ERR14667668/ERR14667668_1.fastq.gz --r2 ${E}ERR14667668/ERR14667668_2.fastq.gz
check "Parse swapped"           "SWAPPED|Parse Evercode WT  (confidence HIGH)" --r1 ${E}SRR21740391/SRR21740391_2.fastq.gz --r2 ${E}SRR21740391/SRR21740391_1.fastq.gz
check "Drop-seq SRR1853178"     "Drop-seq  (confidence MEDIUM)" --r1 ${E}SRR1853178/SRR1853178_1.fastq.gz --r2 ${E}SRR1853178/SRR1853178_2.fastq.gz
W=$R/wild/cache
check "inDrops ERR10146537"     "SWAPPED|inDrops  (confidence HIGH)" --r1 $W/ERR10146537_1.fastq.gz --r2 $W/ERR10146537_2.fastq.gz
check "ScaleBio SRR28437723"    "ScaleBio RNA  (confidence HIGH)" --r1 $W/SRR28437723_1.fastq.gz --r2 $W/SRR28437723_2.fastq.gz
check "Visium (negative)"       "no known layout" --r1 ${E}SRR24710534/SRR24710534_1.fastq.gz --r2 ${E}SRR24710534/SRR24710534_2.fastq.gz

# --- JSON output: one object a pipeline can act on. Each check runs the tool with --format json
# and asserts on the parsed object (python expression over `d`).
json_check() {
  local name=$1 expr=$2; shift 2
  if "$B" --format json "$@" 2>/dev/null | python3 -I -c "import json,sys; d=json.load(sys.stdin); sys.exit(0 if ($expr) else 1)"; then
    echo "PASS  json: $name"; else echo "FAIL  json: $name: $expr"; fail=1; fi
}
V2=$S/tenx_5p_v2_head/sc5p_v2_hs_PBMC_1k_5gex_S1_L001
json_check "10x v3 match" "d['verdict']=='match' and d['layout']['name']==\"10x 3' v3/v3.1\" and d['cdna_strand']=='forward' and d['files']['barcode_read'].endswith('R1_001.fastq.gz') and d['library']['gene_expression']" \
  --r1 $P/pbmc_1k_v3_S1_L001_R1_001.fastq.gz --r2 $P/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
json_check "5' swapped files untangled" "d['swapped'] and d['files']['barcode_read'].endswith('L001_R1_001.fastq.gz') and d['files']['cdna_read'].endswith('L001_R2_001.fastq.gz') and d['cdna_strand']=='reverse'" \
  --r1 ${V2}_R2_001.fastq.gz --r2 ${V2}_R1_001.fastq.gz
json_check "Parse barcode read is R2" "d['layout']['name']=='Parse Evercode WT' and d['files']['barcode_read'].endswith('_2.fastq.gz') and d['files']['cdna_read'].endswith('_1.fastq.gz')" \
  --r1 ${E}SRR21740391/SRR21740391_1.fastq.gz --r2 ${E}SRR21740391/SRR21740391_2.fastq.gz
json_check "ambiguous v2 without R2" "d['verdict']=='ambiguous' and len(d['candidates'])==2 and d['layout'] is None and d['cdna_strand'] is None" \
  --r1 $S/pbmc_1k_v2_head/pbmc_1k_v2_S1_L001_R1_001.fastq.gz
H2=$R/arc_bench/sra
json_check "VDJ library flagged" "d['verdict']=='match' and not d['library']['gene_expression'] and 'VDJ' in d['library']['kind']" \
  --r1 $H2/SRR31724777/SRR31724777_3.fastq.gz --r2 $H2/SRR31724777/SRR31724777_4.fastq.gz
json_check "missing barcode read" "d['verdict']=='no_match' and 'missing' in d['hint'] and d['layout'] is None" \
  --r1 $R/arc_bench/cache/SRR28707459.fastq.gz
# check_error <name> <wanted stderr substring> <detector args...>: a typed refusal, exit 1.
check_error() {
  local name=$1 want=$2; shift 2
  local err rc
  err=$("$B" "$@" 2>&1 >/dev/null); rc=$?
  if (( rc == 1 )) && [[ "$err" == *"$want"* ]] && [[ "$err" != *panicked* ]]; then
    echo "PASS  $name"
  else
    echo "FAIL  $name: rc=$rc, stderr: ${err[1,200]}"; fail=1
  fi
}
# Kit checks, one row each in scripts/validate_kits.tsv: a new kit adds a row, not a block.
while IFS=$'\t' read -r name want r1 r2 why; do
  case $name in '#'*|name) continue ;; esac
  [ -n "$why" ] || { echo "FAIL  validate_kits.tsv row '$name': missing columns"; fail=1; continue; }
  check "$name" "$want" --r1 "$R/$r1" --r2 "$R/$r2"
done < "${0:A:h}/validate_kits.tsv"
# Long reads, straight from BAM: PacBio's Kinnex tutorial (3' v3) and the heads of Oxford
# Nanopore's GEM-X open data (3' v4 on 293T, 5' v3 on Jurkat), fetched by byte range.
L=$R/longread
check "PacBio Kinnex 3' v3 (BAM)" "10x 3' v3/v3.1, long reads  (confidence HIGH)|3prime:v3" --r1 $L/pacbio_tutorial_ccs.bam --reads 2000
check "ONT GEM-X 3' v4 (BAM)"     "10x 3' v4 (GEM-X), long reads  (confidence HIGH)|3prime:v4" --r1 $L/ont_gemx_293t_3p_v4.head.bam --reads 2000
check "ONT GEM-X 5' v3 (BAM)"     "10x 5' v3 (GEM-X), long reads  (confidence HIGH)|5prime:v3" --r1 $L/ont_gemx_jurkat_5p_v3.head.bam --reads 2000
# A long-read run in two files, as ENA serves some: the first 20,000 reads of design run
# SRR12282461 (ONT, 10x 3' v3) cut in two halves. Both go to the long-read detector, and the
# second is set aside (user test, 12 Oct 2026: it went to the short-read detector before).
check "ONT 3' v3, two files" "Second file, long reads too, set aside|10x 3' v3/v3.1, long reads  (confidence HIGH)" $L/SRR12282461_part1.fq.gz $L/SRR12282461_part2.fq.gz
# A lane's files by position, in a jumbled order: the I1 index file is recognised and set aside,
# and R1/R2 are untangled by swap detection.
check "files by position, I1 set aside" "Index reads, set aside: $P/pbmc_1k_v3_S1_L001_I1_001.fastq.gz|10x 3' v3/v3.1  (confidence HIGH)" \
  $P/pbmc_1k_v3_S1_L001_I1_001.fastq.gz $P/pbmc_1k_v3_S1_L001_R2_001.fastq.gz $P/pbmc_1k_v3_S1_L001_R1_001.fastq.gz
EMPTY_DIR=$(mktemp -d); trap 'rm -rf "$EMPTY_DIR"' EXIT
EMPTY=$EMPTY_DIR/empty.fq.gz; printf '' | gzip > $EMPTY
check_error "empty R1 refused" "has no reads" --r1 $EMPTY --r2 $P/pbmc_1k_v3_S1_L001_R2_001.fastq.gz
check_error "empty R2 refused" "has no reads" --r1 $P/pbmc_1k_v3_S1_L001_R1_001.fastq.gz --r2 $EMPTY
check_error "an index file alone refused" "0 read file(s) are left" $P/pbmc_1k_v3_S1_L001_I1_001.fastq.gz
exit $fail
