#!/usr/bin/env bash
# Copyright 2026 Tussl Ltd
# SPDX-License-Identifier: Apache-2.0

# Downloads every barcode list the layouts use into <dir>, as data/lists.tsv describes them (the
# registry, rendered from src/lists/registry.rs). Each list is checked against its barcode count,
# its barcode length and the sha256 of the file the code was built with; a list that does not
# check is never written, and its fallback source, if it has one ("-" when not), is tried next.
# Lists belong to their vendors (the registry's terms column), so they live outside the repo.
# Some 10x lists (Flex) ship only inside Cell Ranger, which you download from 10x under its
# licence; give its install directory with --cellranger to copy them from there. Without it those
# lists are skipped, and tchem reports the layouts that need them as not checked.
# Usage: scripts/fetch_whitelists.sh <dir> [--cellranger <cellranger install dir>]
set -euo pipefail
usage="usage: fetch_whitelists.sh <dir> [--cellranger <cellranger install dir>]"
dir=${1:?$usage}
cellranger=""
if [ $# -gt 1 ]; then
  [ "$2" = --cellranger ] && [ $# = 3 ] || { echo "$usage" >&2; exit 2; }
  cellranger=$3
  [ -d "$cellranger/lib/python/cellranger/barcodes" ] \
    || { echo "$cellranger: not a Cell Ranger install (no lib/python/cellranger/barcodes)" >&2; exit 2; }
fi
registry="$(dirname "$0")/../data/lists.tsv"
mkdir -p "$dir"
# macOS has shasum, not sha256sum.
if command -v sha256sum > /dev/null; then sha256=sha256sum; else sha256="shasum -a 256"; fi
scratch=$(mktemp -d)
trap 'rm -rf "$scratch"' EXIT

# normalise <recipe> <download> <out>: one barcode per line, gzip decompressed, CRs and empty
# lines dropped; a CSV keeps column 2 below its header and a TSV (no header) its column 2, sorted
# in byte order and deduplicated;
# "slice:S:L" keeps the distinct substrings at bases S to S+L of every line, sorted;
# "quoted" keeps every double-quoted ACGT string in a source file ("quoted:NAME" only those in
# the block that opens with NAME = and ends at the next closing bracket), sorted and deduplicated.
# column2 <separator> <first line>: column 2 of the download from that line on, sorted and
# deduplicated.
column2() {
  tr -d '\r' < "$scratch/plain" | tail -n +"$2" | cut -d"$1" -f2 | grep -v '^$' | LC_ALL=C sort -u
}

normalise() {
  local recipe=$1 download=$2 out=$3
  if [ "$(head -c 2 "$download" | od -An -tx1 | tr -d ' ')" = 1f8b ]; then
    gunzip -c "$download" > "$scratch/plain"
  else
    cp "$download" "$scratch/plain"
  fi
  case $recipe in
    lines) tr -d '\r' < "$scratch/plain" | grep -v '^$' > "$out" || true ;;
    csv-column-2) column2 , 2 > "$out" || true ;;
    tsv-column-2) column2 $'\t' 1 > "$out" || true ;;
    quoted)
      tr -d '\r' < "$scratch/plain" | grep -oE '"[ACGT]+"' | tr -d '"' \
        | LC_ALL=C sort -u > "$out" || true ;;
    quoted:*)
      tr -d '\r' < "$scratch/plain" \
        | awk -v name="${recipe#quoted:}" '$0 ~ "^" name " *= *[[({]" {on=1} on {print} on && /[])}]/ {exit}' \
        | grep -oE '"[ACGT]+"' | tr -d '"' | LC_ALL=C sort -u > "$out" || true ;;
    slice:*)
      local spec=${recipe#slice:}
      tr -d '\r' < "$scratch/plain" | grep -v '^$' \
        | awk -v s="${spec%%:*}" -v l="${spec##*:}" '{print substr($0, s + 1, l)}' \
        | LC_ALL=C sort -u > "$out" || true ;;
    *) echo "registry: unknown recipe $recipe" >&2; exit 1 ;;
  esac
}

# checks <file> <length> <count> <sha256>: prints why the file fails, or nothing.
checks() {
  local file=$1 length=$2 count=$3 sha=$4 n bad got
  n=$(wc -l < "$file" | tr -d ' ')  # macOS's wc pads the count with spaces
  [ "$n" = "$count" ] || { echo "$n barcodes, expected $count"; return; }
  bad=$(awk -v L="$length" 'length($0) != L || $0 !~ /^[ACGT]+$/' "$file" | head -1)
  [ -z "$bad" ] || { echo "a barcode is not $length bases of ACGT: $bad"; return; }
  got=$($sha256 < "$file" | cut -d' ' -f1)
  [ "$got" = "$sha" ] || echo "sha256 $got, expected $sha"
}

failed=0
# row_sound <length> <count> <sha256> <terms>: prints why a registry row is malformed, or nothing.
# A missing column would shift every field after it, so the row is refused before it is used.
row_sound() {
  [[ $1 =~ ^[0-9]+$ ]] || { echo "length '$1' is not a number"; return; }
  [[ $2 =~ ^[0-9]+$ ]] || { echo "count '$2' is not a number"; return; }
  [[ $3 =~ ^[0-9a-f]{64}$ ]] || { echo "sha256 '$3' is not 64 hex digits"; return; }
  [ -n "$4" ] || echo "the terms column is missing"
}

while IFS=$'\t' read -r file length count sha recipe source fallback terms; do
  case $file in '#'*|file) continue ;; esac
  malformed=$(row_sound "$length" "$count" "$sha" "$terms")
  [ -z "$malformed" ] || { echo "registry row $file: $malformed" >&2; exit 1; }
  fetched=""
  case $source in
    cellranger:*)
      if [ -z "$cellranger" ]; then
        echo "$file  skipped: in Cell Ranger only (--cellranger <dir> copies it)"
        continue
      fi
      source="$cellranger/${source#cellranger:}" ;;
  esac
  for url in "$source" "$fallback"; do
    [ "$url" != - ] || continue
    if [ -f "$url" ]; then
      cp "$url" "$scratch/download"
      ok=yes
    elif curl -sSL --fail --retry 2 -o "$scratch/download" "$url"; then
      ok=yes
    else
      ok=""
    fi
    if [ -n "$ok" ]; then
      normalise "$recipe" "$scratch/download" "$scratch/list"
      why=$(checks "$scratch/list" "$length" "$count" "$sha")
    else
      why="download failed"
    fi
    if [ -z "$why" ]; then
      mv "$scratch/list" "$dir/$file"
      echo "$file  $count barcodes  ($terms)"
      fetched=yes
      break
    fi
    echo "$file from $url: $why" >&2
  done
  if [ -z "$fetched" ]; then
    echo "$file: no source gave the registered list" >&2
    failed=1
  fi
done < "$registry"
[ "$failed" = 0 ] || exit 1

# Optional: packed copies load in milliseconds instead of parsing 315 MB of text.
echo "faster loads: tchem --whitelists $dir --pack"
