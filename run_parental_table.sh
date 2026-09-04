#!/bin/bash
# ============================================================
# run_parental_table.sh
# Build the parental expected-frequency + phase + block table
# from the phased VCF, then print sanity checks.
#
# Usage:
#   bash run_parental_table.sh <phased.vcf.gz> [out.tsv]
#
# Requires: bcftools on PATH, python3, and build_parental_table.py
#           in the SAME directory as this script.
# ============================================================
set -euo pipefail

VCF="${1:?usage: run_parental_table.sh <phased.vcf.gz> [out.tsv]}"
OUT="${2:-parental_table.tsv}"

# locate the python builder next to this script
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PY="$HERE/build_parental_table.py"

if [[ ! -f "$PY" ]]; then
  echo "ERROR: build_parental_table.py not found next to this script ($PY)" >&2
  exit 1
fi

echo ">> Extracting genotypes+phase from $VCF ..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' "$VCF" \
  | python3 "$PY" > "$OUT"

echo ">> Wrote $OUT"
echo

echo "=== total rows (incl. header) ==="
wc -l "$OUT"
echo

echo "=== class distribution ==="
tail -n +2 "$OUT" | cut -f7 | sort | uniq -c | sort -rn
echo

echo "=== informative sites (informative==1) ==="
awk -F'\t' 'NR>1 && $12==1' "$OUT" | wc -l

echo "=== block-usable informative sites (informative==1 AND has block_id) ==="
awk -F'\t' 'NR>1 && $12==1 && $9!="."' "$OUT" | wc -l
echo

echo "=== het_x_hom / hom_x_het block-usable sites only ==="
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."' "$OUT" | wc -l

echo "=== distinct het-x-hom blocks ==="
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."{print $8":"$9}' "$OUT" \
  | sort -u | wc -l

echo "=== SNPs-per-block distribution (top 15 biggest blocks) ==="
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."{print $8":"$9}' "$OUT" \
  | sort | uniq -c | sort -rn | head -15
echo

echo "=== how many blocks have >=N SNPs (power tiers) ==="
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."{print $8":"$9}' "$OUT" \
  | sort | uniq -c \
  | awk '{n=$1; for(t=2;t<=10;t++) if(n>=t) c[t]++} END{for(t=2;t<=10;t++) printf "  >=%d SNPs: %d blocks\n", t, c[t]+0}'