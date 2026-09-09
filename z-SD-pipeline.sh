#!/bin/bash
# =============================================================================
# ofav_segdist_pipeline.sh
#
# Segregation-distortion analysis in a reciprocal Orbicella faveolata cross.
# Stage: from a phased family VCF + pool BAMs -> per-pool allele frequencies
#        at parental-informative sites, ready for the block-level TRD test.
#
# Design:
#   Parents:   11-AD, 7-AD
#   Pools:     7x11 direction -> P1, P2, P3   (larval pools, ~100 larvae each)
#              11x7 direction -> PA, PB, PC
#   Reference: GCF_002042975.1_ofav_dov_v1  (RefSeq contig names, NW_/NC_)
#
# Requires: bcftools, angsd (>=0.940), samtools, python3
# Run heavy ANGSD steps via idev/sbatch, not on a login node.
# =============================================================================
set -euo pipefail

# ---- paths (edit these) -----------------------------------------------------
SCRATCH="${SCRATCH:?SCRATCH not set}"
REF="$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna"
PHASED_VCF="$SCRATCH/ofav-wgs-large/7x11_phased.vcf.gz"   # WhatsHap output
BAMDIR="dedup_rg"                                         # holds *-PL.dedup.bam (+ .bai)
WORKDIR="$SCRATCH/ofav-wgs-large"
OUTDIR="poolsANGSD"
P1_NAME="11-AD"                                           # parent sample names in the VCF
P2_NAME="7-AD"
MIN_PARENT_DP=10                                          # harden parental genotype calls
# -----------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p bamlists "$OUTDIR"

# =============================================================================
# STEP 1. Build the informative-sites file from the phased VCF.
#   Informative = biallelic SNP where >=1 parent is heterozygous, both parents
#   confidently genotyped. We fix major/minor to REF/ALT so ANGSD does not
#   re-polarize per site (essential for comparing frequencies across pools).
#   We strip the mitochondrion (NC_*): it is maternally inherited & non-Mendelian.
# =============================================================================

# 1a. Harden parents: biallelic SNPs, both parents above the depth floor.
bcftools view -m2 -M2 -v snps "$PHASED_VCF" \
  | bcftools view -s "$P1_NAME,$P2_NAME" \
  | bcftools filter -e "FMT/DP[0]<$MIN_PARENT_DP | FMT/DP[1]<$MIN_PARENT_DP" \
    -Oz -o parents_hardened.vcf.gz
bcftools index -f parents_hardened.vcf.gz

# 1b. Emit ANGSD sites file (chr pos REF ALT), informative sites only, nuclear only.
#     GT[0]/GT[1] are the two parents in the order given to -s above.
bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' parents_hardened.vcf.gz \
  | grep -v '^NC_' \
  | sort -k1,1 -k2,2n \
  > angsd_sites.nuc.txt

# 1c. Drop any non-ACGT / same-allele rows ANGSD's indexer rejects.
awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  angsd_sites.nuc.txt > angsd_sites.nuc.clean.txt
mv angsd_sites.nuc.clean.txt angsd_sites.nuc.txt

# 1d. Index for ANGSD, and build the region list (contigs that carry sites).
#     NOTE: ANGSD's -sites index is version-keyed. If you change ANGSD version,
#           delete the .idx/.bin and re-run 'angsd sites index'.
rm -f angsd_sites.nuc.txt.idx angsd_sites.nuc.txt.bin
angsd sites index angsd_sites.nuc.txt
cut -f1 angsd_sites.nuc.txt | sort -u > regions.txt

echo ">> informative sites: $(wc -l < angsd_sites.nuc.txt)"

# =============================================================================
# STEP 2. Per-pool allele frequencies + exact per-allele read counts.
#   One ANGSD run per pool, identical settings, shared -sites polarization.
#   -doMajorMinor 3 : take major/minor (REF/ALT) from the sites file.
#   -doMaf 1        : known-major/minor frequency estimator (knownEM = ALT freq).
#   -dumpCounts 3   : per-allele A,C,G,T counts (row-aligned to .pos.gz).
#   No -SNP_pval    : we are NOT discovering SNPs; keep near-fixed sites,
#                     they are the distortion signal.
# =============================================================================
for BAM in "$BAMDIR"/*-PL.dedup.bam; do
  NAME=$(basename "$BAM" .dedup.bam)
  echo "$BAM" > "bamlists/${NAME}.bamlist"
  angsd -bam "bamlists/${NAME}.bamlist" \
        -sites angsd_sites.nuc.txt \
        -rf regions.txt \
        -ref "$REF" \
        -doMajorMinor 3 -doMaf 1 -GL 1 \
        -doCounts 1 -dumpCounts 3 \
        -doDepth 1 -maxDepth 2000 \
        -minMapQ 30 -minQ 20 -baq 1 \
        -P 4 \
        -out "$OUTDIR/freq_${NAME}"
  echo ">> done: $OUTDIR/freq_${NAME}"
done

# =============================================================================
# STEP 3. Parental expected-frequency + phase + block table.
#   build_parental_table.py reads bcftools GT+PS and emits, per site:
#   class, segregating parent, block (PS) id, hapA allele, expected ALT freq.
# =============================================================================
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' "$PHASED_VCF" \
  | python3 build_parental_table.py > parental_table.tsv

echo ">> parental table: $(tail -n +2 parental_table.tsv | wc -l) rows"

# =============================================================================
# NEXT (separate script): join parental_table.tsv to the six freq_*.{mafs,pos,counts}.gz,
# apply a per-site depth floor, fit one transmission ratio t per phased het x hom
# block (beta-binomial on ALT/REF read counts), test t != 0.5, then:
#   - require consistency across the 3 replicates within each cross direction
#   - contrast 7x11 vs 11x7:  same-both-directions = autosomal segregation distortion;
#                             flips-with-direction  = parent-of-origin effect.
# QC: hom x hom-diff sites (expected 0.5, fixed) as a mapping-bias / error check;
#     within-block residual runs as a WhatsHap switch-error flag.
# =============================================================================
echo ">> pool-frequency stage complete."