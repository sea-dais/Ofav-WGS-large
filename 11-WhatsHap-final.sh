#!/bin/bash
# =============================================================================
# ofav_segdist_pipeline.sh   (v2 — joint call, single phasing, het-hardened)
#
# Segregation-distortion analysis in a reciprocal Orbicella faveolata cross,
# at the level of phased parental haplotype blocks.
#
#   Parents : 11-AD, 7-AD
#   Phasing larvae (single) : 2-7x11-LV, 3-7x11-LV, C-11x7-LV
#       -> all three are offspring of the SAME two parents; cross direction is
#          irrelevant for establishing PARENTAL haplotypes, so we phase once
#          using all three (more transmission events = better parental phase).
#   Pools (bulk larvae, ~100 each):
#       7x11 direction -> P1, P2, P3
#       11x7 direction -> PA, PB, PC
#       (pools are independent offspring of the same parents; NOT the phasing larvae)
#   Reference : GCF_002042975.1_ofav_dov_v1  (RefSeq contig names, NW_/NC_)
#
# KEY FIX vs v1: parental heterozygous calls are hardened on minor-allele READ
#   SUPPORT, not just depth. freebayes calls like AD=13,2 (2 ALT reads) are
#   false hets (paralogs / sequencing error); they pass DP>=10 but do not
#   transmit, producing offspring pool ALT freq ~0 at "het x hom-ref" sites and
#   railing the block fitter. We require het parents to have >= MIN_ALT_READS
#   alt reads AND allele balance >= MIN_BALANCE.
#
# VALIDATION: after step 6, check_site_consistency.py MUST show
#   het x hom-ref -> mean ~0.25 and het x hom-alt -> mean ~0.75.
#   If not, tighten the hardening before trusting anything downstream.
#
# Requires: freebayes, bcftools, whatshap, angsd(>=0.940), samtools, python3(numpy,scipy)
# Run heavy steps (freebayes, whatshap, angsd) via idev/sbatch, not login node.
# =============================================================================
set -euo pipefail

# ---- paths / params (edit) --------------------------------------------------
SCRATCH="${SCRATCH:?}"
REF="$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna"
BAMDIR="dedup_rg"
WORKDIR="$SCRATCH/ofav-wgs-large"
OUTDIR="poolsANGSD"

P1_NAME="11-AD"; P2_NAME="7-AD"          # parents (VCF sample order fixed below)
LARVAE=("2-7x11-LV" "3-7x11-LV" "C-11x7-LV")
POOLS="P1-7x11-PL,P2-7x11-PL,P3-7x11-PL,PA-11x7-PL,PB-11x7-PL,PC-11x7-PL"
# -----------------------------------------------------------------------------

cd "$WORKDIR"
mkdir -p vcf ped bamlists "$OUTDIR"

# =============================================================================
# STEP 0. Joint variant call on all 5 samples (2 parents + 3 phasing larvae).
#   Joint calling guarantees a consistent site set and uniform AD fields across
#   samples (needed for the het-hardening filter). Order the BAMs parents-first.
# =============================================================================
idev -p spr -N 1 -n 1 -t 04:00:00 -A IBN21018
conda activate freebayes

fasta_generate_regions.py "$REF".fai 100000 > regions_fb.txt

# use 50 cores on icx node 
freebayes-parallel regions.txt 50 \ 
  -f "$REF" \
  --genotype-qualities \
  "$BAMDIR/${P1_NAME}.dedup.bam" \
  "$BAMDIR/${P2_NAME}.dedup.bam" \
  "$BAMDIR/2-7x11-LV.dedup.bam" \
  "$BAMDIR/3-7x11-LV.dedup.bam" \
  "$BAMDIR/C-11x7-LV.dedup.bam" \
  > vcf/family5_raw.vcf 
### was 2.2G
conda activate ngs-tools
bcftools norm -m -any -f "$REF" vcf/family5_raw.vcf -Oz -o vcf/family5_raw.vcf.gz
### Lines   total/split/joined/realigned/skipped:   3739526/167437/0/518396/0
bcftools index -t vcf/family5_raw.vcf.gz

# =============================================================================
# STEP 1. Phase ONCE with all three larvae (pedigree + read-backed).
#   Cross direction doesn't matter for autosomal parental phasing, so list all
#   three larvae as offspring of the same parents. Parent sexes left unknown (0)
#   so the reciprocal-cross larva doesn't trip sire/dam sex checks.
# =============================================================================
cat > ped/family.ped << EOF
FAM	${P1_NAME}	0	0	0	0
FAM	${P2_NAME}	0	0	0	0
FAM	2-7x11-LV	${P1_NAME}	${P2_NAME}	0	0
FAM	3-7x11-LV	${P1_NAME}	${P2_NAME}	0	0
FAM	C-11x7-LV	${P1_NAME}	${P2_NAME}	0	0
EOF

conda activate WhatsHap
whatshap phase \
  --ped ped/family.ped \
  --reference "$REF" \
  --output family_phased.vcf.gz \
  vcf/family5_raw.vcf.gz \
  "$BAMDIR/${P1_NAME}.dedup.bam" \
  "$BAMDIR/${P2_NAME}.dedup.bam" \
  "$BAMDIR/2-7x11-LV.dedup.bam" \
  "$BAMDIR/3-7x11-LV.dedup.bam" \
  "$BAMDIR/C-11x7-LV.dedup.bam"

### send as job
REF="$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna"
BAMDIR=dedup_rg
P1_NAME=7-AD
P2_NAME=11-AD

> whatshap.cmds
echo "whatshap phase --ped ped/family.ped --reference $REF --output family_phased.vcf.gz vcf/family5_raw.vcf.gz $BAMDIR/${P1_NAME}.dedup.bam $BAMDIR/${P2_NAME}.dedup.bam $BAMDIR/2-7x11-LV.dedup.bam $BAMDIR/3-7x11-LV.dedup.bam $BAMDIR/C-11x7-LV.dedup.bam" >> whatshap.cmds

mkjob.sh -n whatshap -j whatshap.cmds -c 1 -e WhatsHap
sbatch whatshap.slurm
### check progress. tail logs/whatshap.cmds.3469178/out0
idev -p pvc -N 1 -n 1 -t 03:00:00 -A IBN21018

bcftools index -t family_phased.vcf.gz
whatshap stats --tsv=family_phase_stats.tsv family_phased.vcf.gz

# number of sites in phased family 
echo "$(bcftools view -H family_phased.vcf.gz | wc -l) sites"
## 3983347

# parent heterozygous sites: 
# positions where P1 is het (unique)
bcftools view -s "${P1_NAME},${P2_NAME}" family_phased.vcf.gz \
  | bcftools view -H -m2 -M2 -v snps -i 'GT[0]="het"' | cut -f1,2 | sort -u > p1_het.pos

# positions where P2 is het (unique)
bcftools view -s "${P1_NAME},${P2_NAME}" family_phased.vcf.gz \
  | bcftools view -H -m2 -M2 -v snps -i 'GT[1]="het"' | cut -f1,2 | sort -u > p2_het.pos

# positions where BOTH are het (set intersection)
comm -12 p1_het.pos p2_het.pos | wc -l


# =============================================================================
# STEP 2. Harden parental genotypes.
#   (a) biallelic SNPs; (b) both parents DP>=MIN_PARENT_DP;
#   (c) THE FIX: any parent called het must have >=MIN_ALT_READS minor-allele
#       reads AND allele balance >=MIN_BALANCE. Drop the record otherwise.
#   AD indexing after -s: FMT/AD[0:0]=P1 REF, [0:1]=P1 ALT, [1:0]=P2 REF, [1:1]=P2 ALT.
#   NOTE: test the AD expression on a few records first (see CHECK below) — bcftools
#   AD/expression syntax is version-sensitive.
# =============================================================================
count () { bcftools view -H "$1" | cut -f1,2 | sort -u | wc -l; }

# 1. Split multiallelic rows into a consistent biallelic representation,
#    left-align and normalize against the reference
bcftools norm -m- -f $REF family_phased.vcf.gz -Oz -o family_norm.vcf.gz
bcftools index -t family_norm.vcf.gz
echo "$(count family_norm.vcf.gz) positions" # 3851082

# 2. subset to the two parents, trim now absent Alts
bcftools view -s "${P1_NAME},${P2_NAME}" family_norm.vcf.gz \
  | bcftools view -a -Oz -o s2_parents.vcf.gz
bcftools index -t s2_parents.vcf.gz
echo "parents-only:          $(count s2_parents.vcf.gz) positions" # 3851082

# 3. keep only biallelic SNPs that are variant in the parents 
bcftools view -m2 -M2 -v snps -c 1 s2_parents.vcf.gz -Oz -o s3_biallelic_snps.vcf.gz
bcftools index -t s3_biallelic_snps.vcf.gz
echo "stage3 biallelic parent SNPs: $(count s3_biallelic_snps.vcf.gz) positions" # 1440687

# 4. drop positions that have >1 SNP row: 
bcftools view -H s3_biallelic_snps.vcf.gz | cut -f1,2 | sort | uniq -d > multi.pos
echo "  (multiallelic SNP positions removed: $(wc -l < multi.pos))" # 2596
bcftools view s3_biallelic_snps.vcf.gz \
  | awk -F'\t' 'NR==FNR{bad[$1"\t"$2]=1; next} /^#/{print; next} !(($1"\t"$2) in bad){print}' \
      multi.pos - \
  | bcftools view -Oz -o s4_unique_snps.vcf.gz
bcftools index -t s4_unique_snps.vcf.gz
echo "stage4 unique-position SNPs:  $(count s4_unique_snps.vcf.gz) positions" # 1438091

# 5. hardening with DP floor, no missing GT, allele balance: 
MIN_PARENT_DP=2                          # parental depth floor
MAX_PARENT_DP=25                            # reject depth outliers
MIN_ALT_READS=2                           # het parent must have >= minor-allele reads
MIN_BALANCE=0.25                          # het parent minor/(minor+major) >= 0.25

bcftools view \
  -e "FMT/DP[0]<${MIN_PARENT_DP} | FMT/DP[1]<${MIN_PARENT_DP} \
    | FMT/DP[0]>${MAX_PARENT_DP} | FMT/DP[1]>${MAX_PARENT_DP} \
    | GT[0]=\"mis\" | GT[1]=\"mis\" \
    | (GT[0]=\"het\" & ( FMT/AD[0:0]<${MIN_ALT_READS} \
        | FMT/AD[0:1]<${MIN_ALT_READS} \
        | FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])<${MIN_BALANCE} \
        | FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1])>(1-${MIN_BALANCE}) )) \
    | (GT[1]=\"het\" & ( FMT/AD[1:0]<${MIN_ALT_READS} \
        | FMT/AD[1:1]<${MIN_ALT_READS} \
        | FMT/AD[1:1]/(FMT/AD[1:0]+FMT/AD[1:1])<${MIN_BALANCE} \
        | FMT/AD[1:1]/(FMT/AD[1:0]+FMT/AD[1:1])>(1-${MIN_BALANCE}) ))" \
  s4_unique_snps.vcf.gz -Oz -o parents_hardened.vcf.gz
bcftools index -t parents_hardened.vcf.gz
echo "stage5 hardened final:        $(count parents_hardened.vcf.gz) positions" # 305327



# CHECK (run manually, don't trust the filter blind):
   bcftools query -i 'GT[0]="het"|GT[1]="het"' \
     -f '%CHROM\t%POS[\t%SAMPLE=%GT:%AD]\n' parents_hardened.vcf.gz | head
#   -> every het shown should have healthy ALT read counts, no X,2 survivors.

# =============================================================================
# STEP 3. Informative-sites file (ANGSD) + region list.
#   >=1 parent het, nuclear only (drop NC_ mito), REF/ALT fixed as major/minor.
# =============================================================================
bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' parents_hardened.vcf.gz \
  | grep -v '^NC_XXXXXXXX\.1'  \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n > angsd_sites.nuc.txt

rm -f angsd_sites.nuc.txt.idx angsd_sites.nuc.txt.bin
angsd sites index angsd_sites.nuc.txt
cut -f1,2 angsd_sites.nuc.txt | sort -u > angsd_sites.regions.txt
echo ">> informative sites: $(wc -l < angsd_sites.nuc.txt)"
## >> informative sites: 62170
# =============================================================================
# STEP 4. Parental expected-freq + phase + block table (from HARDENED VCF).
# =============================================================================
bcftools query -T angsd_sites.regions.txt \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' parents_hardened.vcf.gz \
  | python3 build_parental_table.py > parental_table.tsv

echo ">> parental table rows: $(tail -n +2 parental_table.tsv | wc -l)"
## >> parental table rows: 62170
# =============================================================================
# STEP 5. Per-pool ANGSD: frequency (knownEM=ALT freq) + exact allele counts.
#   -doMajorMinor 3 (REF/ALT from sites), -doMaf 1, -dumpCounts 3 (A,C,G,T),
#   no -SNP_pval (keep near-fixed sites = the signal).
# =============================================================================
# Check depths: 
# positions file for samtools (CHROM<TAB>POS)
cut -f1,2 angsd_sites.nuc.txt > sites.pos.txt

# raw depth at those sites (no MAPQ/BAQ filter)
samtools depth -a -b sites.pos.txt dedup_rg/P1-7x11-PL.dedup.bam \
  | awk '{print $3}' | sort -n \
  | awk '{a[NR]=$1} END{
      print "n covered:", NR;
      print "median:  ", a[int(NR/2)];
      print "90th pct:", a[int(NR*0.9)];
      print "max:     ", a[NR];
      s=0; for(i=1;i<=NR;i++)s+=a[i]; print "mean:    ", s/NR}'
#n covered: 62170
#median:   21
#90th pct: 81
#max:      773
#mean:     35.1323

### use pool_depth_report.sh to summarize the other pools 

# Run ANGSD on pools
for BAM in "$BAMDIR"/*-PL.dedup.bam; do
  NAME=$(basename "$BAM" .dedup.bam)
  echo "$BAM" > "bamlists/${NAME}.bamlist"
  angsd -bam "bamlists/${NAME}.bamlist" -sites angsd_sites.nuc.txt \
        -rf regions.txt -ref "$REF" \
        -doMajorMinor 3 -doMaf 1 -GL 1 \
        -doCounts 1 -dumpCounts 3 -doDepth 1 -maxDepth 2000 \
        -minMapQ 30 -minQ 20 -baq 1 -P 4 \
        -out "$OUTDIR/freq_${NAME}"
done

mkdir -p bamlists logs

> angsd.cmds
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
SITES=$SCRATCH/ofav-wgs-large/angsd_sites.nuc.txt
mkdir -p bamlists poolsANGSD
for BAM in dedup_rg/*-PL.dedup.bam; do
  NAME=$(basename "$BAM" .dedup.bam)
  echo "$BAM" > bamlists/${NAME}.bamlist
  echo "angsd -bam bamlists/${NAME}.bamlist -sites $SITES -rf regions.txt -ref $REF -doMajorMinor 3 -doMaf 1 -GL 1 -doCounts 1 -dumpCounts 3 -doDepth 1 -maxDepth 2000 -minMapQ 30 -minQ 20 -baq 1 -P 4 -out poolsANGSD/freq_${NAME}" >> angsd.cmds
done

wc -l angsd.cmds        # must equal 6

mkjob.sh -n angsd -j angsd.cmds -c 4 -e angsd094
sbatch angsd.slurm


# =============================================================================
# STEP 6. VALIDATE the foundation before any distortion inference.
#   Expect het x hom-ref -> ~0.25, het x hom-alt -> ~0.75. If the means are far
#   off, the parental hardening is still letting false hets through — tighten
#   MIN_ALT_READS / MIN_BALANCE and rebuild from step 2.
# =============================================================================
for POOL in P1-7x11-PL PB-11x7-PL; do
  echo "=== consistency: $POOL ==="
  python3 check_site_consistency.py parental_table.tsv "$OUTDIR" "$POOL" "$MIN_POOL_READS"
done

# =============================================================================
# STEP 7. Block-level transmission-ratio fit (stage A) + aggregation (stage B).
#   Only run once STEP 6 validates. See fit_block_t.py; aggregation across the
#   3 replicates per direction and the 7x11-vs-11x7 contrast comes next.
# =============================================================================
python3 fit_block_t.py \
  --parental parental_table.tsv --pooldir "$OUTDIR" \
  --pools "$POOLS" --minreads "$MIN_POOL_READS" --min-snps 2 \
  --out block_t_perpool.tsv


samtools collate -@ 8 -O mapped/P1-7x11-PL.sorted.bam \
  | samtools fixmate -m -@ 8 - - \
  | samtools sort -@ 8 - \
  | samtools markdup -s -d 2500 -f P1_markdup_stats.txt -@ 8 - /dev/null
cat P1_markdup_stats.txt
[bam_sort_core] merging from 6 files and 8 in-memory blocks...
COMMAND: samtools markdup -s -d 2500 -f P1_markdup_stats.txt -@ 8 - /dev/null
READ: 111,244,662
WRITTEN: 111,244,662
EXCLUDED: 27,330,655
EXAMINED: 83,914,007
PAIRED: 78,300,974
SINGLE: 5,613,033
DUPLICATE PAIR: 43,272,898
DUPLICATE SINGLE: 4,620,196
DUPLICATE PAIR OPTICAL: 20,622,588
DUPLICATE SINGLE OPTICAL: 1,306,097
DUPLICATE NON PRIMARY: 0
DUPLICATE NON PRIMARY OPTICAL: 0
DUPLICATE PRIMARY TOTAL: 47,893,094
DUPLICATE TOTAL: 47,893,094
ESTIMATED_LIBRARY_SIZE: 26,296,085

