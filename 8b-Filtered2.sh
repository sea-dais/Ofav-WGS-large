# Remove qual filter
idev -p skx-dev -N 1 -n 1 -t 02:00:00 -A IBN21018

############ FINAL FILTER ################
bcftools norm -m -any -f $SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna \
    vcf/7x11_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6' -Oz \
      -o vcf/7x11_filtered3.vcf.gz
tabix -p vcf vcf/7x11_filtered3.vcf.gz

# norm -m -any: splits multiallelica into separate biallelic records first 
# -f <ref>: left-alignes and normalizes against the reference 
# view -m2 -M2 -v snps : keep only biallelic SNPs, drop indels
# view -e ... = "mis": drop sites where either parent is missing
# QUAL>=20: site level quality floor, 
# bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6', DP >= 6, 

bcftools query -l vcf/11x7_family.vcf.gz | nl -v0

# set these two from the nl output above:
P1=8   # index of 7-AD in the 11x7 VCF
P2=9   # index of 11-AD in the 11x7 VCF

REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

bcftools norm -m -any -f $REF vcf/11x7_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e "FMT/GT[$P1]=\"mis\" | FMT/GT[$P2]=\"mis\"" \
  | bcftools filter -e "FMT/DP[$P1]<6 | FMT/DP[$P2]<6" -Oz \
      -o vcf/11x7_filtered3.vcf.gz
tabix -p vcf vcf/11x7_filtered3.vcf.gz
bcftools view -H vcf/11x7_filtered3.vcf.gz | wc -l

# Quick check DP < 6, every row: both columns >=6
bcftools query -s 7-AD,11-AD -f '[%DP\t]\n' vcf/11x7_filtered3.vcf.gz | head -10

# For OneMap
# 7x11
bcftools view -s 11-AD,7-AD,2-7x11-LV,3-7x11-LV vcf/7x11_filtered3.vcf.gz \
  -Oz -o vcf/7x11_ind3.vcf.gz && tabix -p vcf vcf/7x11_ind3.vcf.gz
# 11x7 (note: different larvae — B and C)
bcftools view -s 11-AD,7-AD,C-11x7-LV vcf/11x7_filtered3.vcf.gz \
  -Oz -o vcf/11x7_ind3.vcf.gz && tabix -p vcf vcf/11x7_ind3.vcf.gz



scp dmflores@stampede3.tacc.utexas.edu:/scratch/08717/dmflores/ofav-wgs-large/vcf/\*filtered3.vcf.gz .

bcftools view -H 7x11_filtered2.vcf.gz | wc -l #29,994
bcftools view -H 7x11_filtered3.vcf.gz | wc -l #450,060


bcftools view -H vcf/7x11_ind3.vcf.gz | wc -l #450,060
bcftools view -H vcf/11x7_ind3.vcf.gz | wc -l #432,997