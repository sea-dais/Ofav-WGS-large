conda create -n WhatsHap
conda install -c bioconda -c conda-forge whatshap
whatshap --version

bcftools query -l 7x11_ind3.vcf.gz
samtools view -H dedup_rg/11-AD.dedup.bam | grep '^@RG'

mkdir -p ped
cat > ped/7x11.ped << 'EOF'
7x11	11-AD	0	0	1	0
7x11	7-AD	0	0	2	0
7x11	2-7x11-LV	11-AD	7-AD	0	0
7x11	3-7x11-LV	11-AD	7-AD	0	0
EOF

idev -p skx-dev -N 1 -n 1 -t 02:00:00 -A IBN21018
conda activate WhatsHap
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

whatshap phase \
  --ped ped/7x11.ped \
  --reference $REF \
  --output 7x11_phased.vcf.gz \
  vcf/7x11_ind3_dedup.vcf.gz \
  dedup_rg/11-AD.dedup.bam \
  dedup_rg/7-AD.dedup.bam \
  dedup_rg/2-7x11-LV.dedup.bam \
  dedup_rg/3-7x11-LV.dedup.bam

whatshap stats --tsv=7x11_phase_stats.tsv 7x11_phased.vcf.gz
whatshap stats --gtf phased.gtf 7x11_phased.vcf.gz

cat > ped/11x7.ped << 'EOF'
11x7    7-AD        0       0       1       0
11x7    11-AD       0       0       2       0
11x7    C-11x7-LV   7-AD    11-AD   0       0
EOF

whatshap phase \
  --ped ped/11x7.ped \
  --reference $REF \
  --output 11x7_phased.vcf.gz \
  vcf/11x7_individuals.vcf.gz \
  dedup_rg/7-AD.dedup.bam \
  dedup_rg/11-AD.dedup.bam \
  dedup_rg/C-11x7-LV.dedup.bam

whatshap stats --tsv=11x7_phase_stats.tsv 11x7_phased.vcf.gz

bcftools view -h 11x7_phased.vcf.gz 
bcftools view -H 11x7_phased.vcf.gz | head

bcftools view -H 11x7_phased.vcf.gz | grep '|'


# what is longest scaffold
sort -k2,2nr GCF_002042975.1_ofav_dov_v1_genomic.fna.fai | head

tail GCF_002042975.1_ofav_dov_v1_genomic.fna.fai

sort -k2,2nr /scratch/08717/dmflores/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna.fai | head | cut -f1,2
# longest scaffold is 4.7Mb 

sort -k2,2nr /scratch/08717/dmflores/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna.fai | awk '{total+=$2; lengths[NR]=$2} END {half=total/2; sum=0; for(i=1;i<=NR;i++){sum+=lengths[i]; if(sum>=half){print "N50:", lengths[i]; print "Total:", total; print "Scaffolds:", NR; break}}}'

# total records in the VCF
bcftools view -H 7x11_phased.vcf.gz | wc -l
## 413046
# sites where at least one sample's GT is phased
bcftools view -H 7x11_phased.vcf.gz \
  | grep -c '[0-9]|[0-9]'
## 17252


bcftools query -l 7x11_phased.vcf.gz
bcftools view -h 7x11_phased.vcf.gz | grep -E "^##FORMAT"

bcftools query -f "%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:PS:DP]\n" 7x11_phased.vcf.gz | head
zcat 7x11_phased.vcf.gz | grep -v "^#" | grep -c "|"
bcftools view -i 'GT[0]="het" | GT[1]="het"' 7x11_phased.vcf.gz \
  | bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT:%PS]\n' | head


# Informative Sites 
# Stage 1: biallelic SNPs, both parents hardened on depth
bcftools view -m2 -M2 -v snps 7x11_phased.vcf.gz \
  | bcftools view -s 11-AD,7-AD \
  | bcftools filter -e 'FMT/DP[0]<10 | FMT/DP[1]<10' \
  -Oz -o parents_hardened.vcf.gz
bcftools index parents_hardened.vcf.gz

# Stage 2a: BED for samtools depth (at least one parent phased-het)
bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS\n' parents_hardened.vcf.gz \
  | awk 'BEGIN{OFS="\t"} {print $1,$2-1,$2}' \
  | sort -k1,1 -k2,2n > informative_sites.bed

# Stage 2b: matching ANGSD sites file, allele fixed to parental REF/ALT
bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' parents_hardened.vcf.gz \
  | sort -k1,1 -k2,2n > angsd_sites.txt
angsd sites index angsd_sites.txt

awk -F'\t' 'BEGIN{OFS="\t"}
  $3!=$4 && $3 ~ /^[ACGT]$/ && $4 ~ /^[ACGT]$/' \
  angsd_sites.txt > angsd_sites.clean.txt

sort -k1,1 -k2,2n angsd_sites.clean.txt > angsd_sites.sorted.txt
angsd sites index angsd_sites.sorted.txt

wc -l informative_sites.bed
#6552

awk -F'\t' 'BEGIN{OFS="\t"} {print $1, $2-1, $2}' \
  angsd_sites.sorted.txt | sort -k1,1 -k2,2n > informative_sites.bed

wc -l informative_sites.bed
wc -l angsd_sites.sorted.txt   # (or whatever you named the cleaned, indexed file)


grep -v "^NC_" angsd_sites.sorted.txt > angsd_sites.nuc.txt
angsd sites index angsd_sites.nuc.txt

# rebuild BED to match exactly
awk -F'\t' 'BEGIN{OFS="\t"}{print $1,$2-1,$2}' angsd_sites.nuc.txt \
  | sort -k1,1 -k2,2n > informative_sites.bed

wc -l angsd_sites.nuc.txt      # expect 6529

cut -f1 angsd_sites.nuc.txt | sort -u > regions.txt
wc -l regions.txt              # number of distinct nuclear scaffolds