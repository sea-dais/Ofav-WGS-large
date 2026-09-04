# Check structure
bcftools view -i 'GT[0]="het" | GT[1]="het"' 7x11_phased.vcf.gz \
  | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' \
  | head -20

# from a directory where bcftools + python3 are available, with your VCF accessible
bash run_parental_table.sh 7x11_phased.vcf.gz parental_table.tsv

# how many block-usable parental sites overlap the ANGSD sites?
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."{print $1"\t"$2}' parental_table.tsv \
  | sort > /tmp/block_sites.txt
awk '{print $1"\t"$2}' angsd_sites.nuc.txt | sort > /tmp/angsd_sites.txt
comm -12 /tmp/block_sites.txt /tmp/angsd_sites.txt | wc -l

#5216


# restrict parental table to measured sites, then count SNPs per block
awk -F'\t' 'NR==FNR{m[$1"\t"$2]=1; next} 
  FNR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="." && (($1"\t"$2) in m){print $8":"$9}' \
  /tmp/angsd_sites.txt parental_table.tsv \
  | sort | uniq -c \
  | awk '{n=$1; for(t=1;t<=10;t++) if(n>=t) c[t]++} END{for(t=1;t<=10;t++) printf "  >=%d SNPs: %d blocks\n", t, c[t]+0}'

#  >=1 SNPs: 1052 blocks
#  >=2 SNPs: 767 blocks
#  >=3 SNPs: 586 blocks
#  >=4 SNPs: 471 blocks
#  >=5 SNPs: 383 blocks
#  >=6 SNPs: 312 blocks
#  >=7 SNPs: 260 blocks
#  >=8 SNPs: 217 blocks
#  >=9 SNPs: 183 blocks
#  >=10 SNPs: 156 blocks

zcat poolsANGSD/freq_P1-7x11-PL.pos.gz | head
zcat poolsANGSD/freq_P1-7x11-PL.counts.gz | head
zcat poolsANGSD/freq_P1-7x11-PL.mafs.gz | head