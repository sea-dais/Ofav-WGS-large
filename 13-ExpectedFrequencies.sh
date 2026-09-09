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

# Run 
python3 fit_block_t.py \
  --parental parental_table.tsv \
  --pooldir poolsANGSD \
  --pools P1-7x11-PL,P2-7x11-PL,P3-7x11-PL,PA-11x7-PL,PB-11x7-PL,PC-11x7-PL \
  --minreads 5 \
  --min-snps 2 \
  --out block_t_perpool.tsv

# blocks fit per pool
cut -f1 block_t_perpool.tsv | tail -n +2 | sort | uniq -c

# median t_hat (should be ~0.5 if polarization is correct)
tail -n +2 block_t_perpool.tsv | cut -f7 | sort -n | awk '{a[NR]=$1} END{print "median t:", a[int(NR/2)]}'

# strongest candidates by p-value (don't trust single pools yet)
tail -n +2 block_t_perpool.tsv | sort -t$'\t' -k12,12g | head -20

tail -n +2 block_t_perpool.tsv | awk -F'\t' '$2=="7x11"{print $7}' | sort -n | awk '{a[NR]=$1} END{print "7x11 median t:", a[int(NR/2)]}'
tail -n +2 block_t_perpool.tsv | awk -F'\t' '$2=="11x7"{print $7}' | sort -n | awk '{a[NR]=$1} END{print "11x7 median t:", a[int(NR/2)]}'

# pull hom_x_hom_diff sites, join to a pool's counts, look at the ALT-freq distribution
awk -F'\t' 'NR>1 && $7=="hom_x_hom_diff"{print $1"\t"$2"\t"$3"\t"$4}' parental_table.tsv > homdiff_sites.txt
wc -l homdiff_sites.txt

tail -n +2 block_t_perpool.tsv | cut -f7 | \
  awk '{b=int($1*10); h[b]++} END{for(i=0;i<=10;i++) printf "%.1f-%.1f: %d\n", i/10,(i+1)/10,h[i]+0}'

# For het x hom-ref sites (exp 0.25) in one pool, how many have observed ALT freq > 0.5?
# and for het x hom-alt (exp 0.75), how many < 0.5?
awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="."{print $1"\t"$2"\t"$3"\t"$4"\t"$11}' parental_table.tsv > exp_sites.txt



python3 check_sites_consistency.py parental_table.tsv poolsANGSD PB-11x7-PL 5

awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="." && $11=="0.25"{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8}' parental_table.tsv | head -15

awk -F'\t' 'NR>1 && ($7=="het_x_hom"||$7=="hom_x_het") && $9!="." && $11=="0.25"{print $1"_"$2}' parental_table.tsv | head -15 > /tmp/chk_sites.txt

paste <(zcat poolsANGSD/freq_PB-11x7-PL.pos.gz) <(zcat poolsANGSD/freq_PB-11x7-PL.counts.gz) | awk 'NR>1{print $1"_"$2"\t"$1"\t"$2"\tA="$4" C="$5" G="$6" T="$7" dep="$3}' | grep -Ff /tmp/chk_sites.txt | head -15

# hom x hom-diff sites: parents are 0/0 and 1/1, every offspring MUST be het -> ALT freq 0.5
awk -F'\t' 'NR>1 && $7=="hom_x_hom_diff"{print $1"_"$2"\t"$1"\t"$2"\t"$3"\t"$4}' parental_table.tsv | head -2000 > /tmp/homdiff.txt
wc -l /tmp/homdiff.txt

paste <(zcat poolsANGSD/freq_PB-11x7-PL.pos.gz) <(zcat poolsANGSD/freq_PB-11x7-PL.counts.gz) \
  | awk 'NR>1{print $1"_"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7}' \
  | grep -Ff <(cut -f1 /tmp/homdiff.txt) | head -20

# same 15 het x hom-ref sites, but observed in P1 (a 7x11 pool) instead of PB
paste <(zcat poolsANGSD/freq_P1-7x11-PL.pos.gz) <(zcat poolsANGSD/freq_P1-7x11-PL.counts.gz) \
  | awk 'NR>1{print $1"_"$2"\t"$1"\t"$2"\tA="$4" C="$5" G="$6" T="$7" dep="$3}' \
  | grep -Ff /tmp/chk_sites.txt | head -15

# For 7-AD-segregating sites: observed ALT freq in P1 (7x11) vs PB (11x7)
awk -F'\t' 'NR>1 && $8=="7-AD" && $9!="." && $11=="0.25"{print $1"_"$2}' parental_table.tsv | head -500 > /tmp/seg7.txt

echo "=== 7-AD-segregating sites in P1 (7x11) ==="
paste <(zcat poolsANGSD/freq_P1-7x11-PL.pos.gz) <(zcat poolsANGSD/freq_P1-7x11-PL.counts.gz) \
  | awk 'NR>1{print $1"_"$2"\t"$4+$5+$6+$7}' | grep -Ff /tmp/seg7.txt | head -8

echo "=== 7-AD-segregating sites in PB (11x7) ==="
paste <(zcat poolsANGSD/freq_PB-11x7-PL.pos.gz) <(zcat poolsANGSD/freq_PB-11x7-PL.counts.gz) \
  | awk 'NR>1{print $1"_"$2"\t"$4+$5+$6+$7}' | grep -Ff /tmp/seg7.txt | head -8

python3 check_sites_consistency.py parental_table.tsv poolsANGSD P1-7x11-PL 5


# For the parents, at het x hom sites, what's the het parent's allele balance and depth?
# Pull from the phased VCF the AD (allele depth) for the het parent at these sites.
bcftools query -f '%CHROM\t%POS[\t%SAMPLE=%GT:%AD:%DP]\n' \
  -r NW_018148504.1:9104-9146 7x11_phased.vcf.gz