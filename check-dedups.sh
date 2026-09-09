# one pool, stats only — assumes mapped/ still has the pre-dedup BAM
samtools collate -@ 8 -O mapped/P1-7x11-PL.sorted.bam \
  | samtools fixmate -m -@ 8 - - \
  | samtools sort -@ 8 - \
  | samtools markdup -s -d 2500 -f P1_markdup_stats.txt -@ 8 - /dev/null
cat P1_markdup_stats.txt


mosdepth --no-per-base -t 4 PB-11x7-PL PB-11x7-PL.dedup.bam

mosdepth --no-per-base -t 4 11-AD 11-AD.dedup.bam


# list positions that still have duplicates
bcftools view -H family_phased.vcf.gz | cut -f1,2 | sort | uniq -d > dup.pos
wc -l dup.pos

# show the full rows (CHROM POS REF ALT + genotypes) for the first few
bcftools view -H family_phased.vcf.gz \
  | awk 'NR==FNR{d[$1"\t"$2]=1; next} ($1"\t"$2) in d' dup.pos - \
  | cut -f1-5,10- \
  | head -40
  
