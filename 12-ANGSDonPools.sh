# ANGSD
conda create -y -n angsd094 -c bioconda -c conda-forge angsd=0.940
conda activate angsd094
angsd --version 2>&1 | head -1     # confirm it's 0.940, not 0.935 again
#-------------
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

## for testing. 
angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -rf regions.txt -ref $REF \
  -doMajorMinor 3 -doMaf 1 -GL 1 \
  -doDepth 1 -doCounts 1 -maxDepth 2000 \
  -minMapQ 30 -minQ 20 -baq 1 -P 4 \
  -out poolsANGSD/freq_P1-7x11-PL
### test 2
angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -r NW_018148503.1: -ref $REF \
  -doMajorMinor 3 -doMaf 1 -GL 1 \
  -minMapQ 30 -minQ 20 \
  -only_proper_pairs 0 -remove_bads 0 \
  -out test_relaxed
zcat test_relaxed.mafs.gz | tail -n +2 | wc -l
# test 3
angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -r NW_018148503.1: -ref "$REF" \
  -doMajorMinor 3 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 \
  -only_proper_pairs 0 -remove_bads 0 \
  -doCounts 1 -dumpCounts 2 \
  -out test_nofilter
zcat test_nofilter.mafs.gz | tail -n +2 | wc -l
zcat test_nofilter.pos.gz | head        # shows per-site depth
## check the bams
# Are there ANY reads on this contig in the BAM?
samtools view -c dedup_rg/P1-7x11-PL.dedup.bam NW_018148503.1

# Is the BAM even indexed and coordinate-sorted? (ANGSD needs both)
samtools view -H dedup_rg/P1-7x11-PL.dedup.bam | grep '^@HD'
ls -la dedup_rg/P1-7x11-PL.dedup.bam.bai

# Where are the sites on this contig?
awk '$1=="NW_018148503.1"{print $2}' angsd_sites.nuc.txt | head
awk '$1=="NW_018148503.1"{print $2}' angsd_sites.nuc.txt | wc -l

# Where do reads actually cover? Depth at those exact site positions:
awk 'BEGIN{OFS="\t"} $1=="NW_018148503.1"{print $1,$2}' angsd_sites.nuc.txt > /tmp/p1sites.txt
samtools depth -a -b <(awk 'BEGIN{OFS="\t"}$1=="NW_018148503.1"{print $1,$2-1,$2}' angsd_sites.nuc.txt) \
  dedup_rg/P1-7x11-PL.dedup.bam | awk '$3>0' | head

ls -la angsd_sites.nuc.txt*     # check .idx/.bin timestamps vs .txt
rm -f angsd_sites.nuc.txt.idx angsd_sites.nuc.txt.bin
angsd sites index angsd_sites.nuc.txt
# test again: 
angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -r NW_018148503.1: -ref "$REF" \
  -doMajorMinor 1 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_dmm1
zcat test_dmm1.mafs.gz | tail -n +2 | wc -l

# test with no sites
angsd -bam bamlists/P1-7x11-PL.bamlist -r NW_018148503.1: -ref "$REF" \
  -doDepth 1 -doCounts 1 -maxDepth 2000 -minMapQ 0 -minQ 0 \
  -only_proper_pairs 0 -remove_bads 0 -out test_nosites
cat test_nosites.depthGlobal

# Need to reorder sites...
# Get the contig order from the reference .fai
cut -f1 "$REF".fai > ref_order.txt
head ref_order.txt

# Check: does your sites file order match?
cut -f1 angsd_sites.nuc.txt | uniq | head
# Check delimiter? 
head angsd_sites.nuc.txt | cat -A | head    # -A shows tabs (^I) and line-ends ($)
awk -F'\t' '{print NF}' angsd_sites.nuc.txt | sort | uniq -c   # column count distribution

# test sites only: 
angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -ref "$REF" \
  -doMajorMinor 3 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_sitesonly
zcat test_sitesonly.mafs.gz | tail -n +2 | wc -l

angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
  -ref "$REF" -doMajorMinor 3 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_freshidx
zcat test_freshidx.mafs.gz | tail -n +2 | wc -l

ls -la "$REF".fai
cut -f1 "$REF".fai | head
cut -f1 "$REF".fai | wc -l

samtools view -H dedup_rg/P1-7x11-PL.dedup.bam | awk '/^@SQ/{sub(/SN:/,"",$2);print $2}' | sort > /tmp/bam_contigs.txt
cut -f1 "$REF".fai | sort > /tmp/ref_contigs.txt
wc -l /tmp/bam_contigs.txt /tmp/ref_contigs.txt
comm -3 /tmp/bam_contigs.txt /tmp/ref_contigs.txt | head
# test 2col
cut -f1,2 angsd_sites.nuc.txt > sites_2col.txt
angsd sites index sites_2col.txt
angsd -bam bamlists/P1-7x11-PL.bamlist -sites sites_2col.txt \
  -r NW_018148504.1: -ref "$REF" -doMajorMinor 1 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_2col
zcat test_2col.mafs.gz | tail -n +2 | wc -l

samtools depth -r NW_018148504.1:9104-9104 dedup_rg/P1-7x11-PL.dedup.bam
printf 'NW_018148504.1\t9104\n' > one_site.txt
angsd sites index one_site.txt
angsd -bam bamlists/P1-7x11-PL.bamlist -sites one_site.txt \
  -r NW_018148504.1:9104-9104 -ref "$REF" \
  -doMajorMinor 1 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_one
zcat test_one.mafs.gz | tail -n +2 | wc -l


printf 'NW_018148504.1\t9104\n' > one_site.txt
angsd sites index one_site.txt 2>&1 | tee index.log
angsd -bam bamlists/P1-7x11-PL.bamlist -sites one_site.txt \
  -r NW_018148504.1:9104-9104 -ref "$REF" \
  -doMajorMinor 1 -doMaf 1 -GL 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_one 2>&1 | tee run.log

angsd -bam bamlists/P1-7x11-PL.bamlist -sites one_site.txt \
  -r NW_018148504.1:9104-9104 -ref "$REF" \
  -doCounts 1 -dumpCounts 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_sc
zcat test_sc.pos.gz

angsd -bam bamlists/P1-7x11-PL.bamlist \
  -r NW_018148504.1:9100-9110 -ref "$REF" \
  -doCounts 1 -dumpCounts 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_ns
zcat test_ns.pos.gz

angsd -bam bamlists/P1-7x11-PL.bamlist -sites one_site.txt \
  -ref "$REF" -doCounts 1 -dumpCounts 1 \
  -minMapQ 0 -minQ 0 -only_proper_pairs 0 -remove_bads 0 \
  -out test_norflag
zcat test_norflag.pos.gz



##############
zcat freq_P2-7x11-PL.mafs.gz | head
zcat freq_P2-7x11-PL.mafs.gz | tail -n +2 | wc -l

for f in freq_*.mafs.gz; do echo -n "$f: "; zcat "$f" | tail -n +2 | wc -l; done

for f in freq_*.depthGlobal; do
  awk -v n="$f" '{t=0;w=0;for(i=1;i<=NF;i++){t+=$i;w+=(i-1)*$i}; printf "%-32s mean=%.1f n=%d\n",n,(t?w/t:0),t}' "$f"
done

### DEPTH
for f in poolsANGSD/freq_*.depthGlobal; do
  awk -v n="$f" '{t=0;w=0;for(i=1;i<=NF;i++){t+=$i;w+=(i-1)*$i}; printf "%-32s mean=%.1f n=%d\n",n,(t?w/t:0),t}' "$f"
done
#poolsANGSD/freq_P1-7x11-PL.depthGlobal mean=23.1 n=6163
#poolsANGSD/freq_P2-7x11-PL.depthGlobal mean=22.6 n=6101
#poolsANGSD/freq_P3-7x11-PL.depthGlobal mean=22.9 n=6141
#poolsANGSD/freq_PA-11x7-PL.depthGlobal mean=23.6 n=6153
#poolsANGSD/freq_PB-11x7-PL.depthGlobal mean=28.2 n=6251
#poolsANGSD/freq_PC-11x7-PL.depthGlobal mean=24.0 n=6174

zcat poolsANGSD/freq_P1-7x11-PL.counts.gz | head
zcat poolsANGSD/freq_P1-7x11-PL.pos.gz | head
# confirm same length (row-aligned):
echo "counts: $(zcat poolsANGSD/freq_P1-7x11-PL.counts.gz | wc -l)  pos: $(zcat poolsANGSD/freq_P1-7x11-PL.pos.gz | wc -l)"