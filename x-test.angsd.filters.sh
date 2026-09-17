

for b in 0 1 2; do
  angsd -bam bamlists/P1-7x11-PL.bamlist -sites angsd_sites.nuc.txt \
        -rf regions.txt -ref "$REF" -doMajorMinor 3 -doMaf 1 -GL 1 \
        -doCounts 1 -dumpCounts 3 -baq $b -minMapQ 30 -minQ 20 -P 4 \
        -out poolsANGSD/TESTbaq${b}
  echo "baq=$b -> sites: $(zcat poolsANGSD/TESTbaq${b}.mafs.gz | tail -n +2 | wc -l)"
done


# for a given pool floor, how many phased blocks still have >=2 informative sites?
# join informative sites' PS blocks against pool depth, count sites/block per threshold
# (rough version using the parents' PS + this pool's depth)

# 1. site -> PS block map (phased het sites only)
bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS[\t%PS]\n' parents_hardened.vcf.gz \
  | awk -F'\t' '{ps=($3!="."?$3:($4!="."?$4:"")); if(ps!="") print $1"_"$2"\t"ps}' \
  | sort > site2block.txt

# 2. pool depth per site
samtools depth -a -b sites.pos.txt -q 20 -Q 30 dedup_rg/P1-7x11-PL.dedup.bam \
  | awk '{print $1"_"$2"\t"$3}' | sort > site2depth.txt

# 3. join, then for each threshold count blocks with >=2 passing sites
join site2block.txt site2depth.txt \
  | awk '{blk=$2; dp=$3;
      for(t=6;t<=10;t+=2) if(dp>=t) cnt[t"_"blk]++}
     END{
      for(t=6;t<=10;t+=2){b=0; for(k in cnt) if(k ~ "^"t"_" && cnt[k]>=2) b++;
        print "floor="t"  blocks with >=2 sites: "b}}'


samtools depth -a -b sites.pos.txt dedup_rg/P2-7x11-PL.dedup.bam \
  | awk '{print $3}' | sort -n \
  | awk '{a[NR]=$1} END{
      print "n covered:", NR;
      print "median:  ", a[int(NR/2)];
      print "90th pct:", a[int(NR*0.9)];
      print "max:     ", a[NR];
      s=0; for(i=1;i<=NR;i++)s+=a[i]; print "mean:    ", s/NR}'


# depth distribution in bins, this pool
samtools depth -a -b sites.pos.txt -q 20 -Q 30 dedup_rg/P2-7x11-PL.dedup.bam \
  | awk '{d=$3;
      if(d==0)b0++;
      else if(d<5)b1++;
      else if(d<10)b2++;
      else if(d<20)b3++;
      else if(d<50)b4++;
      else if(d<100)b5++;
      else b6++}
    END{
      print "0:      "b0;
      print "1-4:    "b1;
      print "5-9:    "b2;
      print "10-19:  "b3;
      print "20-49:  "b4;
      print "50-99:  "b5;
      print ">=100:  "b6}'

## P2-7x11-PL.dedup.bam
n covered: 62170
median:   19
90th pct: 75
max:      711
mean:     32.2423



idev -p spr -N 1 -n 1 -t 08:00:00 -A IBN21018
conda activate ANGSD

REF="$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna"

angsd -bam bamlists/P1-7x11-PL.bamlist \
      -sites angsd_sites.nuc.txt \
      -rf angsd_sites.regions.txt \
      -ref "$REF" \
      -doMajorMinor 3 -doMaf 1 -GL 1 \
      -doCounts 1 -dumpCounts 3 \
      -baq 1 -minMapQ 30 -minQ 20 -P 4 \
      -out poolsANGSD/TESTsites



