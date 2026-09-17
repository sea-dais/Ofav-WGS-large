idev -p skx-dev -N 1 -n 1 -t 02:00:00 -A IBN21018

REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
ls -lh "$REF"          # confirm it exists
ls -lh parental_table.tsv poolsANGSD/   # confirm your inputs are here

# build hom x hom sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' s4_unique_snps.vcf.gz \
  | awk -F'\t' '($5=="0/0" && $6=="1/1") || ($5=="1/1" && $6=="0/0"){print $1"\t"$2"\t"$3"\t"$4}' \
  | grep -v '^NC_' \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n > homdiff_sites.txt
wc -l homdiff_sites.txt

rm -f homdiff_sites.txt.idx homdiff_sites.txt.bin
angsd sites index homdiff_sites.txt
cut -f1 homdiff_sites.txt | sort -u > regions_homdiff.txt
wc -l regions_homdiff.txt


REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
ls -lh "$REF"

conda activate angsd094
angsd -bam bamlists/P1-7x11-PL.bamlist -sites homdiff_sites.txt \
  -rf regions_homdiff.txt -ref "$REF" \
  -doMajorMinor 3 -doMaf 1 -GL 1 -doCounts 1 -dumpCounts 3 \
  -minMapQ 30 -minQ 20 -baq 1 -P 4 -out homdiff_P1

# Check: 
python3 - << 'EOF'
import gzip
BASE={'A':0,'C':1,'G':2,'T':3}
want={}
for L in open("homdiff_sites.txt"):
    c,p,r,a=L.split()[:4]; want[(c,int(p))]=(r,a)
freqs=[]
with gzip.open("homdiff_P1.pos.gz","rt") as pf, gzip.open("homdiff_P1.counts.gz","rt") as cf:
    pf.readline(); cf.readline()
    for pl,cl in zip(pf,cf):
        pc=pl.split(); key=(pc[0],int(pc[1]))
        if key not in want: continue
        r,a=want[key]; cnt=[int(x) for x in cl.split()]
        n=cnt[BASE[r]]+cnt[BASE[a]]
        if n>=5: freqs.append(cnt[BASE[a]]/n)
freqs.sort()
if freqs:
    n=len(freqs)
    print(f"CALIBRATION hom x hom-diff (must be 0.5): n={n} mean={sum(freqs)/n:.3f} median={freqs[n//2]:.3f}")
    h=[0]*10
    for v in freqs: h[min(int(v*10),9)]+=1
    print("hist:"," ".join(map(str,h)))
else:
    print("no sites depth>=5 -- check the output prefix matches your ANGSD -out")
EOF
######
conda activate ANGSD
angsd sites index homdiff_sites.txt

# need clean 
awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  homdiff_sites.txt > homdiff_sites.clean.txt
mv homdiff_sites.clean.txt homdiff_sites.txt

rm -f homdiff_sites.txt.idx homdiff_sites.txt.bin
angsd sites index homdiff_sites.txt
wc -l homdiff_sites.txt

# Run Angsd with a larval pool 
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

angsd -bam bamlists/PB-11x7-PL.bamlist -sites homdiff_sites.txt \
  -rf regions_homdiff.txt -ref "$REF" \
  -doMajorMinor 3 -doMaf 1 -GL 1 -doCounts 1 -dumpCounts 3 \
  -minMapQ 30 -minQ 20 -baq 1 -P 4 -out homdiff_PB

#

python3 - << 'EOF'
import gzip
BASE={'A':0,'C':1,'G':2,'T':3}
want={}
for L in open("homdiff_sites.txt"):
    c,p,r,a = L.split()[:4]
    want[(c,int(p))] = (r,a)
freqs=[]
with gzip.open("homdiff_PB.pos.gz","rt") as pf, gzip.open("homdiff_PB.counts.gz","rt") as cf:
    pf.readline(); cf.readline()
    for pl,cl in zip(pf,cf):
        pc=pl.split(); key=(pc[0],int(pc[1]))
        if key not in want: continue
        r,a = want[key]
        cnt=[int(x) for x in cl.split()]
        refn,altn = cnt[BASE[r]], cnt[BASE[a]]
        n = refn+altn
        if n>=5: freqs.append(altn/n)
freqs.sort()
if freqs:
    n=len(freqs)
    print(f"hom x hom-diff (FORCED het, must be 0.5): n={n} mean={sum(freqs)/n:.3f} median={freqs[n//2]:.3f}")
    h=[0]*10
    for v in freqs: h[min(int(v*10),9)]+=1
    print("hist [0..1 in 0.1 bins]:", " ".join(map(str,h)))
else:
    print("no hom-diff sites with depth>=5 — check the join / that ANGSD produced output")
EOF



# Result: 
CALIBRATION hom x hom-diff (must be 0.5): n=167176 mean=0.135 median=0.000
hist: 110893 11905 11642 7790 6705 6045 4905 2518 2716 2057

# At a homxhom diff site, the hom-ALT parent should show 100% alt reads in its own bams 
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
for P in 11-AD 7-AD; do
  echo "dedup_rg/${P}.dedup.bam" > bamlists/${P}.bamlist
  angsd -bam bamlists/${P}.bamlist -sites homdiff_sites.txt \
    -rf regions_homdiff.txt -ref "$REF" \
    -doMajorMinor 3 -doMaf 1 -GL 1 -doCounts 1 -dumpCounts 3 \
    -minMapQ 30 -minQ 20 -baq 1 -P 4 -out homdiff_${P}
done

# compute parents ALT fraction at these sites 
python3 - << 'EOF'
import gzip
BASE={'A':0,'C':1,'G':2,'T':3}
want={}
for L in open("homdiff_sites.txt"):
    c,p,r,a=L.split()[:4]; want[(c,int(p))]=(r,a)

for parent in ["11-AD","7-AD"]:
    freqs=[]
    with gzip.open(f"homdiff_{parent}.pos.gz","rt") as pf, gzip.open(f"homdiff_{parent}.counts.gz","rt") as cf:
        pf.readline(); cf.readline()
        for pl,cl in zip(pf,cf):
            pc=pl.split(); key=(pc[0],int(pc[1]))
            if key not in want: continue
            r,a=want[key]; cnt=[int(x) for x in cl.split()]
            n=cnt[BASE[r]]+cnt[BASE[a]]
            if n>=5: freqs.append(cnt[BASE[a]]/n)
    freqs.sort()
    if freqs:
        n=len(freqs)
        print(f"{parent}: n={n} mean_ALT={sum(freqs)/n:.3f} median={freqs[n//2]:.3f}")
        h=[0]*10
        for v in freqs: h[min(int(v*10),9)]+=1
        print(f"   hist: {' '.join(map(str,h))}")
    else:
        print(f"{parent}: no sites depth>=5")
EOF

11-AD: n=2020 mean_ALT=0.132 median=0.000
   hist: 1754 0 0 0 0 0 0 0 0 266
7-AD: n=4018 mean_ALT=0.109 median=0.000
   hist: 3574 0 1 2 2 1 1 0 3 434


python3 - << 'EOF'
import gzip
BASE={'A':0,'C':1,'G':2,'T':3}
want={}
for L in open("homdiff_sites.txt"):
    c,p,r,a=L.split()[:4]; want[(c,int(p))]=(r,a)
# bin by depth
import collections
byfreq_lowdepth=[]; byfreq_highdepth=[]
with gzip.open("homdiff_P1.pos.gz","rt") as pf, gzip.open("homdiff_P1.counts.gz","rt") as cf:
    pf.readline(); cf.readline()
    for pl,cl in zip(pf,cf):
        pc=pl.split(); key=(pc[0],int(pc[1]))
        if key not in want: continue
        r,a=want[key]; cnt=[int(x) for x in cl.split()]
        n=cnt[BASE[r]]+cnt[BASE[a]]
        if n<5: continue
        f=cnt[BASE[a]]/n
        if n>=20: byfreq_highdepth.append(f)
        else: byfreq_lowdepth.append(f)
for name,arr in [("depth 5-19",byfreq_lowdepth),("depth >=20",byfreq_highdepth)]:
    if arr:
        arr.sort(); n=len(arr)
        h=[0]*10
        for v in arr: h[min(int(v*10),9)]+=1
        print(f"{name}: n={n} mean={sum(arr)/n:.3f} median={arr[n//2]:.3f}")
        print(f"   hist: {' '.join(map(str,h))}")
EOF

depth 5-19: n=159276 mean=0.134 median=0.000
   hist: 106269 11170 10906 7210 6320 5718 4679 2371 2622 2011
depth >=20: n=7900 mean=0.157 median=0.000
   hist: 4624 735 736 580 385 327 226 147 94 46

# at hom x hom-diff sites, the single larvae MUST be heterozygous.
# check what fraction actually are, in the family VCF.
bcftools view -T homdiff_sites.txt family_norm.vcf.gz \
  | bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
  | awk -F'\t' '{
      # columns after pos: C-11x7-LV, 3-7x11-LV, 2-7x11-LV, 7-AD, 11-AD  (per your -l order)
      # larvae are $3,$4,$5 ; parents $6,$7
      for(i=3;i<=5;i++){tot[i]++; if($i ~ /0[|\/]1|1[|\/]0/) het[i]++}
    } END {
      split("C-11x7 3-7x11 2-7x11",names," ");
      for(i=3;i<=5;i++) printf "%s: %d/%d het (%.1f%%)\n", names[i-2], het[i], tot[i], 100*het[i]/tot[i]
    }'

# Try with higher depth floor on parents: 
# 5. hardening with DP floor, no missing GT, allele balance: 
MIN_PARENT_DP=10                          # parental depth floor
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
  s4_unique_snps.vcf.gz -Oz -o parents_hardened2.vcf.gz
bcftools index -t parents_hardened2.vcf.gz
count () { bcftools view -H "$1" | cut -f1,2 | sort -u | wc -l; }
echo "stage5 hardened final:        $(count parents_hardened2.vcf.gz) positions" # 305327


bcftools query -i 'GT[0]="het" | GT[1]="het"' \
  -f '%CHROM\t%POS\t%REF\t%ALT\n' parents_hardened2.vcf.gz \
  | grep -v '^NC_XXXXXXXX\.1'  \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n > angsd_sites2.nuc.txt


angsd sites index angsd_sites2.nuc.txt
cut -f1,2 angsd_sites2.nuc.txt | sort -u > angsd_sites2.regions.txt
echo ">> informative sites: $(wc -l < angsd_sites2.nuc.txt)"
## >> informative sites: 1487

bcftools query -T angsd_sites2.regions.txt \
  -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%PS]\n' parents_hardened2.vcf.gz \
  | python3 build_parental_table.py > parental_table2.tsv

echo ">> parental table rows: $(tail -n +2 parental_table2.tsv | wc -l)"
## >> parental table rows: 1487


# hom x hom-diff from the HARDENED vcf (parents now clean homozygotes)
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' parents_hardened2.vcf.gz \
  | awk -F'\t' '($5=="0/0" && $6=="1/1") || ($5=="1/1" && $6=="0/0"){print $1"\t"$2"\t"$3"\t"$4}' \
  | grep -v '^NC_' \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n > homdiff2_sites.txt
wc -l homdiff2_sites.txt

# THE VALIDATION: single larvae must be ~100% het at these sites now
bcftools view -T homdiff2_sites.txt family_norm.vcf.gz \
  | bcftools query -f '%CHROM\t%POS[\t%GT]\n' \
  | awk -F'\t' '{for(i=3;i<=5;i++){tot[i]++; if($i~/0[|\/]1|1[|\/]0/)het[i]++}}
      END{split("C-11x7 3-7x11 2-7x11",n," "); for(i=3;i<=5;i++)printf "%s: %.1f%% het\n",n[i-2],100*het[i]/tot[i]}'

# Something weird is happening 
bcftools view -T homdiff2_sites.txt family_norm.vcf.gz   | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%SAMPLE=%GT:%AD]\n'   | head -30

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' s4_unique_snps.vcf.gz \
  | awk -F'\t' '
    { split($5,p1,","); split($6,p2,",");   # $5=parent0 AD, $6=parent1 AD (check order!)
      d1=p1[1]+p1[2]; d2=p2[1]+p2[2];
      if(d1<10||d2<10) next;
      f1=p1[2]/d1; f2=p2[2]/d2;
      if( (f1<=0.1 && f2>=0.9) || (f1>=0.9 && f2<=0.1) )
        print $1"\t"$2"\t"$3"\t"$4 }' \
  | grep -v '^NC_' \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n | uniq > homdiff_AD.txt
wc -l homdiff_AD.txt


# Check parent depth
n: 2815164
median: 0
25th: 0
75th: 2.5
90th: 20
mean: 6.84716

printf "both parents >=%d: %d sites\n", t, c[t]+0}'
both parents >=6: 58856 sites
both parents >=8: 26908 sites
both parents >=10: 15816 sites
both parents >=12: 10401 sites
both parents >=14: 7230 sites
both parents >=16: 5398 sites


# hom x hom-diff by AD, both parents >=10
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' s4_unique_snps.vcf.gz \
  | awk -F'\t' '{
      split($5,p1,","); split($6,p2,",");
      d1=p1[1]+p1[2]; d2=p2[1]+p2[2];
      if(d1<10 || d2<10) next;
      f1=p1[2]/d1; f2=p2[2]/d2;
      if((f1<=0.1 && f2>=0.9) || (f1>=0.9 && f2<=0.1))
        print $1"\t"$2"\t"$3"\t"$4
    }' \
  | grep -v '^NC_' \
  | awk -F'\t' 'BEGIN{OFS="\t"} $3!=$4 && $3~/^[ACGT]$/ && $4~/^[ACGT]$/' \
  | sort -k1,1 -k2,2n > homdiff_AD10.txt
wc -l homdiff_AD10.txt

# larvae het check BY AD (larvae = fields 3,4,5), require larva depth >=8
bcftools view -T homdiff_AD10.txt family_norm.vcf.gz \
  | bcftools query -f '%CHROM\t%POS[\t%AD]\n' \
  | awk -F'\t' '{
      split("C-11x7 3-7x11 2-7x11",nm," ");
      for(i=3;i<=5;i++){
        split($i,ad,","); d=ad[1]+ad[2];
        if(d<8) continue;
        f=ad[2]/d; tot[i]++;
        if(f>=0.25 && f<=0.75) het[i]++;
      }
    } END{for(i=3;i<=5;i++) printf "%s: %d/%d het by AD (%.1f%%)\n", nm[i-2], het[i], tot[i], 100*het[i]/tot[i]}'

# Export parent depth: 
# per-parent AD-sum depth at all variant sites, two columns: parent1 parent2
bcftools query -f '[%AD]\t\n' s4_unique_snps.vcf.gz \
  | awk -F'\t' 'BEGIN{OFS="\t"; print "p1_11AD","p2_7AD"} {
      split($1,a,","); split($2,b,",");
      print a[1]+a[2], b[1]+b[2]
    }' > parent_depths.tsv
wc -l parent_depths.tsv
head parent_depths.tsv