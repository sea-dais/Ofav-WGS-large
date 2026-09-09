idev -p spr -N 1 -n 1 -t 02:00:00 -A IBN21018


REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
ls -lh "$REF"    # confirm it exists

# build hom x hom sites
awk -F'\t' 'NR>1 && $7=="hom_x_hom_diff"{print $1"\t"$2"\t"$3"\t"$4}' parental_table.tsv \
  | grep -v '^NC_' | sort -k1,1 -k2,2n > homdiff_sites.txt
wc -l homdiff_sites.txt
rm -f homdiff_sites.txt.idx homdiff_sites.txt.bin

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