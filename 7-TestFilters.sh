conda activate ngs-tools
conda env export --from-history > ngs-tools-env.yml   # --from-history avoids LS6-specific build pins

scp dmflores@ls6.tacc.utexas.edu:$SCRATCH/ngs-tools-env.yml $WORK
conda env create -f $WORK/ngs-tools-env.yml
# Copy vcf
rsync -av dmflores@ls6.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/vcf/ $SCRATCH/ofav-wgs-large/vcf/

# 1. get an skx-dev session (Stampede3)
idev -p skx-dev -t 01:00:00

bgzip vcf/7x11_family.vcf && tabix -p vcf vcf/7x11_family.vcf.gz
bgzip vcf/11x7_family.vcf && tabix -p vcf vcf/11x7_family.vcf.gz

# Depth distribution: 
# per-sample mean/quantile depth — this is what tells you the 10x vs 50x split
bcftools stats -s - vcf/7x11_family.vcf.gz > stats_7x11.txt
grep '^PSC' stats_7x11.txt

# 1. Parent depth distributions — the binding constraint.
#    What does 11-AD actually look like across sites?
bcftools query -s 11-AD -f '[%DP]\n' vcf/7x11_family.vcf.gz \
  | sort -n \
  | awk '{a[NR]=$1} END{
      print "11-AD  n="NR;
      print "  10th:", a[int(NR*0.10)];
      print "  25th:", a[int(NR*0.25)];
      print "  median:", a[int(NR*0.50)];
      print "  75th:", a[int(NR*0.75)];
      print "  90th:", a[int(NR*0.90)]}'



# count sites where both parents have a non-missing genotype
bcftools query -l vcf/7x11_family.vcf.gz | nl -v0
bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' vcf/7x11_family.vcf.gz \
  | grep -vc '^#']

### or by name: 
bcftools view -i 'GT["7-AD"]!="." & GT["11-AD"]!="."' vcf/7x11_family.vcf.gz | grep -vc '^#'
####### 3.46 million sites where both parents are genotyped.

# Require parent depth:
bcftools query -s 7-AD,11-AD -f '[%DP\t]\n' vcf/7x11_family.vcf.gz \
  | awk 'NR%20==0' \
  | awk -F'\t' '{d7=$1;d11=$2; if(d7==".")d7=0; if(d11==".")d11=0;
      for(f=3;f<=8;f++) if(d7>=f&&d11>=f) c[f]++}
      END{for(f=3;f<=8;f++) print "both >= "f"x: "c[f]*20" sites (est)"}'

#both >= 3x: 943960 sites (est)
#both >= 4x: 652160 sites (est)
#both >= 5x: 483840 sites (est)
#both >= 6x: 393360 sites (est)
#both >= 7x: 327540 sites (est)
#both >= 8x: 282300 sites (est)
############# look at GQ: 
# 7-AD and 11-AD GQ distribution, sampled, missing dropped
bcftools query -s 7-AD,11-AD -f '[%GQ\t]\n' vcf/7x11_family.vcf.gz \
  | awk 'NR%20==0' \
  | awk -F'\t' '
      {for(i=1;i<=2;i++){g=$i; if(g!="."&&g!=""){
         lab=(i==1?"7AD":"11AD");
         if(g<5)b[lab"_00-05"]++;
         else if(g<10)b[lab"_05-10"]++;
         else if(g<20)b[lab"_10-20"]++;
         else if(g<30)b[lab"_20-30"]++;
         else b[lab"_30+"]++ }}}
      END{for(k in b) print k, b[k]}' | sort

####### Parent depth vs parent GQ
bcftools query -s 7-AD,11-AD -f '[%DP\t%GQ\t]\n' vcf/7x11_family.vcf.gz \
  | awk 'NR%20==0' \
  | awk -F'\t' '
      # cols: 1=DP7 2=GQ7 3=DP11 4=GQ11
      {dp7=$1;gq7=$2;dp11=$3;gq11=$4;
       # only sites where both genotyped
       if(dp7=="."||dp11=="."||gq7=="."||gq11=="")next;
       dpmin=(dp7<dp11?dp7:dp11);
       gqmin=(gq7<gq11?gq7:gq11);
       dbin=(dpmin<4?"DP<4":dpmin<6?"DP4-5":dpmin<8?"DP6-7":"DP8+");
       gbin=(gqmin<10?"GQ<10":gqmin<20?"GQ10-19":"GQ20+");
       t[dbin"  "gbin]++}
      END{for(k in t) print k, t[k]*20}' | sort
#of the sites at DP 4-5 (which I'd otherwise cut), how many have high GQ (GQ20+) and could be safely kept? and conversely of my DP 6-7 sites, how many are actually low-GQ and should be cut despite passing the depth floor?  

# GQ only where BOTH parents have a non-missing genotype
bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' vcf/7x11_family.vcf.gz \
  | bcftools query -s 7-AD,11-AD -f '%CHROM\t%POS\t[%GT:%GQ\t]\n' \
  | head -20


############# full depth analysis for records
# write the scan as its own script — no nested-quote problems
cat > $SCRATCH/parentdp_scan.sh << 'EOF'
#!/bin/bash
bcftools query -s 7-AD,11-AD -f '[%DP\t]\n' vcf/7x11_family.vcf.gz \
  | awk -F'\t' '
      {d7=$1; d11=$2; if(d7==".")d7=0; if(d11==".")d11=0; total++;
       for(f=3;f<=8;f++){if(d7>=f && d11>=f) c[f]++}}
      END{print "total: " total;
          for(f=3;f<=8;f++) print "both parents >= " f "x: " c[f] " sites"}' \
  > parentdp_result.txt
EOF
chmod +x $SCRATCH/parentdp_scan.sh

# now the .cmds file is trivially clean — one line, one call
> parentdp.cmds
echo "bash $SCRATCH/parentdp_scan.sh" >> parentdp.cmds

$HOME/bin/mkjob.sh -n parentdp -j parentdp.cmds -c 1 -q skx-dev -t 00:30:00 -e ngs-tools

mkdir -p logs
sbatch parentdp.slurm


##################
rsync -av dmflores@stampede3.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/vcf/ $SCRATCH/ofav-wgs-large/vcf/

############ FINAL FILTER ################
bcftools norm -m -any -f $SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna \
    vcf/7x11_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools view -i 'QUAL>=20' \
  | bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6' -Oz \
      -o vcf/7x11_filtered.vcf.gz
tabix -p vcf vcf/7x11_filtered.vcf.gz

# norm -m -any: splits multiallelica into separate biallelic records first 
# -f <ref>: left-alignes and normalizes against the reference 
# view -m2 -M2 -v snps : keep only biallelic SNPs, drop indels
# view -e ... = "mis": drop sites where either parent is missing
# QUAL>=20: site level quality floor, 
# bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6', DP >= 6, 

##########
# only 28996 survive? let's check: 
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

echo -n "1. raw:                        "
bcftools view -H vcf/7x11_family.vcf.gz | wc -l
#1. raw:   11,146,588
echo -n "2. +norm +biallelic SNP:       "
bcftools norm -m -any -f $REF vcf/7x11_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps -H | wc -l
#2. +norm +biallelic SNP: 10,271,315
echo -n "3. +both parents genotyped:    "
bcftools norm -m -any -f $REF vcf/7x11_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' -H | wc -l
##3. +both parents genotyped: 3,374,384
echo -n "4. +QUAL>=20:                  "
bcftools norm -m -any -f $REF vcf/7x11_family.vcf.gz \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools view -i 'QUAL>=20' -H | wc -l
#4. +QUAL>=20:  675,932
echo -n "5. +parent DP>=6 (final):      "
echo "28996"

##########
# QUAL distribution at parent-genotyped sites (sampled for speed)
bcftools norm -m -any -f $SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna vcf/7x11_family.vcf.gz 2>/dev/null \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools query -f '%QUAL\n' \
  | awk 'NR%20==0' \
  | awk '{if($1<1)b["00-01"]++; else if($1<5)b["01-05"]++; else if($1<10)b["05-10"]++;
          else if($1<20)b["10-20"]++; else if($1<50)b["20-50"]++;
          else if($1<100)b["50-100"]++; else b["100+"]++}
         END{for(k in b) print k, b[k]*20}' | sort

#00-01 2335020
#01-05 127000
#05-10 91120
#100+ 495340
#10-20 144780
#20-50 88300
#50-100 92820

# DP>6 without QUAL Filter: 
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
echo -n "DP>=6 only (no QUAL): "
bcftools norm -m -any -f $REF vcf/7x11_family.vcf.gz 2>/dev/null \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6' \
  | grep -vc '^#'

# QUAL After DP>=6 filter 
bcftools norm -m -any -f $REF vcf/7x11_family.vcf.gz 2>/dev/null \
  | bcftools view -m2 -M2 -v snps \
  | bcftools view -e 'FMT/GT[9]="mis" | FMT/GT[10]="mis"' \
  | bcftools filter -e 'FMT/DP[9]<6 | FMT/DP[10]<6' \
  | bcftools query -f '%QUAL\n' \
  | awk '{if($1<1)b["<1"]++; else if($1<20)b["1-20"]++; else b["20+"]++} END{for(k in b) print k,b[k]}' | sort

# how many scaffolds carry markers, and how concentrated?
bcftools view -H vcf/7x11_filtered.vcf.gz | cut -f1 | sort | uniq -c | sort -rn | head -20
echo "scaffolds with >=1 marker:"
bcftools view -H vcf/7x11_filtered.vcf.gz | cut -f1 | sort -u | wc -l