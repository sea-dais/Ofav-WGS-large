conda create -n freebayes2 -c bioconda -c conda-forge freebayes vcflib parallel
conda activate freebayes2

# 1. the genome index (.fai) — freebayes-parallel needs it to make regions
ls ${SCRATCH}/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna.fai
# if missing:  samtools faidx $GENOME_FASTA

# 2. the helper scripts that ship with freebayes
which freebayes-parallel
which fasta_generate_regions.py
###### /work/08717/dmflores/ls6/software/envs/freebayes/bin/fasta_generate_regions.py

# 3. GNU parallel (freebayes-parallel depends on it)

##### conda install conda-forge::parallel
which parallel

export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
# make region chunks (100kb each) from the .fai
awk '{for(i=1;i<=$2;i+=100000) print $1":"i"-"(i+100000<$2 ? i+100000 : $2)}' \
  ${GENOME_FASTA}.fai > regions.txt
head regions.txt
wc -l regions.txt      # should be a lot of lines given 1,933 scaffolds


########################
export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

mkdir -p vcf

## Test how much time is needed: 
# grab the first ~20 regions
head -20 regions.txt > regions_test.txt

idev -t 2:00:00
# run interactively (or a short 30-min dev job) and time it
head -50 regions.txt | shuf | head -50 > regions_test.txt   # (or shuf regions.txt | head -50)

time freebayes-parallel regions_test.txt 24 -f $GENOME_FASTA \
  dedup/11-AD.dedup.bam dedup/7-AD.dedup.bam \
  dedup/2-7x11-LV.dedup.bam dedup/3-7x11-LV.dedup.bam \
  dedup/P1-7x11-PL.dedup.bam dedup/P2-7x11-PL.dedup.bam dedup/P3-7x11-PL.dedup.bam \
  dedup/a7-SP.dedup.bam dedup/b7-SP.dedup.bam \
  dedup/e11-SP.dedup.bam dedup/f11-SP.dedup.bam \
  > vcf/test.vcf

# test2
shuf regions.txt | head -100 > regions_test.txt
time freebayes-parallel regions_test.txt 24 -f $GENOME_FASTA \
  dedup/11-AD.dedup.bam dedup/7-AD.dedup.bam \
  dedup/2-7x11-LV.dedup.bam dedup/3-7x11-LV.dedup.bam \
  dedup/P1-7x11-PL.dedup.bam dedup/P2-7x11-PL.dedup.bam dedup/P3-7x11-PL.dedup.bam \
  dedup/a7-SP.dedup.bam dedup/b7-SP.dedup.bam \
  dedup/e11-SP.dedup.bam dedup/f11-SP.dedup.bam \
  > vcf/test2.vcf
  ## 
grep -vc '^#' vcf/test2.vcf     # variant count — if near 0, estimate is unreliable
wc -l regions.txt                # total regions for the extrapolation

bcftools view vcf/test2.vcf | grep -vc '^#'                          # raw count
bcftools view -e 'QUAL<20' vcf/test2.vcf | grep -vc '^#'             # after a basic QUAL filter
bcftools query -f '%QUAL\n' vcf/test2.vcf | sort -n | head           # low end of the distribution

##############
export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

>fbp
echo "freebayes-parallel regions.txt 24 -f $GENOME_FASTA \
dedup_rg/11-AD.dedup.bam dedup_rg/7-AD.dedup.bam \
dedup_rg/2-7x11-LV.dedup.bam dedup_rg/3-7x11-LV.dedup.bam \
dedup_rg/P1-7x11-PL.dedup.bam dedup_rg/P2-7x11-PL.dedup.bam dedup_rg/P3-7x11-PL.dedup.bam \
dedup_rg/a7-SP.dedup.bam dedup_rg/b7-SP.dedup.bam \
dedup_rg/e11-SP.dedup.bam dedup_rg/f11-SP.dedup.bam \
> vcf/7x11_family.vcf" >> fbp

ls6_launcher_creator.py -j fbp -n fbp \
            -t 05:00:00 -e dmflores@utexas.edu \
            -q vm-small -A IBN21018
sbatch fbp.slurm

# Progress: 
ls -lh vcf/7x11_family.vcf #Finished in 31 minutes on vm-small 

# Check 
bcftools view -h vcf/7x11_family.vcf | tail -1     # header line — should show all 12 sample columns
grep -vc '^#' vcf/7x11_family.vcf                   # count of variant sites
#### 11146588

ls -lh vcf/7x11_family.vcf          # size — 8.6G
head -50 vcf/7x11_family.vcf        # what's actually at the top?


############## 11x7 Family 
export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

>fbp2
echo "freebayes-parallel regions.txt 24 -f $GENOME_FASTA \
dedup_rg/11-AD.dedup.bam dedup_rg/7-AD.dedup.bam \
dedup_rg/C-11x7-LV.dedup.bam \
dedup_rg/PA-11x7-PL.dedup.bam dedup_rg/PB-11x7-PL.dedup.bam dedup_rg/PC-11x7-PL.dedup.bam \
dedup_rg/a7-SP.dedup.bam dedup_rg/b7-SP.dedup.bam \
dedup_rg/e11-SP.dedup.bam dedup_rg/f11-SP.dedup.bam \
> vcf/11x7_family.vcf" >> fbp2

ls6_launcher_creator.py -j fbp2 -n fbp2 \
            -t 02:00:00 -e dmflores@utexas.edu \
            -q normal -A IBN21018
sbatch fbp2.slurm

######## on stampede3 
rsync -av dmflores@ls6.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/dedup_rg $SCRATCH

conda activate freebayes2
conda env export --from-history > freebayes2-env.yml   # --from-history avoids LS6-specific build pins

scp dmflores@ls6.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/freebayes2-env.yml $WORK/
conda env create -f $WORK/freebayes2-env.yml

scp dmflores@ls6.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/regions.txt .


######run 
export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna

> fbp2.cmds
echo "freebayes-parallel regions.txt 24 -f $GENOME_FASTA \
dedup_rg/11-AD.dedup.bam dedup_rg/7-AD.dedup.bam \
dedup_rg/C-11x7-LV.dedup.bam \
dedup_rg/PA-11x7-PL.dedup.bam dedup_rg/PB-11x7-PL.dedup.bam dedup_rg/PC-11x7-PL.dedup.bam \
dedup_rg/a7-SP.dedup.bam dedup_rg/b7-SP.dedup.bam \
dedup_rg/e11-SP.dedup.bam dedup_rg/f11-SP.dedup.bam \
> vcf/11x7_family.vcf" >> fbp2.cmds
mkdir -p vcf

mkjob.sh -n fbp2 -j fbp2.cmds -c 24 -t 02:00:00 -q pvc -e freebayes2
sbatch fbp2.slurm

# Progress: 
ls -lh vcf/11x7_family.vcf #Finished in 31 minutes on vm-small 
