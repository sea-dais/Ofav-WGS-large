####### Trim Reads ############
conda activate cutadapt

mkdir -p trimmed
>trim
for file in files/*R1_001.fastq.gz; do
base=$(basename "$file")
echo "cutadapt -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed/${base/R1_001.fastq.gz/R2.trim.gz} $file ${file/R1_001.fastq.gz/R2_001.fastq.gz}" >> trim; done

ls6_launcher_creator.py -j trim -n trim -t 24:00:00 -e dmflores@utexas.edu -q normal -A IBN21018
sbatch trim.slurm

####### Prepare Genome ########
# *skip if already indexed* 
idev -t 00:30:00 
# bowtie index genome 
export GENOME_FASTA=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna
bowtie2-build $GENOME_FASTA $SCRATCH/OfavGenome/ofav_index

####### Map ###################
export GENOME_INDEX=$SCRATCH/OfavGenome/ofav_index

mkdir mapped
>map
for file in trimmed/*R1.trim.gz; do
base=$(basename "$file")
sample=${base/_S*/}
echo "bowtie2 -x $GENOME_INDEX -1 $file -2 ${file/R1.trim.gz/R2.trim.gz} -S mapped/${sample}.sam --un-conc-gz mapped/${sample}_unaligned.fastq.gz -p 16" >> map; done

ls6_launcher_creator.py -j map -n map -a IBN21018 -e dmflores@utexas.edu -t 02:00:00 -q development
sbatch map.slurm

# all of this is to assign read numbers and alignment rates to each samples
grep 'overall alignment rate' map.e3374949 > align_rates
grep 'reads; of these:' map.e3374949 > read_num


conda create -n seqkit -c bioconda -c conda-forge seqkit
conda activate seqkit

seqkit stats trimmed/*.trim.gz > trimreadcounts
seqkit stats files/*.fastq.gz trimmed/*.trim.gz > readcounts_compare



##################
##################
##################
idev -t 00:30:00 -p development
conda activate fastp          # or module load fastp

mkdir -p trimmed2
fastp \
  -i files/e11-SP_S26_L001_R1_001.fastq.gz \
  -I files/e11-SP_S26_L001_R2_001.fastq.gz \
  -o trimmed2/e11-SP_S26_L001_R1.trim.gz \
  -O trimmed2/e11-SP_S26_L001_R2.trim.gz \
  --detect_adapter_for_pe --trim_poly_g --length_required 100 -q 20 \
  -j trimmed2/e11-SP.json -h trimmed2/e11-SP.html

zcat trimmed2/e11-SP_S26_L001_R2.trim.gz | head -20
# should NOT see CTGTCTCTTATACACATCT or long GGGGGG runs anymore

# Re try with cutadapt
cutadapt \
  -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG \
  -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT \
  -q 20 -m 100 --nextseq-trim=20 \
  -o trimmed3/e11-SP_S26_L001_R1.trim.gz \
  -p trimmed3/e11-SP_S26_L001_R2.trim.gz \
  files/e11-SP_S26_L001_R1_001.fastq.gz files/e11-SP_S26_L001_R2_001.fastq.gz

zcat trimmed3/e11-SP_S26_L001_R2.trim.gz | head -20
# want: NO AGATGTGTATAAGAGACAG at read starts, NO CTGTCTCTTATACACATCT at ends


cutadapt \
  -g file:forward.fasta -G file:forward.fasta \
  -a file:reverse.fasta -A file:reverse.fasta \
  -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
  -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 \
  -o trimmed4/e11-SP_S26_L001_R1.trim.gz \
  -p trimmed4/e11-SP_S26_L001_R2.trim.gz \
  files/e11-SP_S26_L001_R1_001.fastq.gz files/e11-SP_S26_L001_R2_001.fastq.gz

zcat trimmed4/e11-SP_S26_L001_R2.trim.gz | head -20
# want clean ends: no CTGTCTCTTATACACATCT, no AGATCGGAAGAGC

bowtie2 -x $SCRATCH/OfavGenome/ofav_index \
  -1 trimmed4/e11-SP_S26_L001_R1.trim.gz \
  -2 trimmed4/e11-SP_S26_L001_R2.trim.gz \
  -S mapped/e11-SP_v4.sam -p 16

samtools flagstat mapped/e11-SP_v4.sam