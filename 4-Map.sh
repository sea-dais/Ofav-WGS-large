rsync -av dmflores@ls6.tacc.utexas.edu:$SCRATCH/OfavGenome $SCRATCH

############### Check time to map
ls -lhS trimmed/*R1.trim.gz | head    # largest samples first
# P3-7x11-PL_S12_L001_R1.trim.gz is 3.2G was 4.2G before trimming 
idev -p spr -N 1 -n 1 -t 02:00:00
# on the node:
export GENOME_INDEX=$SCRATCH/OfavGenome/ofav_index

conda activate mapping # see x-mappingEnv.sh 
time bowtie2 -x $GENOME_INDEX \
  -1 trimmed/P3-7x11-PL_S12_L001_R1.trim.gz -2 trimmed/P3-7x11-PL_S12_L001_R2.trim.gz \
  -S mapped/test.sam --un-conc-gz mapped/test_unaligned.fastq.gz -p 8
## check progress: 
ls -lh mapped/test.sam #35G
## compare -p 8 vs -p 16 
# same subset, 8 threads
time bowtie2 -x $GENOME_INDEX -u 2000000 \
  -1 trimmed/P3-7x11-PL_S12_L001_R1.trim.gz -2 trimmed/P3-7x11-PL_S12_L001_R2.trim.gz \
  -p 8 -S /dev/null

# same subset, 16 threads
time bowtie2 -x $GENOME_INDEX -u 2000000 \
  -1 trimmed/P3-7x11-PL_S12_L001_R1.trim.gz -2 trimmed/P3-7x11-PL_S12_L001_R2.trim.gz \
  -p 16 -S /dev/null

####### Map ###################
export GENOME_INDEX=$SCRATCH/OfavGenome/ofav_index

mkdir -p mapped
for file in trimmed/*R1.trim.gz; do
  case "$file" in *Undetermined*) continue ;; esac    # skip the reject bin
  base=$(basename "$file")
  sample=${base/_S*/}
  echo "bowtie2 -x $GENOME_INDEX -1 $file -2 ${file/R1.trim.gz/R2.trim.gz} --un-conc-gz mapped/${sample}_unaligned.fastq.gz -p 16 | samtools sort -@ 2 -m 1G -o mapped/${sample}.sorted.bam -" >> map
done
wc -l map      # should now be 30, not 31

mkjob.sh -n map -j map -c 18 -t 08:00:00 -e mapping
sbatch map.slurm
# Progress: 
ls -lh mapped*
ls -l logs/map.3438310/success* | wc
###########
for b in mapped/*.sorted.bam; do samtools quickcheck "$b" || echo "BAD: $b"; done

samtools quickcheck -v mapped/a7-SP.sorted.bam
samtools view -H mapped/a7-SP.sorted.bam  
samtools flagstat mapped/a7-SP.sorted.bam

samtools view -H mapped/a7-SP.sorted.bam | grep '@SQ' | wc -l


for b in mapped/*.sorted.bam; do
  echo "=== $b ==="
  samtools flagstat "$b" 2>/dev/null | grep -E 'primary mapped|properly paired'
done

samtools stats mapped/a7-SP.sorted.bam 2>/dev/null | grep ^SN | grep -i insert