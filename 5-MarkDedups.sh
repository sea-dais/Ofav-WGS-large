rsync -av dmflores@stampede3.tacc.utexas.edu:$SCRATCH/ofav-wgs-large/mapped $SCRATCH/ofav-wgs-large

# Deduplicate aligned reads: mark PCR/optical duplicates so that
# downstream analyses (variant calling, coverage) treat each fragment as
# an independent observation.

#   collate  -> group reads by name so mates are adjacent (required by fixmate)
#   fixmate -m -> add mate coordinate + ms/MC tags that markdup needs
#   sort     -> return to coordinate order
#   markdup -s -> flag/remove duplicates; -s prints duplicate-rate stats
#   index    -> index final BAM for random-access tools
conda activate samtools
mkdir -p dedup

>dedup_job
for file in mapped/*.sorted.bam; do
sample=$(basename "$file" .sorted.bam)
echo "samtools collate -@ 8 -O $file | samtools fixmate -m -@ 8 - - | samtools sort -@ 8 - | samtools markdup -s -@ 8 - dedup/${sample}.dedup.bam && samtools index dedup/${sample}.dedup.bam" >> dedup_job
done

ls6_launcher_creator.py -j dedup_job -n dedup_job \
            -t 02:00:00 -e dmflores@utexas.edu \
            -q development -A IBN21018
sbatch dedup_job.slurm


##** forgot to add rg tags
mkdir -p dedup_rg
>addrg
for file in dedup/*.dedup.bam; do
  sample=$(basename "$file" .dedup.bam)
  echo "samtools addreplacerg -r ID:${sample} -r SM:${sample} -r PL:ILLUMINA -o dedup_rg/${sample}.dedup.bam $file && samtools index dedup_rg/${sample}.dedup.bam" >> addrg
done

ls6_launcher_creator.py -j addrg -n addrg -t 01:00:00 -e dmflores@utexas.edu -q development -A IBN21018
sbatch addrg.slurm
## Verify 
samtools view -H dedup_rg/11-AD.dedup.bam | grep '^@RG'
# expect: @RG  ID:11-AD  SM:11-AD  PL:ILLUMINA





################## on stampede3 if needed. 
# STEP 1: Build the commands file — one line per sample.
> dedup.cmds
for file in $(ls mapped/*.sorted.bam); do
  sample=$(basename "$file" .sorted.bam)
  echo "samtools collate -@ 8 -O $file | samtools fixmate -m -@ 8 - - | samtools sort -@ 8 - | samtools markdup -s -f dedup/${sample}.dupstats -@ 8 - dedup/${sample}.dedup.bam && samtools index dedup/${sample}.dedup.bam" >> dedup.cmds
done
wc -l dedup.cmds          # sanity check: should equal your sample count
mkdir -p dedup logs       # output dir + logs dir (launch.py writes to logs/)

# STEP 2: Generate the Slurm script.
#   -c 8 matches samtools -@ 8; on spr, 112/8 = 14 concurrent tasks.
mkjob.sh -n dedup -j dedup.cmds -c 8 -e mapping

# STEP 3: Submit.
sbatch dedup.slurm