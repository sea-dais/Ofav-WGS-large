# Switch to Stampede3 
idev -p spr -N 1 -n 1 -t 02:00:00 -A IBN21018

rsync -av dmflores@ls6.tacc.utexas.edu:$SCRATCH/ofav-wgs-large $SCRATCH

ls files/*R1_001.fastq.gz | wc -l   # should print 31
ls forward.fasta reverse.fasta       # both should be there

mkdir -p trimmed_test

# pick the first sample
file=$(ls files/*R1_001.fastq.gz | head -1)
base=$(basename "$file")

time cutadapt -j 0 -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed_test/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed_test/${base/R1_001.fastq.gz/R2.trim.gz} "$file" "${file/R1_001.fastq.gz/R2_001.fastq.gz}"

## Idev loop 
mkdir -p trimmed
for file in files/*R1_001.fastq.gz; do
base=$(basename "$file")
cutadapt -j 0 -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed/${base/R1_001.fastq.gz/R2.trim.gz} "$file" "${file/R1_001.fastq.gz/R2_001.fastq.gz}"
done
###* Note -j flag: 
## -j 1 (default) — one core, slowest.
## -j N — use exactly N cores (e.g. -j 14 for 14).
## -j 0 — the useful shortcut: autodetect and use all available cores on the machine.
################################################
ls files/*R1_001.fastq.gz | tail -n 13 | while read -r file; do
  base=$(basename "$file")
  cutadapt -j 0 \
    -g file:forward.fasta -G file:forward.fasta \
    -a file:reverse.fasta -A file:reverse.fasta \
    -a AGATCGGAAGAGC -A AGATCGGAAGAGC \
    -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 \
    -o "trimmed/${base/R1_001.fastq.gz/R1.trim.gz}" \
    -p "trimmed/${base/R1_001.fastq.gz/R2.trim.gz}" \
    "$file" "${file/R1_001.fastq.gz/R2_001.fastq.gz}"
done

# To submit job
mkdir -p trimmed
>trim
for file in files/*R1_001.fastq.gz; do
base=$(basename "$file")
echo "cutadapt -j 14 -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed/${base/R1_001.fastq.gz/R2.trim.gz} $file ${file/R1_001.fastq.gz/R2_001.fastq.gz}" >> trim; done

launcher_creator.py -j trim -n trim -t 06:00:00 -e dmflores@utexas.edu -q spr -A IBN21018 -w 8 -N 1
sbatch trim.slurm


###############################
# On Stampede3
cdw
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $WORK/miniconda3
source $WORK/miniconda3/etc/profile.d/conda.sh

conda init bash
conda config --set auto_activate_base false

conda create -n cutadapt -c conda-forge -c bioconda cutadapt 


# Transfer scRNA data into corral? 
rsync -av /work/08717/dmflores/ls6/amil_sc_2022 /corral-repl/utexas/IBN21018/daisy/

lfs quota -h -g G-824270 /corral-repl


echo "alias sq='squeue -u dmflores'" >> ~/.bashrc
source ~/.bashrc