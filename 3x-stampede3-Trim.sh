> cutadapt.cmds
for file in $(ls files/*R1_001.fastq.gz); do
  base=$(basename "$file")
  echo "cutadapt -j 4 -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed/${base/R1_001.fastq.gz/R2.trim.gz} $file ${file/R1_001.fastq.gz/R2_001.fastq.gz}" >> cutadapt.cmds
done

cat > mylauncher.py << 'EOF'
import pylauncher
pylauncher.ClassicLauncher("cutadapt.cmds", cores=4)
EOF


cat > cutadapt.slurm << 'EOF'
#!/bin/bash
#SBATCH -J cutadapt
#SBATCH -o cutadapt.%j.out
#SBATCH -e cutadapt.%j.err
#SBATCH -p spr
#SBATCH -N 2
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -A IBN21018
#SBATCH --mail-user=dmflores@utexas.edu
#SBATCH --mail-type=all

module load pylauncher
# activate your cutadapt environment here, e.g.:
# source activate cutadapt-env

mkdir -p trimmed
python3 mylauncher.py
EOF

sbatch cutadapt.slurm
import pylauncher
pylauncher.ClassicLauncher("cutadapt.cmds", cores=4)

# This worked so created helper scripts/ wrappers for future jobs x-setUp.sh