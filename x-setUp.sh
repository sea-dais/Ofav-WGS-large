# ============================================================
# STAMPEDE3 pylauncher WORKFLOW (parametric / high-throughput jobs)
# Replaces ls6_launcher_creator.py from Lonestar6.
# ============================================================

# --- ONE-TIME SETUP (do these once, ever) --------------------

# 1. Install paramiko (pylauncher needs it; missing by default).
#    Run on a login node — this is a lightweight exception to
#    the "no compute on login nodes" rule.
module load pylauncher
pip install --user paramiko

# 2. Shared driver — one file for ALL future jobs. Never edit again.
#    Takes the commands file + cores-per-task as arguments.
cat > $HOME/launch.py << 'EOF'
import sys, os, pylauncher
cmds = sys.argv[1]
cores = int(sys.argv[2])
jobid = os.environ.get("SLURM_JOB_ID", "local")
pylauncher.ClassicLauncher(cmds, cores=cores,
                           workdir=f"logs/{cmds}.{jobid}")
EOF

# 3. Slurm-script builder (my version of launcher_creator.py).
#    Flags: -n name  -j cmdfile  -c cores  [-q queue] [-N nodes] [-t time] [-e env]
cat > $HOME/mkjob.sh << 'EOF'
#!/bin/bash
queue=spr; nodes=1; time=04:00:00; env=CHANGE_ME
while getopts n:j:c:q:N:t:e: flag; do
  case $flag in
    n) name=$OPTARG ;;  j) cmds=$OPTARG ;;  c) cores=$OPTARG ;;
    q) queue=$OPTARG ;; N) nodes=$OPTARG ;; t) time=$OPTARG ;;
    e) env=$OPTARG ;;
  esac
done
cat > ${name}.slurm << SLURM
#!/bin/bash
#SBATCH -J ${name}
#SBATCH -o logs/${name}.%j.out
#SBATCH -e logs/${name}.%j.err
#SBATCH -p ${queue}
#SBATCH -N ${nodes}
#SBATCH -n 1
#SBATCH -t ${time}
#SBATCH -A IBN21018
#SBATCH --mail-user=dmflores@utexas.edu
#SBATCH --mail-type=all

module load pylauncher
source \$WORK/miniconda3/etc/profile.d/conda.sh
conda activate ${env}

python3 \$HOME/launch.py ${cmds} ${cores}
SLURM
echo "Wrote ${name}.slurm — sbatch ${name}.slurm"
EOF
chmod +x $HOME/mkjob.sh



# --- PER-JOB WORKFLOW (repeat for each new analysis) ---------
# Work from $SCRATCH. Inputs, outputs, fasta refs all live there.

# STEP 1: Build the commands file — ONE line per sample.
#   The leading `> file.cmds` truncates first so reruns don't stack lines.
#   (This example = cutadapt on paired-end reads; swap in whatever tool.)
> cutadapt.cmds
for file in $(ls files/*R1_001.fastq.gz); do
  base=$(basename "$file")
  echo "cutadapt -j 4 -g file:forward.fasta -G file:forward.fasta -a file:reverse.fasta -A file:reverse.fasta -a AGATCGGAAGAGC -A AGATCGGAAGAGC -n 3 -q 20 -m 100 --nextseq-trim=20 -e 0.2 -o trimmed/${base/R1_001.fastq.gz/R1.trim.gz} -p trimmed/${base/R1_001.fastq.gz/R2.trim.gz} $file ${file/R1_001.fastq.gz/R2_001.fastq.gz}" >> cutadapt.cmds
done
wc -l cutadapt.cmds        # sanity check: should equal your sample count
mkdir -p trimmed           # make the output dir the commands write to

# STEP 2: Generate the Slurm script.
#   -c must match the tool's own thread flag (cutadapt -j 4  ->  -c 4)
#   and divide evenly into node cores (spr=112, skx=48).
$HOME/mkjob.sh -n cutadapt -j cutadapt.cmds -c 4 -e cutadapt-env

# STEP 3: Submit.
sbatch cutadapt.slurm


# --- MONITORING ----------------------------------------------
sq                          # alias for: squeue -u dmflores
tail cutadapt.*.out         # pylauncher stats + cutadapt summaries
tail cutadapt.*.err         # errors if a task failed
ls trimmed/*.trim.gz | wc -l   # confirm all outputs (2 per sample)

# --- NOTES / GOTCHAS -----------------------------------------
# * cutadapt OVERWRITES existing outputs silently. Move old runs
#   aside (mv trimmed trimmed_run1) if you need to keep them.
# * Core math: cores-per-task x concurrent-tasks <= node cores.
#   spr node = 112 cores; -c 4 -> 28 tasks run at once.
# * Only the .cmds file and conda env change between different
#   tools. launch.py and mkjob.sh are permanent.
# * Rerunning the SAME job = just `sbatch cutadapt.slurm`. The
#   3 files persist; nothing to rebuild.
# * pylauncher leaves a pylauncher_tmp_<jobid>/ work dir showing
#   which tasks finished — useful for rerunning only failures.
# ============================================================


mkdir -p $HOME/bin
mv $HOME/mkjob.sh $HOME/launch.py $HOME/bin/
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# make sure to also update on shell script mkjob.sh
###* python3 $HOME/bin/launch.py ....