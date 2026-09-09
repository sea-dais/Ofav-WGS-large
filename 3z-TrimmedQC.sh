################################
####### QC OF TRIMMED READS ########
################################

### Simple code to run fastqc for each file
mkdir trimmedQC

>qc
for file in trimmed/*.trim.gz; do
echo "fastqc -o ./trimmedQC $file" >>qc; done

ls6_launcher_creator.py -j qc -n qc \
            -t 02:00:00 -e dmflores@utexas.edu \
            -q gpu-a100-dev -A IBN21018
sbatch qc.slurm

### View html files that are in the output directory
SOURCE='dmflores@ls6.tacc.utexas.edu:/scratch/08717/dmflores/ofav-wgs-large/trimmedQC'
scp "$SOURCE/*\.html" .
########## looking as expected from past sequenicng
########## mainly good, some poly G tail in the R2 sequences
#*#*#* Results from fastqc can be used to decide what to do in the trimming steps

