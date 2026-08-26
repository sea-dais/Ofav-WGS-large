################################
####### QC OF RAW READS ########
################################

### Simple code to run fastqc for each file
mkdir fastqc

>quality
for file in files/*.fastq.gz; do
echo "fastqc -o ./fastqc $file" >>quality; done

ls6_launcher_creator.py -j quality -n quality \
            -t 02:00:00 -e dmflores@utexas.edu \
            -q development -A IBN21018
sbatch quality.slurm

### View html files that are in the output directory
SOURCE='dmflores@ls6.tacc.utexas.edu:/scratch/08717/dmflores/ofav-wgs/fastqc'
scp "$SOURCE/*\.html" .
########## looking as expected from past sequenicng
########## mainly good, some poly G tail in the R2 sequences
#*#*#* Results from fastqc can be used to decide what to do in the trimming steps

