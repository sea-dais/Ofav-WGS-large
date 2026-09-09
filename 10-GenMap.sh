# Check read lengths. 
conda activate ngs-tools
samtools view P3-7x11-PL.dedup.bam | head -10000 | \
  awk '{print length($10)}' | sort | uniq -c | sort -rn | head

# Install GenMap
conda install -c bioconda genmap
genmap --version   # confirm it's there

idev -p skx-dev -N 1 -n 1 -t 02:00:00 -A IBN21018

# Index reference 
REF=$SCRATCH/OfavGenome/GCF_002042975.1_ofav_dov_v1_genomic.fna      # <-- your actual path
IDX=genmap_index

rm -rf $IDX
genmap index -F $REF -I $IDX

# Compute mappability
K=125            # <-- replace with your read length

genmap map -K $K -E 2 -I $IDX -O mappability_k${K}_e2 -bg -T 4

# -bg writes out a bedgraph; -T 4 uses 4 threads
### output is mappability_k125_e2.bedgraph

scp dmflores@stampede3.tacc.utexas.edu:/scratch/08717/dmflores/mappability_k125_e2.bedgraph .

scp dmflores@ls6.tacc.utexas.edu:/scratch/08717/dmflores/ofav-wgs-large/fastqc/\*html .
