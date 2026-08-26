>download
echo "bs -v download project --name JA26187 \
--extension fastq.gz -o $SCRATCH/ofav-wgs-large/files" > download

ls6_launcher_creator.py -j download -n download \
                        -t 02:00:00 -e dmflores@utexas.edu \
                        -q development -A IBN21018
sbatch download.slurm

cd $SCRATCH/ofav-wgs-large
find files -mindepth 2 -type f -name "*.fastq.gz" -exec mv {} files/ \;
find files -mindepth 1 -type d -empty -delete
