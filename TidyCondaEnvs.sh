/usr/local/etc/taccinfo

#conda remove --name BabrahamRNA --all
#conda remove --name GVA2021 --all
#conda remove --name SATURNenv --all
#conda remove --name freebayes --all
#conda remove --name test --all
#conda remove --name ngstools --all
#conda remove --name seqkit --all
#conda remove --name ncbi_datasets --all
#conda remove --name SAMap --all
#conda remove --name OrthoFinderEnv --all
#conda remove --name samapEnv --all

# 
conda env export --name cutadaptenv > cutadaptenv_backup.yml   # safety net
conda env remove --name cutadaptenv

# Empty envs — nothing inside
#conda env remove --name GVA2-21
#conda env remove --name LDBlock
#conda env remove --name trinityenv

#
conda create --name ngs-tools --clone catch
#conda remove --name catch --all

conda activate ngs-tools
conda install -c bioconda bcftools

# clean up cache space 
conda clean --all #removes unused cached packages, tarbells, and index caches

# see which packages are in each env 
for env in $(conda env list | awk 'NR>2 {print $1}' | grep -v '^#'); do
    echo "=== $env ===" >> all_envs.txt
    conda list --name "$env" >> all_envs.txt
    echo "" >> all_envs.txt
done


scp dmflores@ls6.tacc.utexas.edu:/scratch/08717/dmflores/ofav-wgs-large/all_envs.txt .