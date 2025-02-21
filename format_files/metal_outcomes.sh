#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1:00:00
#SBATCH --array=8,10
#SBATCH --job-name=metal
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=5
#SBATCH --mem-per-cpu=3G

module load StdEnv/2020
module load gcc/9.3.0
module load metal/2011-03-25

project_dir="/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# 2,4,6,8,10
SLURM_ARRAY_TASK_ID2=$(( ${SLURM_ARRAY_TASK_ID} + 1 ))

trait_M=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/gwas_sample_sizes_sexstr_neurodegen.tsv | awk '{print $1}')
trait_F=$(sed -n ${SLURM_ARRAY_TASK_ID2}p ${project_dir}/gwas_sample_sizes_sexstr_neurodegen.tsv | awk '{print $1}')

zcat ${project_dir}/gwas_sumstats/lava_sumstats/gwas_sexstratified/${trait_M}.lava.gz > ${project_dir}/MR/meta_analysis/metal_input/${trait_M}.txt
zcat ${project_dir}/gwas_sumstats/lava_sumstats/gwas_sexstratified/${trait_F}.lava.gz > ${project_dir}/MR/meta_analysis/metal_input/${trait_F}.txt

metal ${project_dir}/MR/meta_analysis/metal_parameters/${trait_M}_${trait_F}.parameters.txt
