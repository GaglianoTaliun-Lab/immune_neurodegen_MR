#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=1-0:00:00
#SBATCH --array=1-22
#SBATCH --job-name=getAF_plink
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=100G

module load StdEnv/2023
module load plink/2.00-20231024-avx2

project_dir="/home/fridald4/projects/def-gsarah/fridald4/als_immune"

plink2 \
   --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr${SLURM_ARRAY_TASK_ID}_v3.bgen 'ref-first' \
   --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample \
   --freq 'cols=chrom,pos,ref,alt1,alt,reffreq,alt1freq,altfreq' \
   --out /scratch/fridald4/als_immune/ukb_allele_freq/chr${SLURM_ARRAY_TASK_ID}