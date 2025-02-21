#!/bin/bash

#SBATCH --account=def-gsarah
#SBATCH --time=8:00:00
#SBATCH --job-name=getAF_chrX
#SBATCH --output=slurm-%x.out
#SBATCH --error=slurm-%x.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=3
#SBATCH --mem-per-cpu=100G

module load StdEnv/2023
module load plink/2.00-20231024-avx2

project_dir="/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# sex-combined
plink2 \
    --bgen /lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chrX_v3.bgen 'ref-first' \
    --sample /home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/chrX/ukb22828_cX_b0_v3_s486556.sample \
    --freq 'cols=chrom,pos,ref,alt1,alt,reffreq,alt1freq,altfreq' \
    --out /scratch/fridald4/als_immune/ukb_allele_freq/chr23