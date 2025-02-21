#!/bin/bash
#SBATCH --account=def-gsarah
#SBATCH --time=4:00:00
#SBATCH --array=2-16
#SBATCH --job-name=ldsc_all
#SBATCH --output=slurm-%x-%a.out
#SBATCH --error=slurm-%x-%a.err
#SBATCH --mail-user=frida.lona-durazo@icm-mhi.org
#SBATCH --mail-type=ALL
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=4G

module load StdEnv/2020
module load python/2.7
module load scipy-stack/2020a

virtualenv --no-download $SLURM_TMPDIR/env
source $SLURM_TMPDIR/env/bin/activate
pip install --upgrade pip
# pip install bitarray
# pip install pybedtools==0.7.10
pip install contextlib2==0.2
pip install importlib_metadata==1.5.0
pip install jsonschema==2.6.0
pip install six==1.16.0+computecanada
pip install bitarray==0.8.0

project_dir="/home/fridald4/projects/def-gsarah/fridald4/als_immune"

pip install -r ${project_dir}/tools/ldsc/CC_requirements.txt

pheno1=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${project_dir}/gwas_sample_sizes_sexstr_neurodegen.tsv | awk '{print $1}')
n_lines=$(wc -l ${project_dir}/gwas_sample_sizes_sexstr_neurodegen.tsv | awk '{print $1}')

for ((i=2;i<=n_lines;i++))
do
	pheno2=$(sed -n ${i}p ${project_dir}/gwas_sample_sizes_sexstr_neurodegen.tsv | awk '{print $1}')
	python ${project_dir}/tools/ldsc/ldsc.py \
		--ref-ld-chr ${project_dir}/reference_data/eur_w_ld_chr/ \
		--out ${project_dir}/results/ldsc_corr/ldsc_output_sexstr/${pheno1}_${pheno2}_rg \
		--rg ${project_dir}/gwas_sumstats/ldsc_sumstats/${pheno1}.sumstats.gz,${project_dir}/gwas_sumstats/ldsc_sumstats/${pheno2}.sumstats.gz \
		--w-ld-chr ${project_dir}/reference_data/eur_w_ld_chr/
done
