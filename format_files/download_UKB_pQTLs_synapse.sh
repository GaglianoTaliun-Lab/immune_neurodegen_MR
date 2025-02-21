
nrow=$(wc -l /home/fridald4/projects/def-gsarah/fridald4/als_immune/gwas_sumstats/pQTLs_UKB/synapse_ids_for_1255_immune_proteins.tsv | awk '{print $1}')
synapse login

for ((i=1; i<=${nrow};i++))
do
    echo ${i}
    synapse_id=$(sed -n ${i}p /home/fridald4/projects/def-gsarah/fridald4/als_immune/gwas_sumstats/pQTLs_UKB/synapse_ids_for_1255_immune_proteins.tsv)
    echo "Processing synapse id" ${synapse_id} "..."
    synapse get -r ${synapse_id}
    echo "Done."
done

# extract tar files:
ls *.tar > tar_files.txt
nrow=$(wc -l /home/fridald4/projects/def-gsarah/fridald4/als_immune/gwas_sumstats/pQTLs_UKB/tar_files.txt | awk '{print $1}')

for ((i=1; i<=${nrow};i++))
do
    file=$(sed -n ${i}p /scratch/fridald4/als_immune/ukb_pqtls/tar_files.txt | awk '{print $1}')
    tar xf /scratch/fridald4/als_immune/ukb_pqtls/$file
    echo ${file} "done."
done