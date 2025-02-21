# Description: get list of synpase ID for immune proteins from UKB pQTLs

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# Read files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Downloaded SQL table from Synapse (filtered to keep only EUR proteins):
synapse_entity <- fread(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "synapse_entity_table.tsv")) %>%
  dplyr::rename(synapse_id = id)

olink_map <- fread(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "olink_protein_map_3k_v1.tsv")) %>%
  mutate(file_name = stringr::str_c(stringr::str_replace_all(UKBPPP_ProteinID, ":", "_"), "_", Panel, ".tar")) %>%
  inner_join(., synapse_entity, by = c("file_name" = "name"))
         
immport_genes <- fread(here(project_dir, "gwas_sumstats", "InnateDB_genes_ImmPort.txt"))

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

olink_immune <- inner_join(immport_genes, olink_map, by = c("ensembl" = "ensembl_id"))

write.table(olink_immune, here(project_dir, "gwas_sumstats", "pQTLs_UKB", "olink_map_with_immune_proteins.tsv"), sep = "\t", row.names = F, quote = F)
write.table(olink_immune %>% dplyr::select(synapse_id), here(project_dir, "gwas_sumstats", "pQTLs_UKB", "synapse_ids_for_1255_immune_proteins.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)



