# Description: compare allele frequencies of IVs used in cis MR main analysis between EUR and FIN.

# Load packages -----------------------------------------------------------

library(dplyr)
library(here)
library(stringr)
library(tidyr)
library(data.table)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# Load data ---------------------------------------------------------------

harmonised_files <- list.files(here(project_dir, "MR", "harmonised_data"), pattern = "data_harmonised_UKB_*",
  full.names = TRUE)

cat("Reading harmonised data...\n")
harmonised_snps <- lapply(harmonised_files, function(x) read.table(x, sep = "\t", header = TRUE) %>% 
    filter(., mr_keep == TRUE) %>% dplyr::select(SNP)) %>%
  rbindlist(.) %>% distinct(SNP)
cat("Done.\n")

cat("Reading allele frequencies...\n")
af_EUR <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_EUR_noFIN.frq")) %>%
  filter(., SNP %in% harmonised_snps$SNP) %>%
  mutate(MAF_EUR = as.numeric(MAF)) %>%
  dplyr::select(CHR, SNP, A1, A2, MAF_EUR)

af_FIN <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_FIN.frq")) %>%
  filter(., SNP %in% harmonised_snps$SNP) %>%
  mutate(MAF_FIN = as.numeric(MAF)) %>%
  dplyr::select(CHR, SNP, A1, A2, MAF_FIN)
cat("Done.\n")

###---------------------------------------- Main

cat("Joining allele frequency files...\n")
afs <- inner_join(af_EUR, af_FIN, by = c("CHR", "SNP", "A1", "A2")) %>%
  mutate(MAF_diff = abs(MAF_EUR - MAF_FIN))
  
afs_maf <- afs %>%
  filter(., MAF_diff > 0.2)
cat("Done.\n")

### ---------------------------------------------- Save files

cat("Writing output...\n")
# write a table of SNPs with MAF difference greater than 0.20:
write.table(afs_maf, here(project_dir, "MR", "MAF_difference_0.20_1KGP_EUR_FIN.txt"), sep = "\t", quote = F, row.names = F)
write.table(afs, here(project_dir, "MR", "MAF_difference_all_1KGP_EUR_FIN.txt"), sep = "\t", quote = F, row.names = F)
cat("Done.\n")