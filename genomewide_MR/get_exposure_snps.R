# script to obtain unique SNPs from clumped exposure files, so to obtain ukb allele frequencies for MR.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
snps_with_AF <- list()
source(here(project_dir, "scripts", "get_AF_plink_fnc.R"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clumped_all <- list.files(here(project_dir, "MR", "instruments", "clumped"), pattern = "clumped_exposure_noHLA_*", full.names = T) %>% sort()
clumped_names <- list.files(here(project_dir, "MR", "instruments", "clumped"), pattern = "clumped_exposure_noHLA_*", full.names = F) %>% sort()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get unique SNPs from all clumped exposures
clumped_files <- lapply(clumped_all, function(x) read.table(x, sep = "\t", header = T))
clumped_snps <- lapply(clumped_all, function(x) read.table(x, sep = "\t", header = T) %>%
                        dplyr::select(SNP))

clumped_snps_uniq <- rbindlist(clumped_snps) %>% distinct()

write.table(clumped_snps_uniq, here(project_dir, "MR", "instruments", "clumped", "ukb_allele_freq", "input_plink_UKB_AF.txt"), sep = "\t", row.names = F, col.names = F, quote = F)
