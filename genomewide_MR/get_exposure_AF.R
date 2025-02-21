# script to insert UKB allele frequencies into the clumped exposure files for MR.

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
scratch_dir <- "/scratch/fridald4/als_immune"
snps_with_AF <- list()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# UKB allele frequencies from plink
all_freq_files <- list.files(here(scratch_dir, "ukb_allele_freq"), pattern = "*.afreq", full.names = T)

all_freqs <- lapply(all_freq_files, function(x) sub("#", "", readLines(x)) %>% 
                    read.table(text = ., sep = "\t", header = T)) %>%
    rbindlist(.) %>%
    rename(SNP = ID)

# clumped exposures files
clumped_files <- list.files(here(project_dir, "MR", "instruments", "clumped"), pattern = "clumped_exposure_noHLA_", full.names = T) %>% sort()
clumped_snps <- lapply(clumped_files, function(x) read.table(x, sep = "\t", header = T))

clumped_names <- list.files(here(project_dir, "MR", "instruments", "clumped"), pattern = "clumped_exposure_noHLA_", full.names = F) %>% sort() %>%
      stringr::str_remove(., ".txt")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Insert AF into each exposure file clumped
clumped_files_af <- lapply(clumped_snps, function(x) left_join(x, all_freqs, by = "SNP") %>%
                           mutate(eaf.exposure =
                                  case_when(
                                        effect_allele.exposure == REF ~ REF_FREQ,
                                        effect_allele.exposure == ALT1 ~ 1 - REF_FREQ,
                                        TRUE ~ NA
                                       )
                                 ) %>%
                                       dplyr::select(
                                        SNP, beta.exposure, se.exposure, effect_allele.exposure,
                                        other_allele.exposure, pval.exposure, exposure,
                                        mr_keep.exposure, pval_origin.exposure, id.exposure, eaf.exposure
                                       )
                           )

# write new clumped files with AF
for (i in 1:length(clumped_files_af)){

    write.table(clumped_files_af[[i]], here(project_dir, "MR", "instruments", "clumped",str_c(clumped_names[[i]], "_with_EAF.txt")),
                sep = "\t", row.names = F, quote = F)

 }