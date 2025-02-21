# Description: post-metal script to include CHR and BP and effect AF into the output summary statistics, from the input summary statistics.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(tidyverse)
library(stringr)
library(data.table)
library(here)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

args <- commandArgs(TRUE)
i <- as.numeric(args[1]) # which trait to process (1,3,5,7,9)

neurodegen_traits <- read.table(here(project_dir, "gwas_sample_sizes_sexstr_neurodegen.tsv"), sep = "\t", header = T)
trait <- neurodegen_traits[i, "trait"] %>% stringr::str_extract(., "[:alpha:]+_") %>% stringr::str_remove(., "_")

trait_sexcomb = stringr::str_c(trait, "_meta_sexcombined")
N_meta = neurodegen_traits %>% filter(., trait == trait_sexcomb) %>% .$n_total
trait_fem = stringr::str_c(trait, "_females")
trait_mal = stringr::str_c(trait, "_males")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

metal_out <- fread(here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(trait, "_meta_sexcombined1.tbl")))
metal_inF <- fread(here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(trait, "_females.txt"))) %>%
    dplyr::select(SNP, CHR, BP, A1, A2) %>%
    dplyr::distinct(., .keep_all = TRUE)
metal_inM <- fread(here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(trait, "_males.txt"))) %>%
    dplyr::select(SNP, CHR, BP, A1, A2) %>%
    dplyr::distinct(., .keep_all = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

metal_out <- metal_out %>%
    mutate(A1 = stringr::str_to_upper(Allele1), A2 = stringr::str_to_upper(Allele2), N = N_meta) %>%
    dplyr::select(
        SNP = MarkerName,
        A1,
        A2,
        BETA = Effect,
        SE = StdErr,
        A1_AF = Freq1,
        P = `P-value`,
        N
    )

metal_input <- full_join(metal_inM, metal_inF, by = c("SNP", "CHR", "BP", "A1", "A2")) %>%
    dplyr::select(-A1, -A2)

output_full <- inner_join(metal_input, metal_out, by = "SNP")

write.table(output_full, here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(trait,"_meta_sexcombined_wCHR_BP.tbl")), sep = "\t", quote = F, row.names = F)