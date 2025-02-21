# post-metal script to include CHR and BP into the output summary statistics, from the input summary statistics.

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
i <- as.numeric(args[1]) # which trait to process (1,3,5,7,9,11,13,15)

leukocyte_traits <- read.table(here(project_dir, "gwas_sample_sizes_sexstr_leukocytes.tsv"), sep = "\t", header = T)
trait <- leukocyte_traits[i, "trait"] %>% stringr::str_extract(., "[:alpha:]+_") %>% stringr::str_remove(., "_")

# females proportion
trait_fem = stringr::str_c(trait, "_females")
N_fem = leukocyte_traits %>% filter(., trait == trait_fem) %>% .$n_total

# males proportion
trait_mal = stringr::str_c(trait, "_males")
N_mal = leukocyte_traits %>% filter(., trait == trait_mal) %>% .$n_total

N_meta = N_fem + N_mal

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read files
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

metal_out <- read.table(here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(trait, "_meta_sexcombined1.tbl")), sep = "\t", header = T)
metal_inF <- read.table(here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(trait, "_females.txt")), sep = "\t", header = T) %>%
    dplyr::select(SNP, CHR, BP, A1, A2)
metal_inM <- read.table(here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(trait, "_males.txt")), sep = "\t", header = T) %>%
    dplyr::select(SNP, CHR, BP, A1, A2)

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
        P = P.value,
        N
    )

metal_input <- inner_join(metal_inF, metal_inM, by = c("SNP", "CHR", "BP", "A1", "A2")) %>%
    dplyr::select(SNP, CHR, BP, A1, A2)

output_full <- inner_join(metal_input, metal_out, by = c("SNP", "A1", "A2"))

write.table(output_full, here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(metalout_names[i],"_meta_sexcombined_wCHR_BP.tbl")), sep = "\t", quote = F, row.names = F)