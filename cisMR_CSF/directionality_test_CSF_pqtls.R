# Directionality test MR

# 1) calculate proportion of variance of liability for outcomes (r)
# 2) harmonise exposures and outcomes
# 3) run directionality_test(). Can provide a column with "r.outcome" in the harmonised data.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

args <- commandArgs(TRUE)
pheno_outcome <- as.character(args[1])

dir_results <- list()

# MR results for 'pheno_outcome':
mr_sign_results <- read.table(here(project_dir, "MR", "results", "main",stringr::str_c("MR_main_results_immune_CSF_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
    filter(., method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    filter(., pval < 0.05)

exposure_files <- list.files(here(project_dir, "MR", "exposures", "CSF_pqtls"), full.names = F, pattern = "clumped_exposure_gene_*")
sample_size_exp=971

N_outcomes <- read.table(here(project_dir, "gwas_sample_sizes_sexstr_neurodegen.tsv"), sep = "\t", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load outcome and obtain r for directionality test, if there are significant results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (nrow(mr_sign_results) > 0) {

    if ((stringr::str_detect(pheno_outcome, "meta_sexcombined") == TRUE)){
        file_name = here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(pheno_outcome,"_wCHR_BP.tbl"))
    } else if ((stringr::str_detect(pheno_outcome, "males") == TRUE)){
        file_name = here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(pheno_outcome,".txt"))
    }

    # get # cases and controls:
    N_outcomes <- N_outcomes %>% filter(., trait == pheno_outcome)
    ncase = N_outcomes[1, "n_cases"]
    ncontrol = N_outcomes[1, "n_controls"]

    outcome_data <- read.table(file_name, sep = "\t", header = T)
    outcome_data <- outcome_data %>%
        filter(., !is.na(A1_AF)) %>% # Note: there are some variants with EAF in PD datasets that do not have SNP ID, and in the ALS datasets all in chr X (we do not have EAF for that CHR from source data - pQTLs do not have chr X data)
        mutate(outcome = pheno_outcome, id.outcome = "random_id", ncase.outcome = ncase, ncontrol.outcome = ncontrol) %>%
        dplyr::select(SNP,
                    effect_allele.outcome = A1,
                    other_allele.outcome = A2,
                    beta.outcome = BETA,
                    se.outcome = SE,
                    eaf.outcome = A1_AF,
                    pval.outcome = P,
                    samplesize.outcome = N,
                    outcome,
                    id.outcome,
                    ncase.outcome,
                    ncontrol.outcome)

    outcome_with_r <- get_r_from_lor(lor = outcome_data$beta.outcome, af = outcome_data$eaf.outcome, 
                                    ncase = outcome_data$ncase.outcome, ncontrol = outcome_data$ncontrol.outcome, prevalence = 0.1) %>%
                    cbind(outcome_data, .) %>%
                    rename(., r.outcome = .) %>%
                    select(SNP, r.outcome)

    # loop across significant (p-value < 0.05) results - outcome is the same only exposure varies:
    for (i in 1:nrow(mr_sign_results)) { 

        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Read exposure
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        
        pheno_exposure = as.character(mr_sign_results[i,"exposure"])
        exposure_data <- read.table(here(project_dir, "MR", "exposures", "CSF_pqtls", stringr::str_c("clumped_exposure_gene_",pheno_exposure,".txt")), sep = "\t", header = T)
                        
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        # Harmonize data for directionality_test()
        # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        dat_harmonised_full <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data, action = 2) %>%
            left_join(., outcome_with_r, by = "SNP") %>%
            mutate(samplesize.exposure = sample_size_exp)

        dir_results[[i]] <- directionality_test(dat_harmonised_full)

    }

    rbindlist(dir_results) %>%
        write.table(., here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("directionality_test_CSF_pqtls_", pheno_outcome, ".txt")), sep = "\t", quote = F)
}