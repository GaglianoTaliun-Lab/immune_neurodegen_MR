# Directionality test MR

# 1) calculate proportion of variance of liability for outcomes (r)
# 2) harmonise exposures and outcomes
# 3) run directionality_test(). Can provide a column with "r.outcome" in the harmonised data.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(TwoSampleMR)
library(MRPRESSO)
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

# colocalized signals:
coloc_sign_results <- try(read.table(here(project_dir, "results", "coloc", "summary_results", stringr::str_c("significant_PPH4_0.7_", pheno_outcome, "_immune_exposures.tsv")), sep = "\t", header = T))

if (inherits(coloc_sign_results, "try-error")) {
    stop("There are no significant coloc results. Exiting. \n", call. = FALSE)
}

exposure_files <- list.files(here(project_dir, "MR", "exposures", "ukb_pqtls"), full.names = F, pattern = "clumped_exposure_UKB_*")
sample_size_exp=33477

N_outcomes <- read.table(here(project_dir, "gwas_sample_sizes_sexstr_neurodegen.tsv"), sep = "\t", header = T)

if (nrow(coloc_sign_results) > 0) {

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

# loop across significant coloc results (PP H4 > 0.70) - outcome is the same only exposure varies:
    for (i in 1:nrow(coloc_sign_results)) { 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read exposure
 # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
        pheno_exposure = as.character(coloc_sign_results[i,"exposure"])
        exposure_data <- read.table(here(project_dir, "MR", "exposures", "ukb_pqtls", stringr::str_c("clumped_exposure_UKB_",pheno_exposure,".txt")), sep = "\t", header = T)
                    
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Harmonize data for directionality_test()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        dat_harmonised_full <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data, action = 2) %>%
            filter(., mr_keep == TRUE)

        if (nrow(dat_harmonised_full) > 4) {

            mrpresso_out <- try(mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", 
                    SdOutcome = "se.outcome", SdExposure = "se.exposure", 
                    OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat_harmonised_full, 
                    NbDistribution = 10000,  SignifThreshold = 0.05))

            if (!inherits(coloc_sign_results, "try-error")) {
                cat("Writing MR PRESSO results between ", pheno_exposure, " (exposure) and ", pheno_outcome, " (outcome)...\n")
                write.table(mrpresso_out, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("MR_PRESSO_pqtls_",pheno_exposure,"_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)
                cat("Done.\n")
            } else {
                cat("There was an error in mrpresso computation. Exiting without writing output file.\n")
            }

        } else {
            cat("There were less than five IVs in harmonised dataset, so MRPRESSO could not be performed.")
        }

    }
}
