# Two sample MR to validate plasma pQTL results in UKB pQTLs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load libraries
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(TwoSampleMR)
library(ieugwasr)
library(dplyr)
library(tidyr)
library(here)
library(stringr)
library(data.table)
library(ggpubr)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

mr_ivw_results = list()
mr_egger_results = list()
mr_results_combined = list()

pheno_outcome <- "PD_meta_sexcombined"
type <- "sex_combined"

cat("Outcome: ", pheno_outcome, "\n")
cat("Type: ", type, "\n")

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
exposure_files <- list.files(here(project_dir, "MR", "exposures", "ukb_pqtls"), full.names = F, pattern = "clumped_exposure_gene_*")

filter_fstat <- data.frame(gene = rep(NA, length(exposure_files)), N_before = rep(NA, length(exposure_files)), N_after = rep(NA, length(exposure_files)))
sample_size_exp=33477

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read exposures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:length(exposure_files)) { 

    gene = exposure_files[i] %>% stringr::str_remove_all(., ".txt") %>% stringr::str_remove_all(., "clumped_exposure_gene_")

    cat(stringr::str_c("Reading exposure for ", gene, ".\n"))
    exposure_data <- read.table(here(project_dir, "MR", "exposures", "ukb_pqtls", exposure_files[i]), sep = "\t", header = T) 
    cat("Done.\n")
    
    # Keep exposures with F-statistics > 10 (script from https://github.com/globalbiobankmeta/multi-ancestry-pwmr/blob/main/R%20code/instrument.R)
    exposure_data <- cbind(exposure_data,fstatistics=1) # formatted exposure data

    filter_fstat[i,1] <- gene
    filter_fstat[i,2] <- nrow(exposure_data)    

    for (s in 1:nrow(exposure_data)){
      
      z <- exposure_data[s,"beta.exposure"]/exposure_data[s,"se.exposure"] # calculate Z-score with beta and SE
      pve <- z^2/(sample_size_exp+z^2) # proportion of variance explained
      # fstatistic: = (N-K-1)*PVE/(1-PVE); where N = sample size; K = N SNPs in that instrument (always = 1):
      exposure_data[s,"fstatistics"] <- (sample_size_exp-2)*pve/(1-pve)
    }

    exposure_data <- exposure_data[exposure_data$fstatistics>10,]
    filter_fstat[i,3] <- nrow(exposure_data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load outcomes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cat(stringr::str_c("Reading outcome...\n"))
    if (type %in% c("sex_combined", "sex_specific")){

        outcome_data <- try(read_outcome_data(
            snps = exposure_data$SNP,
            filename = here(project_dir, "MR", "outcomes", stringr::str_c(pheno_outcome,"_withAF_for_MR.txt")),
            sep = "\t",
            snp_col = "SNP",
            beta_col = "BETA",
            se_col = "SE",
            effect_allele_col = "A1",
            other_allele_col = "A2",
            eaf_col = "eaf",
            ncase_col = "ncase",
            ncontrol_col = "ncontrol",
            pval_col = "P"
        ) %>%
        mutate(outcome = pheno_outcome)
        )

        if  (inherits(outcome_data, "try-error")) {
            cat("No common instruments in the outcome data. Next.\n")
            next
        }

        cat("Done. Kept ", nrow(outcome_data), "instruments out of ", nrow(exposure_data), " in exposure data.\n")

    } else {
        cat("Error reading outcome.\n")
    }

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Harmonize data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (nrow(outcome_data > 0)) {

        cat(stringr::str_c("Harmonizing data for traits in line #", i,"...\n"))
        dat_harmonised <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data, action = 2) 
        # action = 2 means that it will try to infer positive strand alleles using the AF for palindromes.
        cat("Done.\n")
        
        dat_harmonised$SNP %>%
            write.table(., here(project_dir, "MR", "instruments", "pqtls", stringr::str_c("SNPs_for_HP_test_",gene,"_",pheno_outcome,".txt")), sep = "\t", row.names = F, col.names = F, quote = F)
            # HP = horizontal pleiotropy
        
        write.table(dat_harmonised, here(project_dir, "MR", "harmonised_data", stringr::str_c("data_harmonised_UKB_",gene,"_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)
        
        cat(stringr::str_c("There are ", nrow(dat_harmonised), " rows in harmonised data table.\n"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform MR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        cat(stringr::str_c("Performing MR for traits in line #", i,"...\n"))
        if (nrow(dat_harmonised) <= 1) {
            mr_ivw_results[[i]] <- mr(dat_harmonised, method_list = c("mr_wald_ratio"))
        } else {
            mr_ivw_results[[i]] <- mr(dat_harmonised, method_list = c("mr_ivw"))
        }

        mr_egger_results[[i]] <- mr(dat_harmonised, method_list = c("mr_egger_regression"))

        # Include MR Egger regression pvalue in main IVW results only if there are 3 or more genetic instruments
        if(nrow(mr_egger_results[[i]]) > 0) {
            if (mr_egger_results[[i]]$nsnp > 2) {
            mr_ivw_results[[i]]$MREgger_regression_pval <- mr_egger_results[[i]]$pval[i]
            }
        }

        cat("Done.\n")

        cat(stringr::str_c("There are ", nrow(mr_ivw_results[[i]]), " rows in MR results table.\n"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sensitivity tests for nominally significant results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        if (nrow(mr_ivw_results[[i]]) > 0) {
            if (mr_ivw_results[[i]]$method %in% c("Inverse variance weighted", "Wald ratio") & mr_ivw_results[[i]]$pval < 0.05) {
            
                # a) Test for heterogeneity:
                cat("Calculating heterogeneity...\n")

                if (nrow(dat_harmonised) <= 1) {
                mr_ivw_results[[i]]$het_Qpval <- "NA"
                mr_ivw_results[[i]]$heterogeneity_signal <- "NA"
                } else {
                    het <- mr_heterogeneity(dat_harmonised, method_list = c("mr_ivw"))
                    mr_ivw_results[[i]]$het_Qpval <- het$Q_pval[1]
                    mr_ivw_results[[i]]$heterogeneity_signal <- ifelse(het$Q_pval[1] < 0.05, "YES", "NO")
                    write.table(het, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("heterogeneity_MR_",gene,"_",pheno_outcome,"_pqtls.txt")), sep = "\t", row.names = F, quote = F)
                        cat("Done.\n")
                }

                # b) test for horizontal pleiotropy using MR-Egger intercept test:
                if (nrow(dat_harmonised) <= 1) {
                    mr_ivw_results[[i]]$MR_Egger_intercept_pval <- "NA"
                    mr_ivw_results[[i]]$MR_Egger_pleiotropy <- "NA"
                } else {
                    cat("Testing for horizontal pleiotropy...\n")
                    plt <- mr_pleiotropy_test(dat_harmonised)
                    mr_ivw_results[[i]]$MR_Egger_intercept_pval <- plt$pval[1]
                    mr_ivw_results[[i]]$MR_Egger_pleiotropy <- ifelse(plt$pval[1] < 0.05, "YES", "NO")
                    write.table(plt, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("pleiotropy_test_",gene,"_",pheno_outcome,"_pqtls.txt")), sep = "\t", row.names = F, quote = F)
                        cat("Done.\n")
                }

                # c) Leave one out analysis
                cat("Performing leave one out analysis...\n")
                mr_leaveoneout(dat_harmonised) %>%
                    write.table(., here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("LeaveOneOut_",gene,"_",pheno_outcome,"_CSF_pqtls.txt")), sep = "\t", row.names = F, quote = F)
                cat("Done.\n")

                # d) Single SNP analysis
                cat("Performing single SNP analysis...\n")
                sin <- mr_singlesnp(dat_harmonised)
                write.table(sin, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("SingleSNP_",gene,"_",pheno_outcome,"_CSF_pqtls.txt")), sep = "\t", row.names = F, quote = F)
                cat("Done.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

                cat("Generating plots...\n")

                # Scatter plot (slopes are the methods used in results)
                p1 <- mr_scatter_plot(mr_ivw_results[[i]], dat_harmonised)

                p2 <- mr_funnel_plot(mr_singlesnp(dat_harmonised, all_method=c("mr_ivw", "mr_egger_regression", "mr_wald_ratio")))

                p3 <- mr_forest_plot(mr_singlesnp(dat_harmonised, all_method=c("mr_ivw", "mr_egger_regression", "mr_wald_ratio")))

                plots <- ggarrange(p1[[1]],p2[[1]],p3[[1]], nrow=1, ncol=1)

                ggexport(plots, filename=here(project_dir, "MR", "results","figures",stringr::str_c("plots_UKB_pqtls_",gene,"_",pheno_outcome,".pdf")))    

                cat("Done.\n")

            }
        }
        
    } else if (nrow(outcome_data < 1)) {

        cat("There were 0 remaining instruments after loading ", pheno_outcome, "outcome data for gene ", gene, ".\n")
    }
}

cat("All MR tests have been performed. Saving results...\n")

mr_results_df <- rbindlist(mr_ivw_results, fill = TRUE) %>%
    mutate(FDR = p.adjust(pval, method = "BH"))

mr_results_egger <- rbindlist(mr_egger_results, fill = TRUE) %>%
    mutate(FDR = p.adjust(pval, method = "BH"))

write.table(mr_results_df, here(project_dir, "MR", "results", "main", stringr::str_c("MR_main_results_immune_UKBvalidation_pqtls_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)

write.table(mr_results_egger, here(project_dir, "MR", "results", "main", stringr::str_c("MR_Egger_main_results_immune_UKBvalidation_pqtls_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)

filter_fstat %>% mutate(N_instruments_removed = N_before - N_after) %>%
    write.table(., here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("Fstatistics_N_instruments_UKBvalidation_pqtls_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)

cat("Done.\n")

