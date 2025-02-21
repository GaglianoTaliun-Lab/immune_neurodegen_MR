# Two sample MR between ALS, AD and PD (outcomes, sex-stratified and sex-combined) and quantitative measures of leukocytes from UKB (per sex and sex-combined).
# Exposures were clumped in separate script using the UKB as reference.

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

mr_results = list()
mr_results_combined = list()
project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
mr_test_pairs <- read.table(here(project_dir, "MR", "MR_pairs_to_test.txt"), sep = "\t", header = T)
filter_fstat <- data.frame(N_before = rep(NA, nrow(mr_test_pairs)), N_after = rep(NA, nrow(mr_test_pairs)))

cat(stringr::str_c("There are a total of ", nrow(mr_test_pairs)," exposure-outcome pairs to test.\n"))

for (i in 1:nrow(mr_test_pairs)) { 

    pheno_exposure <- mr_test_pairs[i,1]
    pheno_outcome <- mr_test_pairs[i,2]
    type <- mr_test_pairs[i,3] # either sex-combined or sex-specific
    sample_size_exp <- mr_test_pairs[i,4]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read exposures
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cat(stringr::str_c("Reading exposure... \n"))

    exposure_data_as <- read.table(here(project_dir, "MR", "instruments", "clumped", stringr::str_c("clumped_exposure_noHLA_",pheno_exposure,"_",type,"_with_EAF.txt")), sep = "\t", header = T)
    exposure_data_chrx <- try(read.table(here(project_dir, "MR", "instruments", "clumped", stringr::str_c("clumped_exposure_chrX_",pheno_exposure,"_",type,"_with_EAF.txt")), sep = "\t", header = T), silent = TRUE)

    if (class(exposure_data_chrx) != "try-error") {
            exposure_data <- rbind(exposure_data_as, exposure_data_chrx) %>%
                mutate(id.exposure = exposure_data_as$id.exposure[1])
    } else {
        exposure_data <- exposure_data_as %>%
            mutate(id.exposure = exposure_data_as$id.exposure[1])
    }

    cat("Done.\n")

    # Keep exposures with F-statistics > 10 (script from https://github.com/globalbiobankmeta/multi-ancestry-pwmr/blob/main/R%20code/instrument.R)
    exposure_data <- cbind(exposure_data,fstatistics=1) # formatted exposure data

    filter_fstat[i,1] <- nrow(exposure_data)    

    for (s in 1:nrow(exposure_data)){
      
      z <- exposure_data[s,2]/exposure_data[s,3] # calculate Z-score with beta and SE
      pve <- z^2/(sample_size_exp+z^2) # proportion of variance explained
      # fstatistic: = (N-K-1)*PVE/(1-PVE); where N = sample size; K = N SNPs in that instrument (always = 1):
      exposure_data[s,"fstatistics"] <- (sample_size_exp-2)*pve/(1-pve)
    }

    exposure_data <- exposure_data[exposure_data$fstatistics>10,]
    filter_fstat[i,2] <- nrow(exposure_data)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load outcomes
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (type == "sex_combined"){
        file_name = here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(pheno_outcome,"_wCHR_BP.tbl"))
    } else if (type == "sex_specific"){
        file_name = here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(pheno_outcome,".txt"))
    }

    outcome_data <- read_outcome_data(
        snps = exposure_data$SNP,
        filename = file_name,
        sep = "\t",
        snp_col = "SNP",
        beta_col = "BETA",
        se_col = "SE",
        effect_allele_col = "A1",
        other_allele_col = "A2",
        eaf_col = "A1_AF",
        pval_col = "P"
    ) %>%
    filter(., pval.outcome >= 5e-08) %>%
    mutate(outcome = pheno_outcome)
        
    cat("Done.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Harmonize data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cat(stringr::str_c("Harmonizing data for traits in line #", i,"...\n"))
    dat_harmonised <- harmonise_data(exposure_dat = exposure_data, outcome_dat = outcome_data, action = 2) 

    cat("Done.\n")
    
    dat_harmonised$SNP %>%
        write.table(., here(project_dir, "MR", "instruments", stringr::str_c("SNPs_for_HP_test_noHLA_",pheno_exposure,"_",pheno_outcome,".txt")), sep = "\t", row.names = F, col.names = F, quote = F)
        # HP = horizontal pleiotropy
    
    write.table(dat_harmonised, here(project_dir, "MR", "harmonised_data", stringr::str_c("data_harmonised_noHLA_",pheno_exposure,"_",pheno_outcome,".txt")), sep = "\t", row.names = F, quote = F)
    
    cat(stringr::str_c("There are ", nrow(dat_harmonised), " rows in harmonised data table.\n"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Perform MR
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cat(stringr::str_c("Performing MR for traits in line #", i,"...\n"))
    mr_results[[i]] <- mr(dat_harmonised, method_list = c("mr_ivw_mre", "mr_egger_regression", "mr_weighted_median"))
    # other methods that can be included: "mr_weighted_median", "mr_weighted_mode"
    cat("Done.\n")

    cat(stringr::str_c("There are ", nrow(mr_results[[i]]), " rows in MR results table.\n"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SENSITIVITY TESTS
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # a) Test for heterogeneity:
    cat("Calculating heterogeneity...\n")
    het <- mr_heterogeneity(dat_harmonised)
    write.table(het, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("heterogeneity_MR_",pheno_exposure,"_",pheno_outcome,"_excl_conf_excl_HLA.txt")), sep = "\t", row.names = F, quote = F)
    cat("Done.\n")

    # b) test for horizontal pleiotropy using MR-Egger intercept test:
    cat("Testing for horizontal pleiotropy...\n")
    plt <- mr_pleiotropy_test(dat_harmonised)
    write.table(plt, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("pleiotropy_test_",pheno_exposure,"_",pheno_outcome,"_excl_conf_excl_HLA.txt")), sep = "\t", row.names = F, quote = F)
    cat("Done.\n")

    # c) Leave one out analysis
    cat("Performing leave one out analysis...\n")
    mr_leaveoneout(dat_harmonised) %>%
        write.table(., here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("LeaveOneOut_",pheno_exposure,"_",pheno_outcome,"_excl_conf_excl_HLA.txt")), sep = "\t", row.names = F, quote = F)
    cat("Done.\n")

    # d) Single SNP analysis
    cat("Performing single SNP analysis...\n")
    sin <- mr_singlesnp(dat_harmonised)
    write.table(sin, here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("SingleSNP_",pheno_exposure,"_",pheno_outcome,"_excl_conf_excl_HLA.txt")), sep = "\t", row.names = F, quote = F)
    cat("Done.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    cat("Generating plots...\n")

    # Scatter plot (slopes are the methods used in results)
    p1 <- mr_scatter_plot(mr_results[[i]], dat_harmonised)

    p2 <- mr_funnel_plot(mr_singlesnp(dat_harmonised, all_method=c("mr_ivw", "mr_egger_regression")))

    p3 <- mr_forest_plot(mr_singlesnp(dat_harmonised, all_method=c("mr_ivw", "mr_egger_regression")))

    plots <- ggarrange(p1[[1]],p2[[1]],p3[[1]], nrow=1, ncol=1)

    ggexport(plots, filename=here(project_dir, "MR", "results","figures",stringr::str_c("plots_",pheno_exposure,"_",pheno_outcome,"_excl_conf_excl_HLA.pdf")))    

    cat("Done.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Output
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# combine all results
    mr_results_combined[[i]] <- combine_all_mrresults(
        mr_results[[i]],
        het,
        plt,
        sin,
        ao_slc = FALSE,
        Exp = TRUE,
    )

}

cat("All MR tests have been performed. Saving results...\n")
mr_results_df <- rbindlist(mr_results)
mr_results_combdf <- rbindlist(mr_results_combined)
write.table(mr_results_df, here(project_dir, "MR", "results", "main", "MR_main_results_all_pairs_excl_conf_excl_HLA.txt"), sep = "\t", row.names = F, quote = F)
write.table(mr_results_combdf, here(project_dir, "MR", "results", "main", "MR_main_results_combined_all_pairs_excl_conf_excl_HLA.txt"), sep = "\t", row.names = F, quote = F)
write.table(filter_fstat, here(project_dir, "MR", "results", "sensitivity_tests", "Fstatistics_N_instruments.txt"), sep = "\t", row.names = F, quote = F)
cat("Done.\n")

