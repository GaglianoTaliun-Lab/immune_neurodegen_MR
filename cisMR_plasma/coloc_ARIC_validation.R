# Description: run colocalization analysis with ARIC data for validation

# Packages -------------------------------------------------------

library(coloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(BSgenome)
library(colochelpR)
library(vctrs)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

args <- commandArgs(TRUE)
pheno_outcome <- as.character(args[1])

coloc_results_summ <- list()
results_names <- array()

# Coloc priors - less stringent due to already evidence in MR
p1 = 1e-03
p2 = 1e-03
p12 = 1e-04

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load MR results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# directionality test results:
dir_test <- try(read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("directionality_test_pqtls_aric_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
    dplyr::select(exposure, correct_causal_direction))

if (!inherits(dir_test, "try-error")) {

    # MR Egger test results:
    mregger <- try(read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_Egger_main_results_immune_ARICvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        dplyr::select(exposure, MREgger_regression_pval = pval))

    if (inherits(mregger, "try-error")) {
        mr_sign_results <- read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_main_results_immune_ARICvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        filter(., method %in% c("Inverse variance weighted", "Wald ratio")) %>%
        left_join(., dir_test, by = "exposure") %>%
        filter(., pval < 0.005) %>%
        filter(., correct_causal_direction == "TRUE" | is.na(correct_causal_direction)) #%>%
        #filter(., heterogeneity_signal == "NO" | is.na(heterogeneity_signal))

    } else {

        # MR results for pheno_outcome (IVW), keeping:
            # 1) pvalue < 0.005
            # 2) MR Egger intercept not significant (i.e., no pleiotropy)
            # 4) No significant heterogeneity
        mr_sign_results <- read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_main_results_immune_ARICvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
            filter(., method %in% c("Inverse variance weighted", "Wald ratio")) %>%
            left_join(., dir_test, by = "exposure") %>%
            filter(., pval < 0.005) %>%
            filter(., correct_causal_direction == "TRUE" | is.na(correct_causal_direction)) #%>%
            #filter(., heterogeneity_signal == "NO" | is.na(heterogeneity_signal))

        if (nrow(mregger) > 0) {
            mr_sign_results <- mr_sign_results %>%
                left_join(., mregger, by = "exposure") %>%
                filter(., MR_Egger_pleiotropy == "NO" | is.na(MR_Egger_pleiotropy))
        }
    }              

    # load case control sample sizes:
    N_outcomes <- read.table(here(project_dir, "gwas_sample_sizes_sexstr_neurodegen.tsv"), sep = "\t", header = T)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load outcome
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if ((stringr::str_detect(pheno_outcome, "meta_sexcombined") == TRUE)){
            file_name = here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(pheno_outcome,"_wCHR_BP.tbl"))
        } else if ((stringr::str_detect(pheno_outcome, "males") == TRUE)){
            file_name = here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(pheno_outcome,".txt"))
        }

    outcome_data <- read.table(file_name, sep = "\t", header = T) %>%
        filter(., !is.na(A1_AF)) %>%
        mutate(MAF = case_when(
                    A1_AF <= 0.5 ~ A1_AF,
                    A1_AF > 0.5 ~ 1 - A1_AF),
            outcome = pheno_outcome,
            varbeta = SE^2)

    # get # cases and controls:
    N_outcomes <- N_outcomes %>% filter(., trait == pheno_outcome)
    ncase = N_outcomes[1, "n_cases"]
    ncontrol = N_outcomes[1, "n_controls"]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read exposure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (nrow(mr_sign_results) > 0) {

        files <- setNames(list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Zhang2022","EA"), 
                             full.names = T, pattern = "*PHENO1.glm.linear"),
                  nm = list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Zhang2022","EA"), 
                                  full.names = F, pattern = "*PHENO1.glm.linear") %>%
                    stringr::str_remove_all(., ".PHENO1.glm.linear"))

        seq_ids <- read.table(here(project_dir,"gwas_sumstats","pQTLs_Zhang2022", "seqid.txt"),
            sep = "\t", header = T) %>%
            dplyr::rename(gene_name = entrezgenesymbol)

        pqtls_all <- lapply(files, function(x) read.table(x, sep = "\t", header= F, skip = 1, 
                    col.names = c("CHR", "BP", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")) %>%
                    dplyr::mutate(P = as.numeric(P)))

        pqtls_all <- pqtls_all %>%
            dplyr::bind_rows(., .id = "seqid_in_sample") %>%
            left_join(., seq_ids, by = "seqid_in_sample") %>% 
            filter(., gene_name %in% mr_sign_results$exposure)
        colnames(pqtls_all) <- c("seqid_in_sample","CHR","BP","ID","REF","ALT","A1","A1_FREQ","TEST","OBS_CT",
                                "BETA","SE","T_STAT","P","ERRCODE","uniprot_id","gene_name", "chromosome_name", "transcription_start_site")
        
        for (i in 1:nrow(mr_sign_results)) {

            gene = mr_sign_results[i, "exposure"]
            cat("The gene being assessed is ", gene, ".\n")

            exposure_gene_data <- pqtls_all %>%
                filter(., gene_name == gene) %>%
                dplyr::mutate(N = 7213, varbeta = SE^2,
                              MAF = case_when(
                                A1_FREQ > 0.5 ~ 1-A1_FREQ,
                                A1_FREQ <= 0.5 ~ A1_FREQ
                              )) %>%
                filter(., str_detect(ID, "rs")) %>%
                filter(., !is.na(P)) %>%
                filter(., !is.na(MAF) | !is.na(ID)) %>%
                filter(., MAF > 0) %>%
                distinct(., ID, .keep_all = TRUE) %>%
                dplyr::select(SNP = ID,
                              beta = BETA,
                              varbeta,
                              pvalues = P,
                              MAF,
                              N)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run coloc
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # outcome
            df1 <- outcome_data %>%
                        filter(., SNP %in% exposure_gene_data$SNP) %>%
                        filter(., MAF > 0) %>%
                        filter(., BETA != 0) %>%
                        filter(., !is.na(BETA)) %>%
                        distinct(., SNP, .keep_all = TRUE) %>%
                        mutate(pvalues = as.numeric(P),
                            s = ncase/N)

            # exposure (pQTLs)
            df2 <- exposure_gene_data %>%
                        filter(., SNP %in% df1$SNP)

            coloc_results <- coloc.abf(dataset1 = list(type = "cc",
                                                    snp = df1$SNP,
                                                    beta = df1$BETA,
                                                    varbeta = df1$varbeta,
                                                    pvalues = df1$pvalues,
                                                    MAF = df1$MAF,
                                                    N = df1$N,
                                                    s = df1$s[1]),
                                    dataset2 = list(type = "quant",
                                                    snp = df2$SNP,
                                                    beta = df2$beta,
                                                    varbeta = df2$varbeta,
                                                    pvalues = df2$pvalues,
                                                    MAF = df2$MAF,
                                                    N = df2$N),
                                    p1 = p1, p2 = p2, p12 = p12
            )

            coloc_results_summ[[i]] <- coloc_results$summary
            results_names[i] <- gene

        }

        coloc_results_summ <- setNames(coloc_results_summ, nm = results_names)

        coloc_summary_out <- coloc_results_summ %>% dplyr::bind_rows(., .id = "exposure")
            write.table(coloc_summary_out, here(project_dir, "results", "coloc", "summary_results", str_c(pheno_outcome, "_aric_validation_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        mr_sign_results %>%
            left_join(., coloc_summary_out, by = "exposure") %>%
            dplyr::rename(coloc_nsnps = nsnps) %>%
            write.table(., here(project_dir, "MR", "results", "validation", stringr::str_c("MR_significant_main_results_immune_ARICvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", row.names = F, quote = F)

    } else { cat("There are no significant MR results. Ending without outputting anything.")}

} else { cat("There are no directionality results, therefore no analysis to perform.\n")}