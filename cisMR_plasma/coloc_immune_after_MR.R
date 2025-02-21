# Description: run colocalization analysis for MR with UKB plasma pQTLs significant results

# Packages -------------------------------------------------------

library(coloc)
library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Set arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

args <- commandArgs(TRUE)
pheno_outcome <- as.character(args[1])

coloc_results_summ <- list()
coloc_results_res <- list()
results_names <- array()

# Coloc priors - less stringent due to already evidence in MR
p1 = 1e-03
p2 = 1e-03
p12 = 1e-04

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load MR results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# directionality test results:
dir_test <- try(read.table(here(project_dir, "MR", "results", "sensitivity_tests", stringr::str_c("directionality_test_pqtls_UKB_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
    dplyr::select(exposure, correct_causal_direction))

# if there is no error in trying to read file (i.e., it does exist, then continue):
if (!inherits(dir_test, "try-error")) {

    # MR Egger test results:
    mregger <- read.table(here(project_dir, "MR", "results", "main", stringr::str_c("MR_Egger_main_results_immune_pqtls_UKB_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        dplyr::select(exposure, MREgger_regression_pval = pval)

    # MR results for pheno_outcome (IVW), keeping:
        # 1) FDR < 0.05
        # 2) MR Egger intercept not significant (i.e., no pleiotropy)
        # 3) Directionality test = TRUE
        # 4) No significant heterogeneity
    mr_sign_results <- read.table(here(project_dir, "MR", "results", "main",stringr::str_c("MR_main_results_immune_pqtls_UKB_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        filter(., method %in% c("Inverse variance weighted", "Wald ratio")) %>%
        left_join(., dir_test, by = "exposure") %>%
        left_join(., mregger, by = "exposure") %>%
        filter(., FDR < 0.05) %>%
        filter(., correct_causal_direction == "TRUE" | is.na(correct_causal_direction)) %>%
        filter(., MR_Egger_pleiotropy == "NO" | is.na(MR_Egger_pleiotropy)) %>%
        filter(., heterogeneity_signal == "NO" | is.na(heterogeneity_signal))

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
    # Read exposure (same as in formatting pQTLs for MR, but without FDR filtering)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    # Olink map with 1,248 immune-related proteins
    olink_map <- as.data.frame(fread(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "olink_map_with_immune_proteins.tsv"))) %>%
        mutate(pheno_id = stringr::str_c(HGNC.symbol , "_", OlinkID,"_",UniProt))

    # loop across significant results (as defined above) - outcome is the same only exposure varies:

    if (nrow(mr_sign_results) > 0) {

        for (i in 1:nrow(mr_sign_results)) { 

            pheno_exposure = as.character(mr_sign_results[i,"exposure"])

            olink_map_sub <- olink_map %>%
                filter(., pheno_id == pheno_exposure)

            chr = olink_map_sub[1, "chr"]
            gene_start = olink_map_sub[1, "gene_start"]
            gene_end = olink_map_sub[1, "gene_end"]
            name_id = olink_map_sub[1, "UKBPPP_ProteinID"]
            Panel = olink_map_sub[1, "Panel"]
            directory_name = olink_map_sub[1, "file_name"] %>% stringr::str_remove(., ".tar")
            file_name = stringr::str_c("discovery_chr", chr, "_", name_id, ":", Panel, ".gz")

            rsid_map <- fread(stringr::str_c("/scratch/fridald4/als_immune2/ukb_pqtls/olink_rsid_map_mac5_info03_b0_7_chr", chr, "_patched_v2.tsv.gz")) %>%
                dplyr::select(rsid, GENPOS = POS38)

            exposure_gene_data <- fread(stringr::str_c("/scratch/fridald4/als_immune2/ukb_pqtls/", directory_name, "/", file_name)) %>%
                mutate(pvalues = as.numeric(10^-LOG10P),
                    MAF = case_when(
                    A1FREQ <= 0.5 ~ A1FREQ,
                    A1FREQ > 0.5 ~ 1 - A1FREQ),
                    varbeta = SE^2) %>%
                filter(., !is.na(pvalues)) %>%
                filter(CHROM == chr & GENPOS >= gene_start - 500000 & GENPOS <= gene_end + 500000) %>%
                left_join(., rsid_map, by = "GENPOS") %>%
                dplyr::select(SNP = rsid,
                            effect_allele = ALLELE1,
                            other_allele = ALLELE0,
                            MAF,
                            BETA,
                            varbeta,
                            pvalues,
                            N) %>%
                distinct(., SNP, .keep_all = TRUE) %>%
                filter(., !is.na(MAF)) %>%
                filter(., !is.na(SNP)) %>%
                filter(., MAF > 0)
            
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Run coloc
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            # outcome (either ALS, AD or PD)
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
                                                    beta = df2$BETA,
                                                    varbeta = df2$varbeta,
                                                    pvalues = df2$pvalues,
                                                    MAF = df2$MAF,
                                                    N = df2$N),
                                    p1 = p1, p2 = p2, p12 = p12
            )

            coloc_results_summ[[i]] <- coloc_results$summary
            coloc_results_res[[i]] <- coloc_results$results
            results_names[i] <- pheno_exposure
            
        }

        coloc_results_summ <- setNames(coloc_results_summ, nm = results_names)
        coloc_results_res <- setNames(coloc_results_res, nm = results_names)

        coloc_results_out <- coloc_results_summ %>% dplyr::bind_rows(., .id = "exposure")
        write.table(coloc_results_out, here(project_dir, "results", "coloc", "summary_results", str_c(pheno_outcome, "_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        sign_coloc <- coloc_results_summ %>% dplyr::bind_rows(., .id = "exposure") %>%
            filter(., PP.H4.abf >= 0.7)
        write.table(sign_coloc, here(project_dir, "results", "coloc", "summary_results", str_c("significant_PPH4_0.7_", pheno_outcome, "_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        coloc_results_res %>% dplyr::bind_rows(., .id = "exposure") %>%
            write.table(., here(project_dir, "results", "coloc", "all_results", str_c(pheno_outcome, "_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        coloc_results_res %>% dplyr::bind_rows(., .id = "exposure") %>%
            filter(., exposure %in% sign_coloc$exposure) %>%
            write.table(., here(project_dir, "results", "coloc", "all_results", str_c("significant_", pheno_outcome, "_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        mr_sign_results %>%
            left_join(., coloc_results_out, by = "exposure") %>%
            dplyr::rename(coloc_nsnps = nsnps) %>%
            write.table(., here(project_dir, "MR", "results", "main", stringr::str_c("MR_significant_main_results_immune_pqtls_UKB_", pheno_outcome, ".txt")), sep = "\t", row.names = F, quote = F)

    } else { cat("There are no significant MR results. Ending without outputting anything.")}

} else { cat("There are no directionality results, therefore no analysis to perform.\n")}