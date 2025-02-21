# Description: run colocalization analysis with deCODE data for validation

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
dir_test <- try(read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("directionality_test_pqtls_decode_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
    dplyr::select(exposure, correct_causal_direction))

if (!inherits(dir_test, "try-error")) {

    # MR Egger test results:
    mregger <- try(read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_Egger_main_results_immune_deCODEvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        dplyr::select(exposure, MREgger_regression_pval = pval))

    if (inherits(mregger, "try-error")) {
        mr_sign_results <- read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_main_results_immune_deCODEvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
        filter(., method %in% c("Inverse variance weighted", "Wald ratio")) %>%
        left_join(., dir_test, by = "exposure") %>%
        filter(., pval < 0.005) %>%
        filter(., correct_causal_direction == "TRUE" | is.na(correct_causal_direction)) #%>%
        #filter(., heterogeneity_signal == "NO" | is.na(heterogeneity_signal))

    } else {

        # MR results for pheno_outcome (IVW), keeping:
            # 1) pvalue < 0.005
            # 2) MR Egger intercept not significant (i.e., no pleiotropy)
            # 4) Remove signals with heterogeneneity
        mr_sign_results <- read.table(here(project_dir, "MR", "results", "validation", stringr::str_c("MR_main_results_immune_deCODEvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", header = T) %>%
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
# Read downloaded files and keep cis-pQTLs
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if (nrow(mr_sign_results) > 0) {

        eaf_decode <- fread(here(project_dir, "gwas_sumstats", "pQTLs_deCODE", "assocvariants.annotated.txt.gz")) %>%
            dplyr::select(Name, effectAllele, otherAllele)
        cat("Loaded deCODE EAFs.\n")

        pqtl_files <- list.files(here(project_dir, "gwas_sumstats", "pQTLs_deCODE"), pattern = "*.txt.gz", full.names = T) %>%
            stringr::str_subset(., "/home/fridald4/projects/def-gsarah/fridald4/als_immune/gwas_sumstats/pQTLs_deCODE/assocvariants.annotated.txt.gz", negate = T)

        genes_to_test <- as.array(mr_sign_results$exposure)

        # Load gtf file
        cat("Loading GTF file...\n")
        gtf_file <- read.table(here(project_dir, "reference_data", "Homo_sapiens.GRCh38.v46.gtf"),
                            sep = "\t", header = F, skip = 5) %>% .[,c(1:5,9)]
        colnames(gtf_file) = c("chr","source","type","start","end","info") 

        gtf_file <- gtf_file %>%
        mutate(ensembl_id = stringr::str_extract(info, "ENS[:alpha:][:alnum:]+"),
            gene_name = stringr::str_extract(info, "gene_name [:graph:]+") %>%
            stringr::str_remove_all(., "gene_name ") %>%
            stringr::str_remove_all(., ";") %>%
            stringr::str_remove_all(., '"')
        ) %>%
        dplyr::select(.,
                source, type, chr, start, end, ensembl_id, gene_name) %>%
        filter(., type == "gene") %>%
        filter(., gene_name %in% genes_to_test)
        cat("Done.\n")
        
        for (i in 1:nrow(mr_sign_results)) {

            gene = mr_sign_results[i, "exposure"]

            gtf_file_subset <- gtf_file %>% filter(., gene_name == gene)
            chr = as.numeric(gtf_file_subset[1, "chr"] %>% stringr::str_remove(., "chr"))
            start = gtf_file_subset[1, "start"]
            end = gtf_file_subset[1, "end"]

            gene_fullname = pqtl_files %>% stringr::str_subset(., gene)

            cat("The gene being assessed is ", gene, " and the pQTLs file name is ", basename(gene_fullname), ".\n")

            exposure_gene_data <- fread(gene_fullname[1]) %>%
                mutate(Chrom = stringr::str_remove(Chrom, "chr"), varbeta = SE^2) %>%
                filter(., !is.na(Pval)) %>%
                filter(., Chrom == chr & Pos >= start-500000 & Pos <= end+500000) %>%
                filter(., !is.na(ImpMAF) | !is.na(rsids)) %>%
                filter(., ImpMAF > 0) %>%
                distinct(., rsids, .keep_all = TRUE) %>%
                dplyr::select(-effectAllele, -otherAllele) %>%
                right_join(., eaf_decode, by = "Name") %>%
                dplyr::select(SNP = rsids,
                              beta = Beta,
                              varbeta,
                              pvalues = Pval,
                              MAF = ImpMAF,
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
        write.table(coloc_summary_out, here(project_dir, "results", "coloc", "summary_results", str_c(pheno_outcome, "_decode_validation_immune_exposures.tsv")), sep = "\t", row.names = F, quote = F)

        mr_sign_results %>%
            left_join(., coloc_summary_out, by = "exposure") %>%
            dplyr::rename(coloc_nsnps = nsnps) %>%
            write.table(., here(project_dir, "MR", "results", "validation", stringr::str_c("MR_significant_main_results_immune_deCODEvalidation_pqtls_", pheno_outcome, ".txt")), sep = "\t", row.names = F, quote = F)

    } else { cat("There are no significant MR results. Ending without outputting anything.")}

} else { cat("There are no directionality results, therefore no analysis to perform.\n")}