# Description: regional plots of colocalization results with H4 > 0.7 between neurodegenerative diseases and CSF pQTLs

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load packages
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(vctrs)
library(BSgenome)
library(colochelpR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
dbsnp_155 <- SNPlocs.Hsapiens.dbSNP155.GRCh37

args <- commandArgs(TRUE)
pheno_outcome <- as.character(args[1])

source("/scratch/fridald4/als_immune/ggmirror_regional.R")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read coloc significant coloc results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Reading coloc file...\n")
coloc_file <- here(project_dir, "results", "coloc", "summary_results", stringr::str_c("significant_PPH4_0.7_", pheno_outcome, "_CSF_immune_exposures.tsv"))

coloc <- try(read.table(coloc_file, sep = "\t", header = T), silent = TRUE)

if (class(coloc) == "try-error"){
  stop("There is no coloc file for ", pheno_outcome, ". Exiting without any output.\n")
}

cat("Done.\n")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read genome-wide outcome data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Reading outcome data...\n")
if ((stringr::str_detect(pheno_outcome, "meta_sexcombined") == TRUE)){
  file_name = here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(pheno_outcome,"_wCHR_BP.tbl"))
} else if ((stringr::str_detect(pheno_outcome, "males") == TRUE)){
  file_name = here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(pheno_outcome,".txt"))
}

outcome_data <- read.table(file_name, sep = "\t", header = T) %>%
  dplyr::select(SNP, CHR, BP, P) %>%
  mutate(P = case_when(
    P >= 1e-200 ~ as.numeric(P),
    P < 1e-200 ~ as.numeric(1e-200)))
cat("Done.\n")

# loop across significant coloc results - outcome is the same only exposure varies:

if (nrow(coloc) > 0) {

  cat("Start looping across exposures for ", pheno_outcome, "...\n")
  for (i in 1:nrow(coloc)) { 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read exposure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    exposure_files <- setNames(list.files(here(project_dir, "gwas_sumstats","pQTLs_Yang2021","cis_pqtls_for_coloc"), pattern = "*.tsv", full.names = TRUE),
                                nm = list.files(here(project_dir, "gwas_sumstats","pQTLs_Yang2021","cis_pqtls_for_coloc"), pattern = "*.tsv", full.names = FALSE) %>%
                                stringr::str_remove(., ".tsv"))

    pheno_exposure = as.character(coloc[i,"exposure"])

    exposure_gene_data <- fread(exposure_files[[pheno_exposure]]) %>%
      dplyr::select(-CHR, -BP) %>%
      filter(., stringr::str_detect(SNP, "rs")) %>%
      colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
      tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
      mutate(P = case_when(
        pvalues >= 1e-200 ~ as.numeric(pvalues),
        pvalues < 1e-200 ~ as.numeric(1e-200))) %>%
      dplyr::select(SNP, CHR, BP, P)

    topSNP_exposure <- dplyr::arrange(exposure_gene_data, P) %>% .[1,"SNP"]

    outcome_data_sub <- outcome_data %>%
      filter(., SNP %in% exposure_gene_data$SNP)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regional plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    p <- gmirror_v2(top=exposure_gene_data, bottom=outcome_data_sub, tline=5e-08, bline=5e-08,
      toptitle=pheno_exposure, bottomtitle=pheno_outcome, opacity= 0.8, log10 = TRUE, annotate_snp = topSNP_exposure,
      highlight_snp = topSNP_exposure, highlight_p = c(5e-08,5e-08), highlighter_snp="yellow", highlighter_p="#999999", freey = TRUE)

    cat("Saving plot for ", pheno_outcome, " and ", pheno_exposure, "...\n")  
    ggsave(filename = here(project_dir, "results", "coloc", "regional_plots", stringr::str_c("regional_plot_CSF_",pheno_outcome, "_", pheno_exposure,".jpg")),
      plot = p, units = "cm", height = 30, width = 30)
    cat("Done.\n")

  }
}