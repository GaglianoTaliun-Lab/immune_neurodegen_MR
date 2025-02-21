# Description: regional plots of colocalization results with H4 > 0.7 between neurodegenerative diseases and plasma pQTLs

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
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Arguments
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

args <- commandArgs(TRUE)
pheno_outcome <- as.character(args[1])

source("/scratch/fridald4/als_immune/ggmirror_regional.R")

# Olink map with 1,248 immune-related proteins
olink_map <- as.data.frame(fread(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "olink_map_with_immune_proteins.tsv"))) %>%
  mutate(pheno_id = stringr::str_c(HGNC.symbol , "_", OlinkID,"_",UniProt))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Read coloc significant coloc results
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Reading coloc file...\n")
coloc_file <- here(project_dir, "results", "coloc", "summary_results", stringr::str_c("significant_PPH4_0.7_", pheno_outcome, "_immune_exposures.tsv"))

coloc <- try(read.table(coloc_file, sep = "\t", header = T) %>%
  mutate(gene_name = stringr::str_extract(exposure, "[:alnum:]+_") %>% stringr::str_remove(., "_")), silent = TRUE)

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

    pheno_exposure <- as.character(coloc[i,"exposure"])
    gene_tested <- as.character(coloc[i,"gene_name"])

    cat("Gene to test:", gene_tested, ".\n")

    olink_map_sub <- olink_map %>%
      filter(., pheno_id == pheno_exposure)

    chr = olink_map_sub[1, "chr"]
    gene_start = olink_map_sub[1, "gene_start"]
    gene_end = olink_map_sub[1, "gene_end"]
    name_id = olink_map_sub[1, "UKBPPP_ProteinID"]
    Panel = olink_map_sub[1, "Panel"]
    directory_name = olink_map_sub[1, "file_name"] %>% stringr::str_remove(., ".tar")
    file_name_exp = stringr::str_c("discovery_chr", chr, "_", name_id, ":", Panel, ".gz")

    rsid_map <- fread(stringr::str_c("/scratch/fridald4/als_immune/ukb_pqtls/olink_rsid_map_mac5_info03_b0_7_chr", chr, "_patched_v2.tsv.gz")) %>%
      dplyr::select(rsid, GENPOS_hg38 = POS38, GENPOS_hg19 = POS19)

    exposure_gene_data <- fread(stringr::str_c("/scratch/fridald4/als_immune/ukb_pqtls/", directory_name, "/", file_name_exp)) %>%
      mutate(P = case_when(
        LOG10P <= -log10(1e-200) ~ as.numeric(10^-LOG10P),
        LOG10P > -log10(1e-200) ~ as.numeric(1e-200))) %>%
      filter(CHROM == chr & GENPOS >= gene_start - 500000 & GENPOS <= gene_end + 500000) %>%
      left_join(., rsid_map, by = c("GENPOS" = "GENPOS_hg38")) %>%
      dplyr::select(SNP = rsid, CHR = CHROM, BP = GENPOS_hg19, P)

    topSNP_exposure <- dplyr::arrange(exposure_gene_data, P) %>% .[1,"SNP"]

    outcome_data_sub <- outcome_data %>%
      filter(., SNP %in% exposure_gene_data$SNP)

    #write.table(exposure_gene_data, here(project_dir, "results", "coloc", "regional_plots", stringr::str_c("regional_plot_exposure_data_",pheno_outcome, "_", pheno_exposure,".txt")), sep = "\t", quote = F, row.names = F)
    #write.table(outcome_data_sub, here(project_dir, "results", "coloc", "regional_plots", stringr::str_c("regional_plot_outcome_data_",pheno_outcome, "_", pheno_exposure,".txt")), sep = "\t", quote = F, row.names = F)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Regional plots
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    p <- gmirror_v2(top=exposure_gene_data, bottom=outcome_data_sub, tline=5e-08, bline=5e-08,
      toptitle=gene_tested, bottomtitle=pheno_outcome, opacity= 0.8, log10 = TRUE, annotate_snp = topSNP_exposure,
      highlight_snp = topSNP_exposure, highlight_p = c(5e-08,5e-08), highlighter_snp="yellow", highlighter_p="#999999", freey = TRUE)

    cat("Saving plot for ", pheno_outcome, " and ", pheno_exposure, "...\n")  
    ggsave(filename = here(project_dir, "results", "coloc", "regional_plots", stringr::str_c("regional_plot_",pheno_outcome, "_", pheno_exposure,".jpg")),
      plot = p, units = "cm", height = 30, width = 30)
    cat("Done.\n")

  }
}