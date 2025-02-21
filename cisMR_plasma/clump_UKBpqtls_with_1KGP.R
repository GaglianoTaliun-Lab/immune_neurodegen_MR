# UKB pQTLs clumping using the EUR 1000G LD reference.

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
i <- as.numeric(args[1]) # loop on number of proteins (i.e., exposures) to clump

exposure_files <- list.files(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "cis_pqtls"), pattern = ".tsv", full.names = F)

# confounder SNPs to remove to avoid horizontal pleiotropy:
BMI_rsids <- read.table(here(project_dir, "MR", "instruments", "BMI_rsids.txt"), sep = "\t", header = F)
smoking_rsids <- read.table(here(project_dir, "MR", "instruments", "smoking_rsids.txt"), sep = "\t", header = F)
confounder_snps <- rbind(BMI_rsids, smoking_rsids)

source("/home/fridald4/projects/def-gsarah/fridald4/als_immune/scripts/ld_clump_1KGP_fnc.R")

pheno_exposure_file <- exposure_files[i]
pheno_exposure <- pheno_exposure_file %>% stringr::str_remove(., ".tsv")
pheno_gene <- pheno_exposure_file %>% stringr::str_extract(., pattern = "[:alnum:]+_") %>% stringr::str_remove(., "_")

# extended HLA coords in GRCh38 (exposure genetic coordinates are in the same genome version)
HLA_start = 25726063
HLA_end = 33410226
HLA_genes <- read.table(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "list_HLA_prot_coding_genes.tsv"), header = T, sep = "\t")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load exposure (pQTLs) and clump
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if (pheno_gene %in% HLA_genes$prot_coding_genes){
    stop("Protein is within the MHC region. Skipping it without an output.\n")
}

cat(stringr::str_c("Formatting exposure for protein #", i, ":\n", pheno_exposure, ".\n"))

exposure_data <- read.table(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "cis_pqtls", pheno_exposure_file), sep = "\t", header = T) %>%
    filter(., CHR == 6 & BP >= HLA_start & BP <= HLA_end)

if (nrow(exposure_data) > 0){

    chr = exposure_data$CHR[1]

    cat("Clumping ", pheno_exposure, "(chromosome ", chr, ").\n")
        
        exposure_clumped <- exposure_data %>%
            format_data(., type = "exposure") %>%
            rename(rsid = SNP,
                pval = pval.exposure,
                trait_id = exposure,
                ) %>%
            ld_clump_1KGP(., plink_bin = "/home/fridald4/projects/def-gsarah/fridald4/als_immune/tools/plink2", 
                            bfile = stringr::str_c("/home/fridald4/projects/def-gsarah/fridald4/als_immune/reference_data/TwoSampleMR/EUR_hg38_chr",chr)) %>%
            filter(., !(rsid %in% confounder_snps$V1)) %>%
            rename(SNP = rsid,
                pval.exposure = pval,
                exposure = trait_id)

    if (nrow(exposure_clumped >= 1)){
        write.table(exposure_clumped, here(project_dir, "MR", "exposures", "ukb_pqtls", stringr::str_c("clumped_exposure_UKB_",pheno_exposure,".txt")), sep = "\t", row.names = F, quote = F)
    }
}