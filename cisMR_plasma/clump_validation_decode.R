# Exposure clumping using the 1KGP EUR LD reference.

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

source("/home/fridald4/projects/def-gsarah/fridald4/als_immune/scripts/ld_clump_1KGP_fnc.R")

# confounder SNPs to remove to avoid horizontal pleiotropy:
BMI_rsids <- read.table(here(project_dir, "MR", "instruments", "BMI_rsids.txt"), sep = "\t", header = F)
smoking_rsids <- read.table(here(project_dir, "MR", "instruments", "smoking_rsids.txt"), sep = "\t", header = F)
confounder_snps <- rbind(BMI_rsids, smoking_rsids)

# extended HLA coords in GRCh38 (exposure genetic coordinates are in the same genome version) (Reference: https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh38.p13)
HLA_start = 28510120
HLA_end = 33480577

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load exposures and clump
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#pqtls_list <- list.files(here(project_dir, "gwas_sumstats", "pQTLs_deCODE", "cis_pqtls"), pattern = "*.tsv", full.names = TRUE)
pqtls_list <- list.files(here(project_dir, "gwas_sumstats", "pQTLs_deCODE", "cis_pqtls"), pattern = "deCODE_pqtls_CTSB*", full.names = TRUE)
exposure_data_list <- lapply(pqtls_list, function(x) fread(x) %>%
    dplyr::select(
        CHR,
        BP,
        SNP,
        beta,
        se,
        effect_allele,
        other_allele,
        eaf,
        pval,
        Phenotype
    ))

cat("There are", length(exposure_data_list), "pQTL files.\n")

for (i in 1:length(exposure_data_list)){
    exposure_data <- as.data.frame(exposure_data_list[[i]])
    chr=exposure_data$CHR[1]

    if (chr == 6){
        snps_HLA <- exposure_data %>%
            filter(., BP >= HLA_start & BP <= HLA_end) %>%
            dplyr::select(SNP)
    }

    cat("Clumping chromosome ", chr, ".\n")
        
    exposure_out <- exposure_data %>%
        format_data(.) %>%
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
        
    if (nrow(exposure_out >= 1)) {
        if (chr != 6){
            write.table(exposure_out, here(project_dir, "MR", "exposures", "decode_pqtls", stringr::str_c("clumped_exposure_gene_", exposure_data$Phenotype[1],".txt")), sep = "\t", row.names = F, quote = F)
        } else {
            exposure_out %>% filter(., !(SNP %in% snps_HLA)) %>%
                write.table(., here(project_dir, "MR", "exposures", "decode_pqtls", stringr::str_c("clumped_exposure_gene_", exposure_data$Phenotype[1],".txt")), sep = "\t", row.names = F, quote = F)
        }
    }
}