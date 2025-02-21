# Exposure clumping using the EUR 1KGP as reference.

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
sign_genes <- read.table(here(project_dir, "gwas_sumstats", "pQTLs_Zhang2022", "list_sign_genes.txt"), header = F, sep = "\t")

# extended HLA coords in GRCh37 (exposure genetic coordinates are in the same genome version)
# HLA gene list from: https://hla.alleles.org/genes/index.html
HLA_genes <- read.table(here(project_dir, "MR","HLA_genes.txt"), sep = "\t", header = T)
HLA_start = 25726291
HLA_end = 33377673

palindromes = c("AT", "TA", "GC", "CG")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load exposures and clump
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:nrow(sign_genes)) {

    gene = sign_genes[i,1]
    if (gene %in% HLA_genes$name) next

    pqtls <- read.table(here(project_dir, "gwas_sumstats", "pQTLs_Zhang2022", "pqtls_for_MR", stringr::str_c(gene, ".tsv")), header = T, sep = "\t")
    chr = as.numeric(pqtls[1, "CHR"])

    if (chr == 6){
        snps_HLA <- pqtls %>%
            filter(., BP >= HLA_start & BP <= HLA_end) %>%
            dplyr::select(SNP)
    }

    cat("Clumping for gene ", gene, ".\n")
    # pqtls <- pqtls %>%
    #    mutate(alleles = stringr::str_c(effect_allele,other_allele)) %>%
    #    subset(., !(alleles %in% palindromes))

    # case where all instruments left were plaindromes, then skip - but I commented the lines above, so all should have nrow >= 1.
    if (nrow(pqtls >= 1)) {

        pqtls_clumped <- format_data(pqtls, type = "exposure") %>%
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

        if (nrow(pqtls_clumped >= 1)) {
            if (chr != 6){
                write.table(pqtls_clumped, here(project_dir, "MR", "exposures", "aric_pqtls", stringr::str_c("clumped_exposure_gene_",gene,".txt")), sep = "\t", row.names = F, quote = F)
            } else {
                pqtls_clumped %>% filter(., !(SNP %in% snps_HLA)) %>%
                    write.table(., here(project_dir, "MR", "exposures", "aric_pqtls", stringr::str_c("clumped_exposure_gene_",gene,".txt")), sep = "\t", row.names = F, quote = F)
    
            }
        }
    }
}

