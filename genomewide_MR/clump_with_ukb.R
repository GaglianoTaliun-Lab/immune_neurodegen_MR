# Exposure clumping using the UKB LD reference.

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
i <- as.numeric(args[1])  # arrays 1-9 (unique number of exposures)

mr_test_pairs <- read.table(here(project_dir, "MR", "MR_pairs_to_test.txt"), sep = "\t", header = T) %>%
  dplyr::select(exposure, type) %>%
  distinct()

palindromes = c("AT", "TA", "GC", "CG")

# extended HLA coords in GRCh37 (exposure genetic coordinates are in the same genome version)
HLA_start = 25726291
HLA_end = 33377673

# confounder SNPs to remove to avoid horizontal pleiotropy:
BMI_rsids <- read.table(here(project_dir, "MR", "instruments", "BMI_rsids.txt"), sep = "\t", header = F)
smoking_rsids <- read.table(here(project_dir, "MR", "instruments", "smoking_rsids.txt"), sep = "\t", header = F)
confounder_snps <- rbind(BMI_rsids, smoking_rsids)

exposure_list <- list()

source("/home/fridald4/projects/def-gsarah/fridald4/als_immune/scripts/ld_clump_ukb_fnc.R")

pheno_exposure <- mr_test_pairs[i,1]
type <- mr_test_pairs[i,2] # either sex-combined or sex-specific

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load exposures and clump
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat(stringr::str_c("Formatting exposure for trait in line #", i,":\n",pheno_exposure,"\n"))

if (type == "sex_combined"){

    exposure_data <- read.table(here(project_dir, "MR", "meta_analysis", "metal_output", stringr::str_c(pheno_exposure,"_wCHR_BP.tbl")), sep = "\t", header = T) %>%
        filter(pvalue <= 5e-09) %>%
    filter(CHR %in% c(1:22)) %>%
    filter(., !(CHR == 6 & BP >= HLA_start & BP <= HLA_end)) %>%
        mutate(
            Phenotype = pheno_exposure,
            effect_allele = stringr::str_to_upper(Allele1), 
            other_allele = stringr::str_to_upper(Allele2),
            alleles = stringr::str_c(effect_allele,other_allele)
        ) %>%
        subset(., !(alleles %in% palindromes)) %>%
        dplyr::select(
            CHR,
            BP,
            SNP = MarkerName,
            beta = Effect,
            se = StdErr,
            effect_allele,
            other_allele,
            pval = pvalue,
            Phenotype
        )

    unique_chrs <- unique(exposure_data$CHR)

    for (chr in unique_chrs){

	cat("Clumping chromosome ", chr, ".\n")
        
        exposure_list[[chr]] <- exposure_data %>%
            filter(., CHR == chr) %>%
            format_data(.) %>%
            rename(rsid = SNP,
                pval = pval.exposure,
                trait_id = exposure,
                ) %>%
            ld_clump_ukb(., plink_bin = "/home/fridald4/projects/def-gsarah/fridald4/als_immune/plink2", 
                        bgen = stringr::str_c("/lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr",chr,"_v3.bgen"),
            sample = "/home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample") %>%
            filter(., !(rsid %in% confounder_snps$V1)) %>%
            rename(SNP = rsid,
                pval.exposure = pval,
                exposure = trait_id)

    }

} else if (type == "sex_specific"){

    exposure_data <- read.table(here(project_dir, "MR", "meta_analysis", "metal_input", stringr::str_c(pheno_exposure,".txt")), sep = "\t", header = T) %>%
        filter(P <= 5e-09) %>%
    filter(CHR %in% c(1:22)) %>%
    filter(., !(CHR == 6 & BP >= HLA_start & BP <= HLA_end)) %>%
        mutate(
            Phenotype = pheno_exposure,
            alleles = stringr::str_c(A1,A2)
        ) %>%
        subset(., !(alleles %in% palindromes)) %>%
        dplyr::select(
            CHR,
            SNP,
            beta = BETA,
            se = SE,
            effect_allele = A1,
            other_allele = A2,
            pval = P,
            Phenotype
        )

    unique_chrs <- unique(exposure_data$CHR)

    for (chr in unique_chrs){

        cat("Clumping chromosome ", chr, ".\n")

        exposure_list[[chr]] <- exposure_data %>%
            filter(., CHR == chr) %>%
            format_data(.) %>%
            rename(rsid = SNP,
                pval = pval.exposure,
                trait_id = exposure,
                ) %>%
    ld_clump_ukb(., plink_bin = "/home/fridald4/projects/def-gsarah/fridald4/als_immune/plink2",
                        bgen = stringr::str_c("/lustre03/project/6003113/uk_biobank/genotypes/bgen/ukb_chr",chr,"_v3.bgen"),
                        sample = "/home/fridald4/projects/def-gsarah/ukb/imputation_basket_2014072/ukb22828_c5_b0_v3_s487239.sample") %>%
            filter(., !(rsid %in% confounder_snps$V1)) %>%
            rename(SNP = rsid,
                pval.exposure = pval,
                exposure = trait_id)

    }
}        

exposure_data <- rbindlist(exposure_list) %>%
    mutate(id.exposure = exposure_list[[1]]$id.exposure[1]) %>%
    write.table(., here(project_dir, "MR", "instruments", "clumped", stringr::str_c("clumped_exposure_noHLA_",pheno_exposure,"_",type,".txt")), sep = "\t", row.names = F, quote = F)

cat(stringr::str_c("There are ", nrow(exposure_data), " instruments for harmonizing.\n"))
