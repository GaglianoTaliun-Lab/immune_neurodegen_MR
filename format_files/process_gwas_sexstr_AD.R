# process GWAS summary statistics for LAVA (sex-stratified analyses) - to run as a normal job on compute canada

# AD
# effect allele = ALLELE1
# coordinates are in GRCh38, so lifted them down to GRCh37 to be aligned with the other traits.

# Summary statistics must have the following columns:
# SNP / ID / SNPID_UKB/ SNPID / MarkerName / RSID / RSID_UKB: SNP IDs
# A1 / ALT: effect allele
# A2 / REF: reference allele
# N / NMISS / OBS_CT / N_analyzed: number of samples
# Z / T / STAT / Zscore: if provided, no p-values or coefficients are needed; otherwise, please provide both:
# B / BETA / OR / logOdds: effect size coefficients
# P: p-values
# summary statistics should be in GRCh37

library(here)
library(stringr)
library(devtools)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(colochelpR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
lava_ext = ".lava.gz"
dbsnp_155 <- SNPlocs.Hsapiens.dbSNP155.GRCh37

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD males no UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_FinnGen_Males_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
ADnp_m <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(ADnp_m, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("ADnp_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD females no UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_FinnGen_Females_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
ADnp_f <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(ADnp_f, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("ADnp_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD meta no UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_FinnGen_sex_het_case_control.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
ADnp_meta_sexcombined <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(ADnp_meta_sexcombined, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("ADnp_meta_sexcombined", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD males with UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_UKB_FinnGen_Males_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
AD_m <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(AD_m, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("AD_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD females with UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_UKB_FinnGen_Females_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
AD_f <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(AD_f, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("AD_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD meta with UKB |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "ADGC_ADSP_UKB_FinnGen_sex_het_case_control_reduced_bias.gwama.clean.gen090.exclude_APOE_region.shared_var")
 
AD_meta_sexcombined <- fread(file_path) %>%
   dplyr::mutate(effect_allele = ALLELE1,
                 other_allele = ALLELE0) %>%
   dplyr::select(-CHR, -BP) %>%
   filter(., stringr::str_detect(SNP, "rs")) %>%
   colochelpR::convert_rs_to_loc(., "SNP", dbsnp_155) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA,
     SE,
     P,
     N = N_incl,
     A1_AF = A1FREQ,
     variant_ID = posID)

fwrite(AD_meta_sexcombined, file = here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified", stringr::str_c("AD_meta_sexcombined", lava_ext)), sep = "\t")