# process chromosome X GWAS summary statistics for LAVA (sex-stratified analyses)

# ALS (GRCh37), PD (GRCh37), AD (GRCh38)
# output does not include CHR and BP (not needed for LAVA).

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
library(biomaRt)
library(BSgenome)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(colochelpR)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

home_dir="/home/fridald4/projects/def-gsarah"
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
lava_ext = ".lava.gz"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# Read ALS X chr file (includes info for both males and females)
file_path_Xchr <- file.path(home_dir, "als_sumstats_Byrne2023/xchr_Byrne_2023_strata_geno0.05_maf0.01_20pcs_fishers_info_0.3.01.xstrat.logistic.gz")
ALS_X <- fread(file_path_Xchr) 

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | ALS males |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ALS_m_X <- ALS_X %>%
  dplyr::mutate(CHR = "X",
                N = 21991) %>%
  dplyr::select(
    CHR,
    BP,
    A1,
    A2,
    BETA = BETA_M,
    SE = SE_M,
    P = P_M,
    N) %>%
    colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
    mutate(CHR = 23)

# then need to remove positions with more than one rsid:
ALS_m_X <- ALS_m_X %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP) %>%
 dplyr::select(
 SNP,
 A1,
 A2,
 BETA,
 SE,
 P,
 N) %>% as.data.frame()

fwrite(ALS_m_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ALSX_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | ALS females |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ALS_f_X <- ALS_X %>%
  dplyr::mutate(CHR = "X",
                N = 19576) %>%
  dplyr::select(
    CHR,
    BP,
    A1,
    A2,
    BETA = BETA_F,
    SE = SE_F,
    P = P_F,
    N) %>%
    colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
    mutate(CHR = 23)

# then need to remove positions with more than one rsid:
ALS_f_X <- ALS_f_X %>%
 dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP)) %>%
 dplyr::group_by(CHR_BP) %>%
 dplyr::filter(!any(row_number() > 1))  %>%
 dplyr::ungroup() %>%
 dplyr::select(-CHR_BP) %>%
 dplyr::select(
 SNP,
 A1,
 A2,
 BETA,
 SE, 
 P,
 N) %>% as.data.frame()
 
fwrite(ALS_f_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ALSX_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | PD males w/ proxies |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified", "GCST90104086_buildGRCh37.tsv")

PD_m_X <- fread(file_path) %>%
  mutate(BETA = log(male_odds_ratio), N = 261231) %>%
  dplyr::select(
  SNP = variant_id,
  A1 = effect_allele,
  A2 = other_allele,
  BETA,
  SE = male_odds_ratio_se, 
  P = `male_p-val`,
  N)

fwrite(PD_m_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PDX_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | PD females w/ proxies |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified", "GCST90104085_buildGRCh37.tsv")

PD_f_X <- fread(file_path) %>%
  mutate(BETA = log(female_odds_ratio), N = 302908) %>%
  dplyr::select(
  SNP = variant_id,
  A1 = effect_allele,
  A2 = other_allele,
  BETA,
  SE = female_odds_ratio_se, 
  P = `female_p-val`,
  N)

fwrite(PD_f_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PDX_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD males w/ proxies |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "GCST90444376.tsv")

AD_m_X <- fread(file_path) %>%
  dplyr::select(
    SNP = rs_id,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error, 
    P = p_value,
    N = n)

  fwrite(AD_m_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ADX_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | AD females w/ proxies |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "AD_sexstratified", "GCST90444375.tsv")

AD_f_X <- fread(file_path) %>%
  dplyr::select(
    SNP = rs_id,
    A1 = effect_allele,
    A2 = other_allele,
    BETA = beta,
    SE = standard_error, 
    P = p_value,
    N = n)

  fwrite(AD_f_X, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ADX_females", lava_ext)), sep = "\t")