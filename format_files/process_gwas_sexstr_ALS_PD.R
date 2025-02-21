# process GWAS summary statistics for LAVA (sex-stratified analyses) - to run as a normal job on compute canada

# ALS, PD

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
source(here(project_dir, "scripts", "get_other_allele.R"))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | ALS males |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# file_path <- file.path(home_dir, "als_sumstats_Byrne2023/male_ALS_sumstats_Byrne_2023_SAIGE_LOCO.gz")
# # fread(file_path, nrows = 6)
 
# ALS_m <- fread(file_path) %>%
#    dplyr::mutate(CHR = as.factor(chr),
#                  effect_allele = stringr::str_to_upper(A1),
#                  other_allele = stringr::str_to_upper(A2)) %>%
#    dplyr::select(
#      SNP = rsid,
#      CHR,
#      BP = pos,
#      A1 = effect_allele,
#      A2 = other_allele,
#      BETA = Effect,
#      SE = StdErr,
#      P = `P-value`,
#      A1_AF = Freq1,
#      N = Neff) %>%
#    dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# head(ALS_m)

# # include X chromosome (in different file)
# file_path_Xchr <- file.path(home_dir, "als_sumstats_Byrne2023/xchr_Byrne_2023_strata_geno0.05_maf0.01_20pcs_fishers_info_0.3.01.xstrat.logistic.gz")
# # fread(file_path, nrows = 6)

# ALS_X <- fread(file_path_Xchr) 

# ALS_m_X <- ALS_X %>%
#   dplyr::mutate(CHR = "X",
#                 N = 21991, A1_AF = NA) %>%
#   dplyr::select(
#     CHR,
#     BP,
#     A1,
#     A2,
#     BETA = BETA_M,
#     SE = SE_M,
#     P = P_M,
#     A1_AF,
#     N) %>%
#     colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
#     mutate(CHR = 23) %>%
#   dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# tail(ALS_m_X)

# # then need to remove positions with more than one rsid:
# ALS_m_X <-
#  ALS_m_X %>%
#  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP),
#                SNP = case_when(
# 		is.na(SNP) ~ variant_ID,
# 		TRUE ~ SNP)) %>%
#  dplyr::group_by(CHR_BP) %>%
#  dplyr::filter(!any(row_number() > 1))  %>%
#  dplyr::ungroup() %>%
#  dplyr::select(-CHR_BP) %>%
#  dplyr::select(
#  SNP,
#  CHR,
#  BP,
#  A1,
#  A2,
#  BETA,
#  SE,
#  P,
#  A1_AF,
#  N,
#  variant_ID) %>% as.data.frame()

# head(ALS_m_X)

# ALS_m_all <- rbind(ALS_m, ALS_m_X)

# fwrite(ALS_m_all, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ALS_males", lava_ext)), sep = "\t")

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #                                               | ALS females |
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# file_path <- file.path(home_dir, "als_sumstats_Byrne2023/female_ALS_sumstats_Byrne_2023_SAIGE_LOCO.gz")
# fread(file_path, nrows = 6)
 
# ALS_f <- fread(file_path) %>%
#    dplyr::mutate(CHR = as.factor(chr),
#                  effect_allele = stringr::str_to_upper(A1),
#                  other_allele = stringr::str_to_upper(A2)) %>%
#    dplyr::select(
#      SNP = rsid,
#      CHR,
#      BP = pos,
#      A1 = effect_allele,
#      A2 = other_allele,
#      BETA = Effect,
#      SE = StdErr,
#      P = `P-value`,
#      A1_AF = Freq1,
#      N = Neff) %>%
#    dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# # include X chromosome (in different file read above)
# ALS_f_X <- ALS_X %>%
#   dplyr::mutate(CHR = "X",
#                 N = 19576, A1_AF = NA) %>%
#   dplyr::select(
#     CHR,
#     BP,
#     A1,
#     A2,
#     BETA = BETA_F,
#     SE = SE_F,
#     P = P_F,
#     A1_AF,
#     N) %>%
#     colochelpR::convert_loc_to_rs(., dbsnp_144) %>%
#     mutate(CHR = 23) %>%
#   dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# # then need to remove positions with more than one rsid:
# ALS_f_X <-
#  ALS_f_X %>%
#  dplyr::mutate(CHR_BP = stringr::str_c(CHR, ":", BP),
#                SNP = case_when(
#                 is.na(SNP) ~ variant_ID,
#                 TRUE ~ SNP)) %>%
#  dplyr::group_by(CHR_BP) %>%
#  dplyr::filter(!any(row_number() > 1))  %>%
#  dplyr::ungroup() %>%
#  dplyr::select(-CHR_BP) %>%
#  dplyr::select(
#  SNP,
#  CHR,
#  BP,
#  A1,
#  A2,
#  BETA,
#  SE, 
#  P,
#  A1_AF,
#  N,
#  variant_ID) %>% as.data.frame()

# ALS_f_all <- rbind(ALS_f, ALS_f_X)
 
# fwrite(ALS_f_all, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("ALS_females", lava_ext)), sep = "\t")

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #                                               | PD males w/ proxies |
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified/MALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt.gz")
# fread(file_path, nrows = 6)

# PD_m <- fread(file_path) %>%
#    separate(., MarkerName, c("CHR", "BP"), sep = ":") %>%
#    dplyr::mutate(CHR = as.factor(CHR),
#                  effect_allele = stringr::str_to_upper(Allele1),
#                  other_allele = stringr::str_to_upper(Allele2),
#                  N = 20956 + 89660) %>%
#    dplyr::select(
#      SNP = ID,
#      CHR,
#      BP,
#      A1 = effect_allele,
#      A2 = other_allele,
#      BETA = Effect,
#      SE = StdErr,
#      P = `P-value`,
#      A1_AF = Freq1,
#      N) %>%
#    dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# fwrite(PD_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PD_males", lava_ext)), sep = "\t")

# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# #                                               | PD females w/ proxies |
# #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified/FEMALE_PD_filtered_sumstats_no_multi_allelics_RSID.txt.gz")
# fread(file_path, nrows = 6)

# PD_f <- fread(file_path) %>%
#    separate(., MarkerName, c("CHR", "BP"), sep = ":") %>%
#    dplyr::mutate(CHR = as.factor(CHR),
#                  effect_allele = stringr::str_to_upper(Allele1),
#                  other_allele = stringr::str_to_upper(Allele2),
#                  N = 13420 + 90662) %>%
#    dplyr::select(
#      SNP = ID,
#      CHR,
#      BP,
#      A1 = effect_allele,
#      A2 = other_allele,
#      BETA = Effect,
#      SE = StdErr,
#      P = `P-value`,
#      A1_AF = Freq1,
#      N) %>%
#    dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

# fwrite(PD_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PD_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | PD males w/out UKB at all |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified/MALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz")
fread(file_path, nrows = 6)

PD_m_np <- fread(file_path) %>%
   separate(., MarkerName, c("CHR", "BP"), sep = ":") %>%
   dplyr::mutate(CHR = as.factor(CHR),
                 effect_allele = stringr::str_to_upper(Allele1),
                 other_allele = stringr::str_to_upper(Allele2),
                 N = 12054 + 11999) %>%
   dplyr::select(
     SNP = ID,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA = Effect,
     SE = StdErr,
     P = `P-value`,
     A1_AF = Freq1,
     N) %>%
   dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(PD_m_np, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PDnp_males", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | PD females w/out UKB at all |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

file_path <- file.path(project_dir, "gwas_sumstats", "PD_sex_stratified/FEMALE_PD_filtered_sumstats_NO_UKB_AT_ALL_no_multi_allelics_RSID.txt.gz")
fread(file_path, nrows = 6)

PD_f_np <- fread(file_path) %>%
   separate(., MarkerName, c("CHR", "BP"), sep = ":") %>%
   dplyr::mutate(CHR = as.factor(CHR),
                 effect_allele = stringr::str_to_upper(Allele1),
                 other_allele = stringr::str_to_upper(Allele2),
                 N = 7384 + 12389) %>%
   dplyr::select(
     SNP = ID,
     CHR,
     BP,
     A1 = effect_allele,
     A2 = other_allele,
     BETA = Effect,
     SE = StdErr,
     P = `P-value`,
     A1_AF = Freq1,
     N) %>%
   dplyr::mutate(variant_ID = paste(CHR,BP,A1,A2,sep=":"))

fwrite(PD_f_np, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("PDnp_females", lava_ext)), sep = "\t")

