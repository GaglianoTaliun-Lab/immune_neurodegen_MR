# process GWAS summary statistics for LAVA (sex-stratified analyses) - to run as a normal job on compute canada

# UKB autoimmune (quantitative)

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
library(BSgenome.Hsapiens.UCSC.hg19)
library(colochelpR)
library(rutils)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
lava_ext = ".lava.gz"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37
source(here(project_dir, "scripts", "get_other_allele.R"))
hg19 <- BSgenome.Hsapiens.UCSC.hg19

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Neutrophil percentage (code 30200) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# males and females in same file:
# file_path <- file.path(project_dir, "gwas_sumstats", "autoimmune_sex_stratified", "beta.sexcomparison_30200-0.0.tsv.gz")
# fread(file_path, nrows = 6)

# neu <- fread(file_path) 

# # males
# neu_m <- neu %>%
#    dplyr::mutate(N = 200368, A2 = NA) %>%
#    filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
#    dplyr::select(
#      SNP,
#      A1 = ALLELE,
#      A2,
#      BETA = `NBETA_M-30200-0.0`,
#      SE = `NSE_M-30200-0.0`,
#      P = `BPV_M-30200-0.0`,
#      N) %>%
#    colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
#    tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
#    get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
#    dplyr::mutate(CHR = case_when(
#                    CHR == "X" ~ as.numeric(23),
#                    TRUE ~ as.numeric(CHR)),
#                  variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
#    dplyr::select(
#      SNP,
#      CHR,
#      BP,
#      A1,
#      A2,
#      BETA,
#      SE,
#      P,
#      N,
#      variant_ID)

# fwrite(neu_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("neutrophper_males", lava_ext)), sep = "\t")

# females
# neu_f <- neu %>%
#    dplyr::mutate(N = 237893, A2 = NA) %>%
#    filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
#    dplyr::select(
#      SNP,
#      A1 = ALLELE,
#      A2,
#      BETA = `NBETA_F-30200-0.0`,
#      SE = `NSE_F-30200-0.0`,
#      P = `BPV_F-30200-0.0`,
#      N) %>%
#    colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
#    tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
#    get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
#    dplyr::mutate(CHR = case_when(                         
#                    CHR == "X" ~ as.numeric(23),
#                    TRUE ~ as.numeric(CHR)),
#                  variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
#    dplyr::select(
#      SNP,
#      CHR,
#      BP, 
#      A1,
#      A2,
#      BETA,
#      SE,
#      P,
#      N, 
#      variant_ID)

# fwrite(neu_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("neutrophper_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Eosinophil percentage (code 30210) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# males and females in same file:
file_path <- file.path(project_dir, "gwas_sumstats", "autoimmune_sex_stratified", "beta.sexcomparison_30210-0.0.tsv.gz")
fread(file_path, nrows = 6)

eos <- fread(file_path) 

# males
eos_m <- eos %>%
   dplyr::mutate(N = 200277, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_M-30210-0.0`,
     SE = `NSE_M-30210-0.0`,
     P = `BPV_M-30210-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1,
     A2,
     BETA,
     SE,
     P,
     N,
     variant_ID)

fwrite(eos_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("eosinophper_males", lava_ext)), sep = "\t")

# females
eos_f <- eos %>%
   dplyr::mutate(N = 237784, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_F-30210-0.0`,
     SE = `NSE_F-30210-0.0`,
     P = `BPV_F-30210-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(                         
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP, 
     A1,
     A2,
     BETA,
     SE,
     P,
     N, 
     variant_ID)

fwrite(eos_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("eosinophper_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Monocyte percentage (code 30190) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# males and females in same file:
file_path <- file.path(project_dir, "gwas_sumstats", "autoimmune_sex_stratified", "beta.sexcomparison_30190-0.0.tsv.gz")
fread(file_path, nrows = 6)

mon <- fread(file_path) 

# males
mon_m <- mon %>%
   dplyr::mutate(N = 200277, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_M-30190-0.0`,
     SE = `NSE_M-30190-0.0`,
     P = `BPV_M-30190-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1,
     A2,
     BETA,
     SE,
     P,
     N,
     variant_ID)

fwrite(mon_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("mono_males", lava_ext)), sep = "\t")

# females
mon_f <- mon %>%
   dplyr::mutate(N = 237784, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_F-30190-0.0`,
     SE = `NSE_F-30190-0.0`,
     P = `BPV_F-30190-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(                         
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP, 
     A1,
     A2,
     BETA,
     SE,
     P,
     N, 
     variant_ID)

fwrite(mon_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("mono_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Basophil percentage (code 30220) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# males and females in same file:
file_path <- file.path(project_dir, "gwas_sumstats", "autoimmune_sex_stratified", "beta.sexcomparison_30220-0.0.tsv.gz")
fread(file_path, nrows = 6)

bas <- fread(file_path) 

# males
bas_m <- bas %>%
   dplyr::mutate(N = 200277, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_M-30220-0.0`,
     SE = `NSE_M-30220-0.0`,
     P = `BPV_M-30220-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1,
     A2,
     BETA,
     SE,
     P,
     N,
     variant_ID)

fwrite(bas_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("basophil_males", lava_ext)), sep = "\t")

# females
bas_f <- bas %>%
   dplyr::mutate(N = 237784, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_F-30220-0.0`,
     SE = `NSE_F-30220-0.0`,
     P = `BPV_F-30220-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(                         
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP, 
     A1,
     A2,
     BETA,
     SE,
     P,
     N, 
     variant_ID)

fwrite(bas_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("basophil_females", lava_ext)), sep = "\t")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                               | Basophil count (code 30160) |
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# males and females in same file:
file_path <- file.path(project_dir, "gwas_sumstats", "autoimmune_sex_stratified", "beta.sexcomparison_30160-0.0.tsv.gz")
fread(file_path, nrows = 6)

bas <- fread(file_path) 

# males
bas_m <- bas %>%
   dplyr::mutate(N = 200277, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_M-30160-0.0`,
     SE = `NSE_M-30160-0.0`,
     P = `BPV_M-30160-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP,
     A1,
     A2,
     BETA,
     SE,
     P,
     N,
     variant_ID)

fwrite(bas_m, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("basophilper_males", lava_ext)), sep = "\t")

# females
bas_f <- bas %>%
   dplyr::mutate(N = 237784, A2 = NA) %>%
   filter(SNP %in% stringr::str_subset(SNP, "^rs")) %>%
   dplyr::select(
     SNP,
     A1 = ALLELE,
     A2,
     BETA = `NBETA_F-30160-0.0`,
     SE = `NSE_F-30160-0.0`,
     P = `BPV_F-30160-0.0`,
     N) %>%
   colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
   tidyr::separate(., loc, c("CHR", "BP"), sep = ":") %>%
   get_other_allele(df = ., SNP_column = "SNP", dbSNP = dbsnp_144, genome = hg19) %>%
   dplyr::mutate(CHR = case_when(                         
                   CHR == "X" ~ as.numeric(23),
                   TRUE ~ as.numeric(CHR)),
                 variant_ID = paste(CHR,BP,A1,A2,sep=":")) %>%
   dplyr::select(
     SNP,
     CHR,
     BP, 
     A1,
     A2,
     BETA,
     SE,
     P,
     N, 
     variant_ID)

fwrite(bas_f, file = here(project_dir, "gwas_sumstats","lava_sumstats","gwas_sexstratified",stringr::str_c("basophilper_females", lava_ext)), sep = "\t")
