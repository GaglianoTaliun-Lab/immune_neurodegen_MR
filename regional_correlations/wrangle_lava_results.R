# Description: wrangle lava results

# Load packages -----------------------------------------------------------

library(dplyr)
library(here)
library(stringr)
library(tidyr)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

args <- commandArgs(TRUE)

# use LAVA results either with PD with proxies (1) or PD without proxies (2):
arg_PD_proxies <- as.numeric(args[1])

# name of LAVA ouput:
if (arg_PD_proxies == 1) {

  datasets <- "ALS_males:ALS_females:PD_males:PD_females:AD_males:AD_females"
  datasets_chrx <- "ALSX_males:ALSX_females:PDX_males:PDX_females:ADX_males:ADX_females"
  loci_file <- "blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile"
  out_prefix <- "neurodegen"

} else if (arg_PD_proxies == 2) {

  datasets <- "ALS_males:ALS_females:PDnp_males:PDnp_females:ADnp_males:ADnp_females.autosomes"
  loci_file <- "blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile"
  out_prefix <- "neurodegen_noUKB"

} else {
  
  stop("Error: need to provide argument 1 (i.e., include PD dataset with proxies) or 2 (i.e., include PD dataset without proxies).")
  
}

# Load data ---------------------------------------------------------------

# read RDS files

if (arg_PD_proxies == 1) {

  bivarA <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".A.bivar.lava.rds")))
  bivarB <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".B.bivar.lava.rds")))

  univarA <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".A.univ.lava.rds")))
  univarB <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".B.univ.lava.rds")))

  bivar_x <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets_chrx,".chrx.bivar.lava.rds")))
  univar_x <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets_chrx,".chrx.univ.lava.rds")))
  
  bivar <- c(bivarA, bivarB, bivar_x)
  univar <- c(univarA, univarB, univar_x)

} else if (arg_PD_proxies == 2) {
  bivarA <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".bivar.lava.A.rds")))
  bivarB <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".bivar.lava.B.rds")))
  bivarC <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".bivar.lava.C.rds")))

  univarA <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".univ.lava.A.rds")))
  univarB <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".univ.lava.B.rds")))
  univarC <- readRDS(here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(datasets,".univ.lava.C.rds")))

  bivar <- c(bivarA, bivarB, bivarC)
  univar <- c(univarA, univarB, univarC)
}

###---------------------------------------- Main

# Univariate -----------------------------------

# get number of univariate tests (number of loci in LAVA locus input)
ntest_univ <- read.table(here(project_dir, "reference_data", loci_file), sep = "\t", header = T) %>% 
  nrow(.)

# create dataframe with all univariate results
univ_all <- bind_rows(univar)

# create dataframe with only significant univariate results
univ_sign <- univ_all %>%
  filter(., p < 0.05/ntest_univ)

# Bivariate ------------------------------------

# get number of bivariate tests
# ntests = 0
# for (i in 1:length(bivar)) {
#   condition = nrow(bivar[[i]]) > 0
#   if (condition == "TRUE" && length(condition) != 0) {
#     ntests <- ntests + nrow(bivar[[i]])
#   }
# }

# remove null loci
bivar <- bivar[!sapply(bivar,is.null)]
bivar_all <- bind_rows(bivar)

# keep only cross trait loci by matching the sex
bivar_all <- bivar_all %>%
  separate(phen1, c("trait1", "sex1"), remove = FALSE, sep = "_") %>%
  separate(phen2, c("trait2", "sex2"), remove = FALSE, sep = "_") %>%
  filter(., sex1 == sex2) %>%
  dplyr::select(-trait1, -trait2, -sex1, -sex2)

ntests = nrow(bivar_all)

# get pvalue bonferroni threshold
pvalue_bivar = 0.05/ntests

# filter out non significant bivariate tests
bivar_significant <- bivar_all %>% filter(., p <= pvalue_bivar)
# bivar_significant <- lapply(bivar, function(x) filter(x, p <= pvalue_bivar))
# bivar_significant <- bind_rows(bivar_significant)

# obtain number of significant bivariate tests
n_sig = nrow(bivar_significant)

### ---------------------------------------------- Save files

# write all results into table
write.table(bivar_all, here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(out_prefix, "_sexstr_bivar_all.tsv")), sep = "\t", quote = F, row.names = F)
write.table(univ_all, here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(out_prefix, "_sexstr_univ_all.tsv")), sep = "\t", quote = F, row.names = F)

# write significant results into table
write.table(bivar_significant, here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(out_prefix, "_sexstr_bivar_significant.tsv")), sep = "\t", quote = F, row.names = F)
write.table(univ_sign, here(project_dir,"results","lava","gwas_sexstr",stringr::str_c(out_prefix, "_sexstr_univ_significant.tsv")), sep = "\t", quote = F, row.names = F)


