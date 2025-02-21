# Description:
# step 1. From the significant LAVA loci, obtain information if these loci are genome-wide significant in the respective GWAS
# step 2. obtain GWAS summary statistics of the most-significant rsids, for traits which genetic correlations pairs had p < 0.05

# Note: this is done once with PD including case proxies (and chromosome X), and a second time with PD without case proxies 
# and no chr X in both cases, since I don't have info for BP!

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Libraries
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Arguments
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
list_output = list()
row_n = 1
args <- commandArgs(TRUE)
arg_PD_proxies <- as.numeric(args[1])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  Read files
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Read files...\n")

if (arg_PD_proxies == 1) {
  sumstats_list <- setNames(
    object = list(
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ALS_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ALS_females.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "PD_females.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "PD_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "AD_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "AD_females.lava.gz"))
    ),
    nm = c("ALS_males", "ALS_females", "PD_females", "PD_males", "AD_males", "AD_females")
  )

  lava_univ_sign <- read.table(here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_sexstr_univ_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05/2589) %>%
    filter(., chr != 23)
  lava_bivar_sign <- read.table(here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_sexstr_bivar_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05) %>%
    filter(., chr != 23)

} else if (arg_PD_proxies == 2) {
  sumstats_list <- setNames(
    object = list(
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ALS_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ALS_females.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "PDnp_females.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "PDnp_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ADnp_males.lava.gz")),
      fread(here(project_dir,"gwas_sumstats", "lava_sumstats", "gwas_sexstratified", "ADnp_females.lava.gz"))
    ),
    nm = c("ALS_males", "ALS_females", "PDnp_females", "PDnp_males", "ADnp_males", "ADnp_females")
  )

  lava_univ_sign <- read.table(here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_noUKB_sexstr_univ_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05/2589) %>%
    filter(., chr != 23)
  lava_bivar_sign <- read.table(here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_noUKB_sexstr_bivar_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05) %>%
    filter(., chr != 23)

}

cat("Done.\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main - step 1.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Starting step 1 for univariate signals...\n")

# for univariate signals:
for (i in 1:nrow(lava_univ_sign)){
  chr=lava_univ_sign[i,2]
  start=lava_univ_sign[i,3]
  stop=lava_univ_sign[i,4]
  trait=lava_univ_sign[i,7]

  tmp_trait <- sumstats_list[[trait]] %>%
    filter(., CHR == chr & BP >= start & BP <= stop)

  min_rsid_trait <- arrange(tmp_trait, P) %>% .[1, "SNP"]
  min_pval_trait <- arrange(tmp_trait, P) %>% .[1, "P"]
  lava_univ_sign[i, "most_sign_SNP_Trait"] <- min_rsid_trait
  lava_univ_sign[i, "pval_most_sign_SNP_Trait"] <- min_pval_trait
}

cat("Starting step 1 for bivariate signals...\n")

# for bivariate signals:
for (i in 1:nrow(lava_bivar_sign)){
  chr=lava_bivar_sign[i,2]
  start=lava_bivar_sign[i,3]
  stop=lava_bivar_sign[i,4]
  t1=lava_bivar_sign[i,7]
  t2=lava_bivar_sign[i,8]

  tmp_t1 <- sumstats_list[[t1]] %>%
    filter(., CHR == chr & BP >= start & BP <= stop) %>%
    filter(., P < 1e-06)
  
  tmp_t2 <- sumstats_list[[t2]] %>%
    filter(., CHR == chr & BP >= start & BP <= stop) %>%
    filter(., P < 1e-06)

  if (nrow(tmp_t1) == 0) {
    lava_bivar_sign[i, "most_sign_SNP_Trait1"] <- "None"
  } else if (nrow(tmp_t1 > 0)) {
    min_rsid_t1 <- arrange(tmp_t1, P) %>% .[1, "SNP"]
    lava_bivar_sign[i, "most_sign_SNP_Trait1"] <- min_rsid_t1
  }

  if (nrow(tmp_t2) == 0) {
    lava_bivar_sign[i, "most_sign_SNP_Trait2"] <- "None"
  } else if (nrow(tmp_t2 > 0)) {
    min_rsid_t2 <- arrange(tmp_t2, P) %>% .[1, "SNP"]
    lava_bivar_sign[i, "most_sign_SNP_Trait2"] <- min_rsid_t2
  }
}

lava_bivar_sign <- lava_bivar_sign %>%
  separate(phen1, c("trait1", "sex1"), remove = FALSE, sep = "_") %>%
  separate(phen2, c("trait2", "sex2"), remove = FALSE, sep = "_") %>%
  filter(., sex1 == sex2) %>%
  dplyr::select(-trait1, -trait2, -sex1, -sex2)

cat("Step 1 is completed.\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Main - step 2.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Starting step 2...\n")
for (i in 1:nrow(lava_bivar_sign)) {

    phen1 = lava_bivar_sign$phen1[i]
    phen2 = lava_bivar_sign$phen2[i]
    chr = lava_bivar_sign$chr[i]
    start = lava_bivar_sign$start[i]
    stop = lava_bivar_sign$stop[i]

    for (j in 1:length(sumstats_list)) {

        if (names(sumstats_list[j]) == phen1 | names(sumstats_list[j]) == phen2) {
	
	    head(sumstats_list[[j]])
            list_output[[row_n]] <- as.data.frame(sumstats_list[[j]]) %>%
                filter(., CHR == chr & BP >= start & BP <= stop) %>%
		mutate(phenotype = names(sumstats_list[j])) %>%
                arrange(., P) %>%
                head(., 1)
            row_n = row_n + 1

        }

    }

}

cat("Step 2 completed.\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Write outputs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

cat("Writing output files...\n")

if (arg_PD_proxies == 1) {

  # step 1.
  write.table(lava_bivar_sign, here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_sexstr_bivar_p_0.05_with_GWS_info.tsv"), sep = "\t", row.names = F, quote = F)
  write.table(lava_univ_sign, here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_sexstr_univar_p_2e-05_with_rsid_pval_info.tsv"), sep = "\t", row.names = F, quote = F)

  # step 2.
  rbindlist(list_output) %>%
    write.table(., here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_GWASsumstats_significant_info.tsv"), sep = "\t", row.names = F, quote = F)

} else if (arg_PD_proxies == 2) {

  # step 1.
  write.table(lava_bivar_sign, here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_noUKB_sexstr_bivar_p_0.05_with_GWS_info.tsv"), sep = "\t", row.names = F, quote = F)
  write.table(lava_univ_sign, here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_noUKB_sexstr_univar_p_2e-05_with_rsid_pval_info.tsv"), sep = "\t", row.names = F, quote = F)

  # step 2.
  rbindlist(list_output) %>%
    write.table(., here(project_dir, "results", "lava", "gwas_sexstr", "neurodegen_noUKB_GWASsumstats_significant_info.tsv"), sep = "\t", row.names = F, quote = F)
}

cat("Done.\n")