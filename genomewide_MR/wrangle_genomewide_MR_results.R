# Include sensitivity tests in genome-wide MR results for Supplementary Table

# Packages -------------------------------------------------------

library(tidyverse)
library(stringr)
library(data.table)
library(here)

# Arguments -------------------------------------------------------

project_dir=here("Documents/research-projects/als_immune/results")

# Read files ----------------------------------------------------

all_results <- read.table(here(project_dir, "genomewide_MR", "main", "MR_main_results_all_pairs_excl_conf_excl_HLA.txt"), sep = "\t", header = T) %>%
  filter(., method %in% c("Inverse variance weighted (multiplicative random effects)", "Weighted median", "MR Egger")) %>%
  dplyr::select(-id.exposure, -id.outcome)

# Main -----------------------------------------------------

ivw_results <- all_results %>%
  filter(., method == "Inverse variance weighted (multiplicative random effects)") %>%
  mutate(FDR = p.adjust(pval, method = "BH"), method = "Inverse variance weighted")

other_results <- all_results %>%
  filter(., method != "Inverse variance weighted (multiplicative random effects)") %>%
  mutate(FDR = "NA")

all_results <- rbind(ivw_results, other_results)

het <- list.files(here(project_dir, "genomewide_MR", "sensitivity"), pattern = "heterogeneity_*", full.names = T)
plei <- list.files(here(project_dir, "genomewide_MR", "sensitivity"), pattern = "pleiotropy_*", full.names = T)

het_df <- lapply(het, function(x) read.table(x, sep = "\t", header = T)) %>%
  rbindlist() %>%
  filter(., method == "Inverse variance weighted") %>%
  dplyr::select(-id.exposure, -id.outcome, method)

plei_df <- lapply(plei, function(x) read.table(x, sep = "\t", header = T)) %>%
  rbindlist() %>%
  mutate(method = "Inverse variance weighted") %>%
  dplyr::select(outcome, exposure, egger_intercept, egger_se = se, egger_pval = pval, method)

full_results <- left_join(all_results, het_df, by = c("outcome", "exposure", "method")) %>%
  left_join(., plei_df, by = c("outcome", "exposure", "method")) %>%
  filter(., se < 1) %>%
  mutate(outcome = case_when(
    outcome == "ALS_females" ~ "ALS females",
    outcome == "ALS_males" ~ "ALS males",
    outcome == "ALS_meta_sexcombined" ~ "ALS sex-combined",
    outcome == "PD_females" ~ "PD females",
    outcome == "PD_males" ~ "PD males",
    outcome == "PD_meta_sexcombined" ~ "PD sex-combined",
    outcome == "AD_females" ~ "AD females",
    outcome == "AD_males" ~ "AD males",
    outcome == "AD_meta_sexcombined" ~ "AD sex-combined",
    outcome == "PDnp_females" ~ "PD (excluding proxy cases) females",
    outcome == "PDnp_males" ~ "PD (excluding proxy cases) males",
    outcome == "PDnp_meta_sexcombined" ~ "PD (excluding proxy cases) sex-combined",
    outcome == "ADnp_females" ~ "AD (excluding UKB) females",
    outcome == "ADnp_males" ~ "AD (excluding UKB) males",
    outcome == "ADnp_meta_sexcombined" ~ "AD (excluding UKB) sex-combined",
  ),
  exposure = case_when(
    exposure == "eosinoph_meta_sexcombined" ~ "Eosinophils (sex-combined)",
    exposure == "eosinoph_females" ~ "Eosinophils (females)",
    exposure == "eosinoph_males" ~ "Eosinophils (males)",
    exposure == "eosinophper_meta_sexcombined" ~ "Eosinophils % (sex-combined)",
    exposure == "eosinophper_females" ~ "Eosinophils % (females)",
    exposure == "eosinophper_males" ~ "Eosinophils % (males)",
    exposure == "lymph_females" ~ "Lymphocytes (females)",
    exposure == "lymph_males" ~ "Lymphocytes (males)",
    exposure == "lymph_meta_sexcombined" ~ "Lymphocytes (sex-combined)",
    exposure == "mono_females" ~ "Monocytes % (females)",
    exposure == "mono_males" ~ "Monocytes % (males)",
    exposure == "mono_meta_sexcombined" ~ "Monocytes % (sex-combined)",
    exposure == "neutroph_males" ~ "Neutrophils (males)",
    exposure == "neutroph_females" ~ "Neutrophils (females)",
    exposure == "neutroph_meta_sexcombined" ~ "Neutrophils (sex-combined)",
    exposure == "neutrophper_males" ~ "Neutrophils % (males)",
    exposure == "neutrophper_females" ~ "Neutrophils % (females)",
    exposure == "neutrophper_meta_sexcombined" ~ "Neutrophils % (sex-combined)",
    exposure == "basophil_males" ~ "Basophils (males)",
    exposure == "basophil_females" ~ "Basophils (females)",
    exposure == "basophil_meta_sexcombined" ~ "Basophils (sex-combined)",
  )) %>%
  arrange(., method, pval) %>% dplyr::select(-Q_df)

write.table(full_results, here(project_dir, "genomewide_MR", "main", "MR_all_results_immune_cell_traits_with_sensitivity_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)
