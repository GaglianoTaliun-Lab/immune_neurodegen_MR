# Description: figure summarizing plasma cisMR results and validation

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

# Read files -----------------------------------------------------

# 1) ALL MR main results files
mr_files <- list.files(here(project_dir, "results", "cis_MR_plasma", "main"), pattern = "MR_main_results_immune_pqtls_UKB*", full.names = T)
mr_results <- lapply(mr_files, function(x) read.table(x, sep = "\t", header = T) %>%
                       separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
                       mutate(exp_trait = stringr::str_c(exposure, "_", Disease))) %>%
  rbindlist(., use.names = TRUE, fill = TRUE)

# 2) MR Egger main results files
mr_egger_files <- list.files(here(project_dir, "results", "cis_MR_plasma", "main"), pattern = "MR_Egger_main_results_immune_pqtls_UKB*", full.names = T)
mr_egger <- lapply(mr_egger_files, function(x) read.table(x, sep = "\t", header = T) %>%
                     separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
                     mutate(exp_trait = stringr::str_c(exposure, "_", Disease),
                            sex = case_when(
                              sex == "females" ~ "Female-specific",
                              sex == "males" ~ "Male-specific",
                              sex == "meta" ~ "Sex-combined"
                            ),
                            exp_trait = stringr::str_c(exposure, "_", outcome) %>% 
                              stringr::str_remove(., "_meta_sexcombined") %>%
                              stringr::str_remove(., "_females") %>%
                              stringr::str_remove(., "_males"),
                            exp_out = stringr::str_c(exposure, "_", outcome))) %>%
  rbindlist() %>%
  dplyr::select(., exp_trait, 
                exp_out, 
                outcome, 
                Disease, 
                sex, 
                exposure,
                pval_egger = pval)

# 3) COLOC files
coloc_files <- setNames(list.files(here(project_dir, "results", "cis_MR_plasma", "coloc"), pattern = "significant_PPH4_0.7_*", full.names = T),
                        nm = list.files(here(project_dir, "results", "cis_MR_plasma","coloc"), pattern = "significant_PPH4_0.7_*", full.names = F) %>%
                          stringr::str_remove_all(., "significant_PPH4_0.7_") %>% stringr::str_remove_all(., "_immune_exposures.tsv")
)

coloc_results <- lapply(coloc_files, function(x) read.table(x, sep = "\t", header = T)) %>%
  vctrs::list_drop_empty() %>%
  dplyr::bind_rows(., .id = "outcome") %>%
  separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
  mutate(sex = case_when(
    sex == "females" ~ "Female-specific",
    sex == "males" ~ "Male-specific",
    sex == "meta" ~ "Sex-combined"
  ),
  PP.H4.abf = round(PP.H4.abf, 2),
  exp_trait = stringr::str_c(exposure, "_", outcome) %>% 
    stringr::str_remove(., "_meta_sexcombined") %>%
    stringr::str_remove(., "_females") %>%
    stringr::str_remove(., "_males"),
  exp_out = stringr::str_c(exposure, "_", outcome))

coloc_files_all <- setNames(list.files(here(project_dir, "results", "cis_MR_plasma", "coloc"), pattern = "*_immune_exposures.tsv", full.names = T) %>% stringr::str_subset(., "significant_PPH4", negate = TRUE),
                            nm = list.files(here(project_dir, "results", "cis_MR_plasma", "coloc"), pattern = "*_immune_exposures.tsv", full.names = F) %>% stringr::str_subset(., "significant_PPH4", negate = TRUE) %>%
                              stringr::str_remove_all(., "_immune_exposures.tsv")
)

coloc_results_all <- lapply(coloc_files_all, function(x) read.table(x, sep = "\t", header = T)) %>%
  vctrs::list_drop_empty() %>%
  dplyr::bind_rows(., .id = "outcome") %>%
  mutate(exp_out = stringr::str_c(exposure, "_", outcome), PP.H4.abf = round(PP.H4.abf, 2)) %>%
  dplyr::select(exp_out, PP.H4.abf)

# 4) Read validation results
aric <- read.table(here(project_dir, "results", "cis_MR_plasma", "validation", "ARIC_validation_results_all.txt"), sep = "\t", header = T) %>%
  dplyr::select(outcome, Gene_symbol = exposure, OR, pvalue, PP.H4.abf) %>%
  separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
  mutate(Disease = stringr::str_remove(Disease, "np"),
         analysis = "V1 (ARIC - SomaScan)",
         sex = case_when(
           sex == "females" ~ "Female-specific",
           sex == "males" ~ "Male-specific",
           sex == "meta" ~ "Sex-combined"
         ),
         sign = case_when(
           pvalue < 0.0045 ~ "sign",
           pvalue >= 0.0045 ~ "not_sign"
         ),
         PP.H4.abf = round(PP.H4.abf,2)) %>%
  dplyr::select(-outcome, -pvalue)

decode <- read.table(here(project_dir, "results", "cis_MR_plasma", "validation", "deCODE_validation_results_all.txt"), sep = "\t", header = T) %>%
  dplyr::select(outcome, Gene_symbol = exposure, OR, pvalue, PP.H4.abf) %>%
  separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
  mutate(Disease = stringr::str_remove(Disease, "np"),
         analysis = "V2 (deCODE - SomaScan)",
         sex = case_when(
           sex == "females" ~ "Female-specific",
           sex == "males" ~ "Male-specific",
           sex == "meta" ~ "Sex-combined"
         ),
         sign = case_when(
           pvalue < 0.0045 ~ "sign",
           pvalue >= 0.0045 ~ "not_sign"
         ),
         PP.H4.abf = round(PP.H4.abf,2)) %>%
  dplyr::select(-outcome, -pvalue)


# Main -----------------------------------------------------------

# note: function of TwoSampleMR "generate_odds_ratios()" can be used instead of computing it manually. It gives same results.

mr_results2 <- mr_results %>%
  mutate(OR = round(exp(b),2), 
         low_CI = round(exp(b-1.96*se),2), 
         upp_CI = round(exp(b+1.96*se),2), 
         OR_CI = str_c(OR, " (", low_CI, "-", upp_CI, ")"),
         sex = case_when(
           sex == "females" ~ "Female-specific",
           sex == "males" ~ "Male-specific",
           sex == "meta" ~ "Sex-combined"
         ),
         exp_out = stringr::str_c(exposure, "_", outcome)) %>%
  filter(., exp_trait %in% coloc_results$exp_trait) %>%
  left_join(., coloc_results, by = c("exp_trait", "exp_out", "outcome", "Disease", "sex", "exposure")) %>%
  left_join(., mr_egger, by = c("exp_trait", "exp_out", "outcome", "Disease", "sex", "exposure")) %>%
  dplyr::select(-PP.H4.abf) %>%
  left_join(., coloc_results_all, by = "exp_out") %>%
  mutate(PP.H4.abf = case_when(
    is.na(PP.H4.abf) ~ "-",
    !is.na(PP.H4.abf) ~ as.character(PP.H4.abf)
  ))

# add tiers according to validation status:
tiers <- fread(here(project_dir, "results", "cis_MR_plasma", "validation", "tiers_forest_plot.txt")) %>%
  separate(exposure, c("Gene_symbol", NA, NA), sep = "_") %>%
  mutate(Disease = stringr::str_remove(Disease, "np"))

# add SDE:
sde_female <- mr_results2 %>%
  filter(sex == "Female-specific") %>%
  dplyr::select(Disease, exposure, beta_female = b, se_female = se)

sde_male <- mr_results2 %>%
  filter(sex == "Male-specific") %>%
  dplyr::select(Disease, exposure, beta_male = b, se_male = se)

sde_test <- inner_join(sde_female, sde_male, by = c("Disease", "exposure")) %>%
  mutate(zscore = (beta_female - beta_male)/(sqrt(se_female^2 + se_male^2)),
         zscore_pval = 2*pnorm(q=abs(zscore), lower.tail = FALSE),
         fdr = p.adjust(zscore_pval, method = "BH"),
         SDE = case_when(
           fdr < 0.05 ~ "**",
           fdr >= 0.05 ~ "-"
         )) %>% 
  separate(exposure, c("Gene_symbol", "Exposure_Olink_ID", "Exposure_Uniprot_ID"), extra = "drop", remove = T) %>%
  dplyr::select(Disease, Gene_symbol, SDE) %>%
  filter(., Disease %in% c("ADnp", "PDnp", "ALS")) %>%
  mutate(Disease = case_when(
    Disease == "ADnp" ~ "AD",
    Disease == "PDnp" ~ "PD",
    Disease == "ALS" ~ "ALS"
  ))

mr_results_plot <- mr_results2 %>%
  filter(., stringr::str_detect(outcome, "AD_", negate = TRUE)) %>%
  filter(., stringr::str_detect(outcome, "PD_", negate = TRUE)) %>%
  dplyr::select(-id.exposure, -id.outcome, -method, -nsnp, -het_Qpval, -heterogeneity_signal,
                -MR_Egger_intercept_pval, -MR_Egger_pleiotropy, -b, -se, -pval, -low_CI, -upp_CI,
                -OR_CI, -nsnps, -PP.H0.abf, -PP.H1.abf, -PP.H2.abf, -PP.H3.abf, -pval_egger) %>%
  mutate(Disease = case_when(
    Disease == "ADnp" ~ "AD",
    Disease == "PDnp" ~ "PD",
    TRUE ~ Disease),
    sign = case_when(
      FDR < 0.05 ~ "sign",
      FDR >= 0.05 ~ "not_sign"
    ),
    analysis = "Primary (UKB - Olink)") %>%
  separate(exposure, c("Gene_symbol", "OlinkID", "UniProtID"), sep = "_") %>%
  dplyr::select(Disease, sex, Gene_symbol, OR, PP.H4.abf, analysis, sign)

plot_all <- rbind(mr_results_plot, aric, decode) %>%
  mutate(ord_group = factor(analysis, 
                            levels = c("Primary (UKB - Olink)","V1 (ARIC - SomaScan)", "V2 (deCODE - SomaScan)")),
         PP.H4.abf = case_when(
           is.na(PP.H4.abf) ~ "-",
           !is.na(PP.H4.abf) ~ PP.H4.abf
         )) %>%
  left_join(., tiers, by = c("Disease", "Gene_symbol")) %>%
  left_join(., sde_test, by = c("Disease", "Gene_symbol")) %>%
  mutate(tier_factor = factor(tier, levels = c(1, 2, 3),
                             labels = c("Tier 1", "Tier 2", "Tier 3")))

# Create tier indicator data as a separate "sex" category
tier_indicators <- plot_all %>%
  filter(analysis == "Primary (UKB - Olink)") %>%
  select(Gene_symbol, Disease, tier) %>%
  distinct() %>%
  mutate(sex = "Tier",
         ord_group = factor("Tier", levels = c("Primary (UKB - Olink)", "V1 (ARIC - SomaScan)", "V2 (deCODE - SomaScan)", "Tier")),
         tier_color = case_when(
           tier == 1 ~ "#88419D",
           tier == 2 ~ "#8C96C6",
           tier == 3 ~ "#BFD3E6",
           TRUE ~ "black"
         ))

# Create SDE indicator data as a separate "sex" category
sde_indicators <- plot_all %>%
  filter(analysis == "Primary (UKB - Olink)") %>%
  select(Gene_symbol, Disease, SDE) %>%
  distinct() %>%
  mutate(sex = "SDE",
         ord_group = factor("SDE", levels = c("Primary (UKB - Olink)", "V1 (ARIC - SomaScan)", "V2 (deCODE - SomaScan)", "SDE")),
         sde_color = "white")

# Add tier/SDE indicators to the main data
# First check what columns are in plot_all
tier_indicators_full <- tier_indicators %>% 
  mutate(OR = NA, 
         PP.H4.abf = NA,
         analysis = "Tier",
         sign = NA,
         SDE = NA,
         tier_factor = NA)

sde_indicators_full <- sde_indicators %>% 
  mutate(OR = NA, 
         PP.H4.abf = NA,
         analysis = "SDE",
         sign = NA,
         tier = NA,
         tier_factor = NA)

# Use bind_rows to put all together
plot_all_with_tiers <- dplyr::bind_rows(plot_all, tier_indicators_full) %>%
  dplyr::bind_rows(., sde_indicators_full)

# Update sex factor levels to include "Tier" and "SDE"
plot_all_with_tiers <- plot_all_with_tiers %>%
  mutate(sex = factor(sex, levels = c("Tier", "SDE", "Female-specific", "Male-specific", "Sex-combined"))) %>%
  arrange(Disease, -tier, desc(Gene_symbol)) %>%
  mutate(Gene_symbol = factor(Gene_symbol, levels = unique(Gene_symbol)))

ggplot(plot_all_with_tiers, aes(x=Gene_symbol, y=ord_group)) +
  # tiers
  geom_point(data = filter(plot_all_with_tiers, sex == "Tier"),
             aes(fill = tier_color),
             shape = 22, size = 8, colour = "black", stroke = 0.5) +
  # SDE
  geom_point(data = filter(plot_all_with_tiers, sex == "SDE"), fill = "white",
             shape = 22, size = 8, colour = "black", stroke = 0.5) +
  # Main plot points for other facets
  geom_point(data = filter(plot_all_with_tiers, !(sex %in% c("Tier", "SDE"))),
             aes(fill = sign), shape = 21, colour = "white", size = 10) +
  geom_point(data = filter(plot_all_with_tiers, !(sex %in% c("Tier", "SDE"))),
             aes(colour = log(OR)), shape = 16, size = 8) +
  theme_bw() +
  scale_fill_manual(values = c("not_sign" = "grey", "sign" = "black",
                               "#88419D" = "#88419D", "#8C96C6" = "#8C96C6", "#BFD3E6" = "#BFD3E6"),
                    name = "Significance",
                    labels = c("not_sign" = "Not significant", "sign" = "Significant"),
                    breaks = c("not_sign", "sign"),
                    guide = guide_legend(override.aes = list(shape = 21, size = 5, stroke = 1.5, fill = "white", colour = c("grey", "black")))) +
  scale_colour_distiller(
    palette = "RdYlBu",
    limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1), 
    na.value = "grey",
    name = "cis-MR effect size") +
  scale_size_continuous(range = c(3, 12), guide = "none") +
  labs(x = "", y = "pQTLs dataset", title = "") +
  theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1, colour = "black"),
        axis.text.y = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 16, colour = "black"),
        axis.title.y = element_text(size = 16, colour = "black"),
        legend.title = element_text(size = 14, colour = "black"),
        legend.text = element_text(size = 12, colour = "black"),
        strip.text.x = element_text(size = 14, colour = "black"),
        strip.text.y = element_text(size = 12, colour = "black"),
        panel.grid.minor = element_blank(),
        strip.placement = "outside",
        plot.caption = element_text(size = 10, hjust = 0),
        axis.ticks.x = element_line(colour = "white")) +
  scale_y_discrete(expand = expansion(mult = c(0.5, 0.5)),
                   labels = function(x) ifelse(x %in% c("Tier", "SDE"), "", x)) +
  geom_text(data = filter(plot_all_with_tiers, !(sex %in% c("Tier", "SDE"))),
            aes(label=PP.H4.abf), size = 2.7, vjust = 0.5) +
  geom_text(data = filter(plot_all_with_tiers, sex == "Tier"),
            aes(label=tier), size = 3.5, vjust = 0.5) +
  geom_text(data = filter(plot_all_with_tiers, sex == "SDE"),
            aes(label=SDE), size = 5, vjust = 0.65) +
  coord_flip() +
  facet_grid(Disease ~ sex, scales = "free", space = "free", switch = "y")

ggsave(here(project_dir, "results", "cis_MR_plasma", "figures", "summary_plot_neurodegen_plasma_pqtls.tif"), width=12, height=8, dpi = 300)
