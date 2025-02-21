# plot global correlations between sex stratified ALS and PD

# Description: script to obtain heat map plots of the global genetic correlations across traits.

# Libraries ----------------------------------------------------------------------------------------

library(tidyverse)
library(here)
library(gghighlight)
library(stringr)
library(reshape2)
library(corrplot)

# Arguments ----------------------------------------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

correlations <- read.table(here(project_dir, "results", "ldsc_neurodegen", "ldsc_correlations_sexstr_neurodegen.txt"), sep = "\t", header = T)

p_nominal = 0.05

n_traits_per_plot = 9

source(here(project_dir, "scripts", "global_corr_function_RHR.R"))

# Wrangle data ----------------------------------------------------------------------------------------

# Including UKB (and UKB proxies) in PD and AD datasets:
correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         p1 =
           case_when(
             p1 == "ALS_males" ~ "ALS (M)",
             p1 == "ALS_females" ~ "ALS (F)",
             p1 == "PD_males" ~ "PD (M)",
             p1 == "PD_females" ~ "PD (F)",
             p1 == "AD_males" ~ "AD (M)",
             p1 == "AD_females" ~ "AD (F)",
             p1 == "AD_meta_sexcombined" ~ "AD (SC)",
             p1 == "PD_meta_sexcombined" ~ "PD (SC)",
             p1 == "ALS_meta_sexcombined" ~ "ALS (SC)"
             ),
         p2 =
           case_when(
             p2 == "ALS_males" ~ "ALS (M)",
             p2 == "ALS_females" ~ "ALS (F)",
             p2 == "PD_males" ~ "PD (M)",
             p2 == "PD_females" ~ "PD (F)",
             p2 == "AD_males" ~ "AD (M)",
             p2 == "AD_females" ~ "AD (F)",
             p2 == "AD_meta_sexcombined" ~ "AD (SC)",
             p2 == "PD_meta_sexcombined" ~ "PD (SC)",
             p2 == "ALS_meta_sexcombined" ~ "ALS (SC)"
             ),
         ord_p1 = factor(p1, levels = c("AD (M)", "AD (F)", "AD (SC)", "ALS (F)", "ALS (M)", "ALS (SC)", "PD (F)","PD (M)", "PD (SC)")),
         ord_p2 = factor(p2, levels = c("AD (M)", "AD (F)", "AD (SC)", "ALS (F)", "ALS (M)", "ALS (SC)", "PD (F)","PD (M)", "PD (SC)"))
  ) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2))

write.table(correlations_plot %>% dplyr::select(-ord_p1, -ord_p2), here(project_dir, "results", "ldsc_neurodegen", "manuscript_SuppTable_with_proxies.txt"), sep = "\t", row.names = F, quote = F)
plot_global_corr(global_corr = correlations_plot, n_phenotypes = n_traits_per_plot)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen.jpg"), width = 45, height = 35, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen.pdf"), width = 45, height = 35, units = "cm")

##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Excluding UKB (and proxies) from PD and AD cohort:
correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         p1 =
           case_when(
             p1 == "ALS_males" ~ "ALS (M)",
             p1 == "ALS_females" ~ "ALS (F)",
             p1 == "PDnp_males" ~ "PD (M)",
             p1 == "PDnp_females" ~ "PD (F)",
             p1 == "ADnp_males" ~ "AD (M)",
             p1 == "ADnp_females" ~ "AD (F)",
             p1 == "ADnp_meta_sexcombined" ~ "AD (SC)",
             p1 == "PDnp_meta_sexcombined" ~ "PD (SC)",
             p1 == "ALS_meta_sexcombined" ~ "ALS (SC)"),
         p2 =
           case_when(
             p2 == "ALS_males" ~ "ALS (M)",
             p2 == "ALS_females" ~ "ALS (F)",
             p2 == "PDnp_males" ~ "PD (M)",
             p2 == "PDnp_females" ~ "PD (F)",
             p2 == "ADnp_males" ~ "AD (M)",
             p2 == "ADnp_females" ~ "AD (F)",
             p2 == "ADnp_meta_sexcombined" ~ "AD (SC)",
             p2 == "PDnp_meta_sexcombined" ~ "PD (SC)",
             p2 == "ALS_meta_sexcombined" ~ "ALS (SC)"),
         ord_p1 = factor(p1, levels = c("AD (M)", "AD (F)", "AD (SC)", "ALS (F)", "ALS (M)", "ALS (SC)", "PD (F)","PD (M)", "PD (SC)")),
         ord_p2 = factor(p2, levels = c("AD (M)", "AD (F)", "AD (SC)", "ALS (F)", "ALS (M)", "ALS (SC)", "PD (F)","PD (M)", "PD (SC)"))
  ) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2))

write.table(correlations_plot %>% dplyr::select(-ord_p1, -ord_p2), here(project_dir, "results", "ldsc_neurodegen", "manuscript_SuppTable_without_UKB_at_all.txt"), sep = "\t", row.names = F, quote = F)
plot_global_corr(global_corr = correlations_plot, n_phenotypes = n_traits_per_plot)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_no_UKB_at_all.jpg"), width = 45, height = 35, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_no_UKB_at_all.pdf"), width = 45, height = 35, units = "cm")

##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Including AD and PD UKB (and proxies), but excluding sex-combined results:

n_traits_per_plot = 6

correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         p1 =
           case_when(
             p1 == "ALS_males" ~ "ALS (M)",
             p1 == "ALS_females" ~ "ALS (F)",
             p1 == "PD_males" ~ "PD (M)",
             p1 == "PD_females" ~ "PD (F)",
             p1 == "AD_males" ~ "AD (M)",
             p1 == "AD_females" ~ "AD (F)"
           ),
         p2 =
           case_when(
             p2 == "ALS_males" ~ "ALS (M)",
             p2 == "ALS_females" ~ "ALS (F)",
             p2 == "PD_males" ~ "PD (M)",
             p2 == "PD_females" ~ "PD (F)",
             p2 == "AD_males" ~ "AD (M)",
             p2 == "AD_females" ~ "AD (F)"
           ),
         ord_p1 = factor(p1, levels = c("AD (M)", "AD (F)", "ALS (F)", "ALS (M)", "PD (F)","PD (M)")),
         ord_p2 = factor(p2, levels = c("AD (M)", "AD (F)", "ALS (F)", "ALS (M)", "PD (F)","PD (M)"))
  ) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2))

plot_global_corr(global_corr = correlations_plot, n_phenotypes = n_traits_per_plot)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noSC.jpg"), width = 45, height = 35, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noSC.pdf"), width = 45, height = 35, units = "cm")

##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Excluding UKB (and proxies) from AD and PD, and excluding sex-combined results:

n_traits = 6
indices=c(1,3,5,8,9,11,14,16,17)
p_thres <- 0.05/n_traits

correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         sex1 = case_when(
           str_detect(p1, "_males") ~ "M",
           str_detect(p1, "_females") ~ "F"),
         sex2 = case_when(
           str_detect(p2, "_males") ~ "M",
           str_detect(p2, "_females") ~ "F"),
         p1 =
           case_when(
             p1 == "ALS_males" ~ "ALS",
             p1 == "ALS_females" ~ "ALS",
             p1 == "PDnp_males" ~ "PD",
             p1 == "PDnp_females" ~ "PD",
             p1 == "ADnp_males" ~ "AD",
             p1 == "ADnp_females" ~ "AD"
           ),
         p2 =
           case_when(
             p2 == "ALS_males" ~ "ALS",
             p2 == "ALS_females" ~ "ALS",
             p2 == "PDnp_males" ~ "PD",
             p2 == "PDnp_females" ~ "PD",
             p2 == "ADnp_males" ~ "AD",
             p2 == "ADnp_females" ~ "AD"
           ),
         ord_p1 = factor(p1, levels = c("AD", "ALS", "PD")),
         ord_p2 = factor(p2, levels = c("AD", "ALS", "PD"))
  ) %>%
  filter(., sex1 == sex2) %>%
  mutate(sign = case_when(
    p < p_thres ~ "**",
    p < 0.05 ~ "*",
    p >= 0.05 ~ "")) %>%
  tidyr::spread(., sex1, rg) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2)) %>%
  mutate(M = case_when(p1 == p2 ~ as.numeric(NA),
                           TRUE ~ M),
         pheno1 = stringr::str_c(p1, " (", sex2, ")"), pheno2 = str_c(p2, " (", sex2, ")"),
         ord_pheno1 = factor(pheno1, levels = c("AD (F)", "AD (M)", "ALS (F)", "ALS (M)", "PD (F)","PD (M)")),
         ord_pheno2 = factor(pheno2, levels = c("AD (F)", "AD (M)", "ALS (F)", "ALS (M)", "PD (F)","PD (M)"))) %>%
  dplyr::rename(rg_females = F, rg_males = M) %>%
  dplyr::slice(indices)

correlations_plot %>%
  mutate(is_upper = str_detect(sex2, "F")) %>%
  mutate(is_diag = (p1 == p2 & sex2 == "F")) %>%
  ggplot(aes(ord_p1, ord_p2 %>% fct_rev())) +
  geom_tile(aes(fill = ifelse(is_upper, rg_females, rg_males)), colour = "black") +
  geom_text(aes(label = ifelse(is_upper, rg_females, rg_males)), size = 10) +
  geom_text(aes(label=sign), size = 10, vjust = 1.95) + #added line
  coord_fixed() +
  labs(x = "", y = "", fill = "Correlation (rg)") +
  scale_fill_gradient2(low = "purple",
                                high = "red3", na.value = "grey",
                                limits = c(-1, 1)) +
  scale_colour_manual(values = c("black", "white")) +
  # geom_text(aes(label=sign), size = 10, vjust = 1.95) + #added line
  guides(colour = "none") +
  theme_classic() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 28, angle = 45, vjust = 1, hjust = 1, colour = "black"),
    axis.text.y = element_text(size = 28, hjust = 0.5, vjust = 0.5, colour = "black"),
    axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 22),
    legend.text = element_text(size = 22)
  )

# write.table(correlations_plot, here(project_dir, "results", "ldsc_neurodegen", "manuscript_SuppTable.txt"), sep = "\t", row.names = F, quote = F)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_no_UKB_at_all_no_sexcomb.jpg"), width = 25, height = 15, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_no_UKB_at_all_no_sexcomb.pdf"), width = 25, height = 15, units = "cm")

#### correlations by sex separately ##### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Excluding UKB (and proxies) from AD and PD and excluding sex-combined results:

n_traits_per_plot = 3

# females:
correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         p1 =
           case_when(
             p1 == "ALS_females" ~ "ALS (F)",
             p1 == "PDnp_females" ~ "PD (F)",
             p1 == "ADnp_females" ~ "AD (F)"),
         p2 =
           case_when(
             p2 == "ALS_females" ~ "ALS (F)",
             p2 == "PDnp_females" ~ "PD (F)",
             p2 == "ADnp_females" ~ "AD (F)"),
         ord_p1 = factor(p1, levels = c("AD (F)", "ALS (F)", "PD (F)")),
         ord_p2 = factor(p2, levels = c("AD (F)", "ALS (F)", "PD (F)"))
  ) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2))

write.table(correlations_plot, here(project_dir, "results", "ldsc_neurodegen", "manuscript_SuppTable_females_no_UKB_at_all.txt"), sep = "\t", row.names = F, quote = F)
plot_global_corr(global_corr = correlations_plot, n_phenotypes = n_traits_per_plot)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noPDproxies_no_UKB_at_all_females_only.jpg"), width = 45, height = 35, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noPDproxies_no_UKB_at_all_females_only.pdf"), width = 45, height = 35, units = "cm")

# males:
correlations_plot <- correlations %>%
  mutate(rg = round(rg, digits = 2),
         p1 =
           case_when(
             p1 == "ALS_males" ~ "ALS (M)",
             p1 == "PDnp_males" ~ "PD (M)",
             p1 == "ADnp_males" ~ "AD (M)"),
         p2 =
           case_when(
             p2 == "ALS_males" ~ "ALS (M)",
             p2 == "PDnp_males" ~ "PD (M)",
             p2 == "ADnp_males" ~ "AD (M)"),
         ord_p1 = factor(p1, levels = c("AD (M)", "ALS (M)", "PD (M)")),
         ord_p2 = factor(p2, levels = c("AD (M)", "ALS (M)", "PD (M)"))
  ) %>%
  arrange(ord_p1, as.factor(ord_p1), as.factor(ord_p2)) %>%
  filter(., !is.na(p1)) %>%
  filter(., !is.na(p2))

write.table(correlations_plot, here(project_dir, "results", "ldsc_neurodegen", "manuscript_SuppTable_males_no_UKB_at_all.txt"), sep = "\t", row.names = F, quote = F)
plot_global_corr(global_corr = correlations_plot, n_phenotypes = n_traits_per_plot)
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noPDproxies_no_UKB_at_all_males_only.jpg"), width = 45, height = 35, units = "cm")
ggsave(here(project_dir,"results","ldsc_neurodegen","global_correlations_sexstr_neurodegen_noPDproxies_no_UKB_at_all_males_only.pdf"), width = 45, height = 35, units = "cm")


