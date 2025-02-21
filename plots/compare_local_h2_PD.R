# Description: compare heritability of PD with and without proxies

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

n_loci = 2495 # 2889 if include chr X

# Read files -----------------------------------------------------

PD_h2 <- read.table(here(project_dir, "results", "lava_neurodegen", "neurodegen_sexstr_univ_all.tsv"), sep = "\t", header = TRUE) %>%
  filter(., phen %in% c("PD_females", "PD_males")) %>%
  mutate(phen = stringr::str_remove(phen, "PD_"),
         p_PD_status = case_when(
           p < 0.05/n_loci ~ "significant",
           p >= 0.05/n_loci ~ "not significant"
         )) %>%
  dplyr::select(
    locus, phen, h2_PD = h2.obs, p_PD = p, p_PD_status
  )

PDnp_h2 <- read.table(here(project_dir, "results", "lava_neurodegen", "neurodegen_noUKB_sexstr_univ_all.tsv"), sep = "\t", header = TRUE) %>%
  filter(., phen %in% c("PDnp_females", "PDnp_males")) %>%
  mutate(phen = stringr::str_remove(phen, "PDnp_"),
         p_PDnp_status = case_when(
           p < 0.05/2589 ~ "significant",
           p >= 0.05/2589 ~ "not significant"
         )) %>%
  dplyr::select(
    locus, phen, h2_PDnp = h2.obs, p_PDnp = p, p_PDnp_status
  )

# Main -----------------------------------------------------

h2 <- inner_join(PD_h2, PDnp_h2, by = c("locus", "phen")) %>%
  # filter(., h2_PDnp < 0.003) %>% # remove outlier values corresponding to the SNCA gene locus (681)
  mutate(p_status = case_when(
    p_PD_status == "not significant" & p_PDnp_status == "not significant" ~ "none significant",
    p_PD_status == "significant" & p_PDnp_status == "significant" ~ "both significant (N= 58)",
    p_PD_status == "not significant" & p_PDnp_status == "significant" ~ "significant without UKB (N = 102)",
    p_PD_status == "significant" & p_PDnp_status == "not significant" ~ "significant with UKB (N = 66)"
  )) %>%
  filter(., p_status != "none significant")

ggplot(h2, aes(x = h2_PD, y = h2_PDnp)) + 
  # geom_abline(intercept = 0, linetype = 2, colour = "grey") +
  geom_point(alpha = 0.8, size = 4, aes(shape = phen, fill = as.factor(p_status))) +
  scale_fill_manual(values = c("orange", "purple", "green")) +
  scale_shape_manual("Sex", values = c(21, 22)) +
  # scale_x_continuous(limits=c(0,0.01)) +
  # scale_y_continuous(limits=c(0,0.01)) +
  guides(fill = guide_legend("Significance status", override.aes = list(shape = 21))) +
  theme_bw() +
  labs(x = "Observed local h2 for PD (including UKB)", y = "Observed local h2 for PD (excluding UKB)") 
  
ggsave(here(project_dir, "results", "lava_neurodegen", "compare_local_h2_PD_with_and_without_UKB.jpg"), width = 20, height = 12, units = "cm")



