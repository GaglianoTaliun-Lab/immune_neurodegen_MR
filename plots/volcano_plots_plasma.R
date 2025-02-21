# Libraries ------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggplot2)
library(cowplot)
library(gridExtra)
library(ggrepel)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

# Read files -----------------------------------------------------

coloc_sign_plasma <- list.files(here(project_dir, "results", "cis_MR_plasma", "coloc"), pattern = "significant_PPH4_0.7_*", full.names = TRUE) %>% sort()

coloc_names <- coloc_sign_plasma  %>%
  stringr::str_remove(., "/Users/frida/Documents/research-projects/als_immune/results/cis_MR_plasma/coloc/significant_PPH4_0.7_") %>% 
  stringr::str_remove(., "_immune_exposures.tsv") %>% stringr::str_replace_all(., "_", " ") %>%
  stringr::str_replace(., "meta sexcombined", "sex-combined")

coloc_names_out <- coloc_sign_plasma  %>%
  stringr::str_remove(., "/Users/frida/Documents/research-projects/als_immune/results/cis_MR_plasma/coloc/significant_PPH4_0.7_") %>% 
  stringr::str_remove(., "_immune_exposures.tsv") %>%
  stringr::str_replace(., "meta_sexcombined", "sex-combined")

loop_names <- coloc_names_out %>% stringr::str_replace(., "sex-combined", "meta_sexcombined")

for (i in 1:length(coloc_sign_plasma)) {
  
  coloc <- read.table(stringr::str_c(project_dir, "/results/cis_MR_plasma/coloc/significant_PPH4_0.7_", loop_names[i], "_immune_exposures.tsv"), sep = "\t", header = T)
  
  MR_results <- read.table(stringr::str_c(project_dir, "/results/cis_MR_plasma/main/MR_main_results_immune_pqtls_UKB_", loop_names[i], ".txt"), sep = "\t", header = T) %>%
    mutate(OR = round(exp(b), 2), 
           OR_dir = case_when(
             OR > 1 & FDR < 0.05 ~ "Increase",
             OR < 1 & FDR < 0.05 ~ "Decrease",
             FDR >= 0.05 ~ "Not significant",
             OR == 1 ~ "Not significant"),
           FDR_significant = case_when(
             FDR < 0.05 ~ "YES",
             FDR >= 0.05 ~ "NO"
           ),
           colocalized = case_when(
             exposure %in% coloc$exposure ~ "YES",
             !(exposure %in% coloc$exposure) ~ "NO"),
           gene_symbol = stringr::str_extract(exposure, "[:alnum:]+_") %>% stringr::str_remove(., "_"),
           shape_plot = case_when(
             exposure %in% coloc$exposure ~ "FDR < 0.05 + colocalized",
             FDR_significant == "YES" & heterogeneity_signal == "NO" & MR_Egger_pleiotropy == "NO" ~ "FDR < 0.05",
             FDR_significant == "YES" & is.na(heterogeneity_signal)  & is.na(MR_Egger_pleiotropy) ~ "FDR < 0.05 + H and/or HP evidence",
             FDR_significant == "YES" & heterogeneity_signal == "YES" & MR_Egger_pleiotropy == "NO" ~ "FDR < 0.05 + H and/or HP evidence",
             FDR_significant == "YES" & heterogeneity_signal == "NO" & MR_Egger_pleiotropy == "YES" ~ "FDR < 0.05 + H and/or HP evidence",
             FDR_significant == "YES" & heterogeneity_signal == "YES" & MR_Egger_pleiotropy == "YES" ~ "FDR < 0.05 + H and/or HP evidence",
             FDR_significant == "YES" & heterogeneity_signal == "YES"  & is.na(MR_Egger_pleiotropy) ~ "FDR < 0.05 + H and/or HP evidence",
             FDR_significant == "YES" & is.na(heterogeneity_signal) & is.na(MR_Egger_pleiotropy) ~ "FDR < 0.05",
             FDR_significant == "YES" & heterogeneity_signal == "NO" & is.na(MR_Egger_pleiotropy) ~ "FDR < 0.05",
             FDR_significant == "YES" & is.na(heterogeneity_signal) & MR_Egger_pleiotropy == "NO" ~ "FDR < 0.05",
             FDR_significant == "NO" ~ "FDR >= 0.05"
           ),
           shape_plotord = factor(shape_plot, 
                                  levels = c("FDR >= 0.05", "FDR < 0.05 + H and/or HP evidence", "FDR < 0.05", "FDR < 0.05 + colocalized")))
  
  sign_genes <- read.table(stringr::str_c(project_dir, "/results/cis_MR_plasma/main/MR_significant_main_results_immune_pqtls_UKB_", loop_names[i], ".txt"), sep = "\t", header = T) %>%
    mutate(OR = round(exp(b), 2), 
           OR_dir = case_when(
             OR > 1 & FDR < 0.05 ~ "increase",
             OR < 1 & FDR < 0.05 ~ "decrease",
             FDR >= 0.05 ~ "ns",
             OR == 1 ~ "ns"),
           gene_symbol = stringr::str_extract(exposure, "[:alnum:]+_") %>% stringr::str_remove(., "_")) %>%
    filter(., exposure %in% coloc$exposure)
  
  ggplot(MR_results, aes(x = OR, y = -log10(pval))) + 
    geom_point(alpha = 0.8, size = 4, aes(shape = shape_plotord, fill = OR_dir)) +
    # geom_point(data = sign_genes, colour = "black", alpha = 1, fill = "red", aes(shape = shape_plot)) +
    geom_vline(xintercept = 1,
               linetype = "dashed") +
    xlim(0.5, max(MR_results$OR)) +
    scale_fill_manual("Direction of odds ratio", values = c("orange", "purple", "grey"), 
                      breaks = c("Decrease", "Increase", "Not significant"), 
                      labels = c("Decrease", "Increase", "Not significant")) +
    scale_shape_manual("Significance", values = c("circle filled", "square filled", "triangle filled", "diamond filled"),
                       breaks = c("FDR >= 0.05", "FDR < 0.05 + H and/or HP evidence", "FDR < 0.05", "FDR < 0.05 + colocalized"),
                       labels = c("FDR >= 0.05", "FDR < 0.05 + H and/or HP evidence", "FDR < 0.05", "FDR < 0.05 + colocalized")) +
    theme(axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20)) +
    guides(fill = guide_legend("Direction of odds ratio", override.aes = list(shape = 21))) +
    theme_bw() +
    labs(x = stringr::str_c("Odds Ratio of ", coloc_names[i]), y = "-log10(p-value)") +
    geom_label_repel(data = sign_genes, aes(label = gene_symbol), fill = "white", alpha= 0.7, force = 2, nudge_y = 0.1, nudge_x = 0.1)
  
  ggsave(here(project_dir, "results", "cis_MR_plasma", "figures", stringr::str_c("volcano_plot_", loop_names[i],".jpg")),
         width = 22, height = 12, units = "cm")
  
}
