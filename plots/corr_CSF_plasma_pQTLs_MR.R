# Description: correlation of MR results between CSF and plasma protein levels.

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

all_files_plasma <- list.files(here(project_dir, "results", "cis_MR_plasma", "main"), pattern = "MR_main_results_immune_pqtls_UKB*", full.names = T)
all_files_csf <- list.files(here(project_dir, "results", "cis_MR_CSF", "main"), pattern = "MR_main_results_immune_CSF_pqtls_*", full.names = T)

df_plasma <- lapply(all_files_plasma, function(x) read.table(x, header = T, sep = "\t") %>%
                      mutate(exposure = stringr::str_remove(exposure, "_[:alnum:]+_[:alnum:]+")) %>%
                        dplyr::select(outcome, exposure, nsnp, b, se, pval, FDR)) %>%
  rbindlist()

df_csf <- lapply(all_files_csf, function(x) read.table(x, header = T, sep = "\t")  %>%
                     dplyr::select(outcome, exposure, nsnp, b, se, pval, FDR)) %>%
  rbindlist()

df_all <- inner_join(df_plasma, df_csf, by = c("outcome", "exposure"), suffix = c(".plasma", ".csf")) %>%
  separate(., outcome, c("outcome", "sex", NA), sep = "_") %>%
  filter(., outcome != "AD") %>%
  filter(., outcome != "PD") %>%
  mutate(sex = case_when(
    sex == "females" ~ "Females",
    sex == "males" ~ "Males",
    sex == "meta" ~ "Sex-combined"
  ),
        outcome = case_when(
          outcome == "ADnp" ~ "AD",
          outcome == "PDnp" ~ "PD",
          outcome == "ALS" ~ "ALS"
        ),
  out_exp_sex = stringr::str_c(outcome, exposure, sex, sep = "_"))

nom_sig <- df_all %>%
  filter(., pval.plasma < 0.05 & pval.csf < 0.05)

cor.test(nom_sig$b.plasma, nom_sig$b.csf)

# Plots -----------------------------------------------------------

# beta comparisons
ggplot(df_all, aes(x = b.plasma, y = b.csf)) +
  geom_point(colour = "black", size = 2.5, shape = 1) +
  geom_point(colour = "#336699", size = 2, alpha = 6/10) +
  geom_point(data=nom_sig, size = 4.5, alpha = 10/10, aes(fill = outcome, shape = sex)) +
  theme(panel.background = element_blank(), axis.line = element_line(color = "black"), legend.key=element_blank()) + 
  labs( x = "Effect size (MR IVW with plasma protein levels)", y = "Effect size (MR IVW with CSF protein levels)") +
  scale_fill_manual("Outcome", values = c("orange", "purple", "green")) +
  scale_shape_manual("Sex", values = c("circle filled", "square filled", "triangle filled")) +
  theme(axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 18),
        legend.position = c(0.2,0.30), 
        legend.box.background = element_rect(color="black", size=0.5)) +
  guides(fill = guide_legend("Outcome", override.aes = list(shape = 21))) +
  geom_hline(yintercept=0, colour = "grey", linetype="dashed") +
  geom_vline(xintercept=0, colour = "grey", linetype="dashed") +
  geom_label_repel(data = nom_sig, aes(label = exposure), fill = "white", force = 2, size = 4.5, max.overlaps = 25)
ggsave(here(project_dir, "results", "comparison_CSF_plasma_MR_results", "beta_comparisons_updated.png"), width = 25, height = 20, units = "cm")


