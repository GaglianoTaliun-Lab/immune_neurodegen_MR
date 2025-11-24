# Description: generate forest plots for MR results

# source: https://rgraphs.com/high-quality-forest-plots-in-r-ggplot2/

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(grid)
library(gridExtra)
library(forcats)
library(vctrs)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read files -----------------------------------------------------

# 1) ALL MR main results files
mr_files <- list.files(here(project_dir, "results", "cis_MR_CSF", "main"), pattern = "MR_main_results_immune_CSF_pqtls*", full.names = T)
mr_results <- lapply(mr_files, function(x) read.table(x, sep = "\t", header = T) %>%
                       separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
                       mutate(exp_trait = stringr::str_c(exposure, "_", Disease))) %>%
  #      dplyr::select(exposure, outcome, nsnp, b, se, pval, FDR, MREgger_regression_pval)) %>%
  rbindlist(., use.names = TRUE, fill = TRUE)

# 2) MR Egger main results files
mr_egger_files <- list.files(here(project_dir, "results", "cis_MR_CSF", "main"), pattern = "MR_Egger_main_results_immune_CSF*", full.names = T)
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

# 3) MR significant results
# mr_sign_results <- mr_results %>%
#  filter(., FDR < 0.05) %>%
#  dplyr::select(., exposure, Disease, sex)

# 4) COLOC files
coloc_files <- setNames(list.files(here(project_dir, "results", "cis_MR_CSF", "coloc"), pattern = "significant_PPH4_0.7_*", full.names = T),
                        nm = list.files(here(project_dir, "results", "cis_MR_CSF","coloc"), pattern = "significant_PPH4_0.7_*", full.names = F) %>%
                          stringr::str_remove_all(., "significant_PPH4_0.7_") %>% stringr::str_remove_all(., "_CSF_immune_exposures.tsv")
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

coloc_files_all <- setNames(list.files(here(project_dir, "results", "cis_MR_CSF", "coloc"), pattern = "*_immune_exposures.tsv", full.names = T) %>% stringr::str_subset(., "significant_PPH4", negate = TRUE),
                        nm = list.files(here(project_dir, "results", "cis_MR_CSF", "coloc"), pattern = "*_immune_exposures.tsv", full.names = F) %>% stringr::str_subset(., "significant_PPH4", negate = TRUE) %>%
                          stringr::str_remove_all(., "_CSF_immune_exposures.tsv")
)

coloc_results_all <- lapply(coloc_files_all, function(x) read.table(x, sep = "\t", header = T)) %>%
  vctrs::list_drop_empty() %>%
  dplyr::bind_rows(., .id = "outcome") %>%
  mutate(exp_out = stringr::str_c(exposure, "_", outcome), PP.H4.abf = round(PP.H4.abf, 2)) %>%
  dplyr::select(exp_out, PP.H4.abf)

######## Generate tables for manuscript ----------------------------------------

all_table_paper <- mr_results %>%
  dplyr::select(outcome, exposure, method, nsnp, b, se, pval, FDR) %>%
  mutate(outcome = stringr::str_replace_all(outcome, "_", " ") %>% stringr::str_replace_all(., "meta sexcombined", "sex-combined"),
         outcome = case_when(
           outcome == "AD females" ~ "AD females (including UKB)",
           outcome == "ADnp females" ~ "AD females",
           outcome == "AD males" ~ "AD males (including UKB)",
           outcome == "ADnp males" ~ "AD males",
           outcome == "AD sex-combined" ~ "AD sex-combined (including UKB)",
           outcome == "ADnp sex-combined" ~ "AD sex-combined",
           outcome == "PD females" ~ "PD females (including UKB)",
           outcome == "PDnp females" ~ "PD females",
           outcome == "PD males" ~ "PD males (including UKB)",
           outcome == "PDnp males" ~ "PD males",
           outcome == "PD sex-combined" ~ "PD sex-combined (including UKB)",
           outcome == "PDnp sex-combined" ~ "PD sex-combined",
           TRUE ~ outcome
         )
  )

write.table(all_table_paper, here(project_dir, "results", "cis_MR_CSF", "main", 
                                  "cisMR_CSF_results_all_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)

mr_sign_files <- list.files(here(project_dir, "results", "cis_MR_CSF", "main"), pattern = "MR_significant_main_results_immune_CSF_pqtls_*", full.names = T)
mr_sign_results <- lapply(mr_sign_files, function(x) read.table(x, sep = "\t", header = T) %>%
                            mutate(exp_out = stringr::str_c(exposure, "_", outcome)) %>%
                            separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
                            mutate(exp_trait = stringr::str_c(exposure, "_", Disease))) %>%
  rbindlist(., use.names = TRUE, fill = TRUE)

sign_table_paper <- mr_sign_results %>%
  dplyr::select(-id.exposure, -id.outcome, -Disease, -sex, -exp_trait) %>%
  mutate(outcome = stringr::str_replace_all(outcome, "_", " ") %>% stringr::str_replace_all(., "meta sexcombined", "sex-combined"),
         outcome = case_when(
           outcome == "AD females" ~ "AD females (including UKB)",
           outcome == "ADnp females" ~ "AD females",
           outcome == "AD males" ~ "AD males (including UKB)",
           outcome == "ADnp males" ~ "AD males",
           outcome == "AD sex-combined" ~ "AD sex-combined (including UKB)",
           outcome == "ADnp sex-combined" ~ "AD sex-combined",
           outcome == "PD females" ~ "PD females (including UKB)",
           outcome == "PDnp females" ~ "PD females",
           outcome == "PD males" ~ "PD males (including UKB)",
           outcome == "PDnp males" ~ "PD males",
           outcome == "PD sex-combined" ~ "PD sex-combined (including UKB)",
           outcome == "PDnp sex-combined" ~ "PD sex-combined",
           TRUE ~ outcome
         )
  ) %>%
  dplyr::select(-exp_out)

write.table(sign_table_paper, here(project_dir, "results", "cis_MR_CSF", "main", 
                                   "cisMR_CSF_results_significant_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)

################################################################################

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

# create levels for plotting:
levels_exp_trait = arrange(mr_results2, Disease, exposure) %>% .$exp_out
mr_results2 <- mr_results2 %>%
  mutate(exp_out_ord = factor(exp_out, levels = levels_exp_trait))

# test for sexually dimorphic effects:
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
         ))

# output = no significant proteins with SDE.

# PLOTS -------------------------------

mr_results_1=mr_results2 %>%
  filter(., stringr::str_detect(Disease, "ADnp", negate = TRUE)) %>%
  filter(., stringr::str_detect(Disease, "PDnp", negate = TRUE)) %>% # set yintercept=9.6
  mutate(Tissue = "CSF",
         dummy_names = case_when(
           sex == "Female-specific" ~ Disease,
           sex == "Male-specific" ~ "",
           sex == "Sex-combined" ~ ""),
         dummy_genes = case_when(
           sex == "Female-specific" ~ exposure,
           sex == "Male-specific" ~ "",
           sex == "Sex-combined" ~ ""))

mr_results_2=mr_results2 %>%
  filter(., stringr::str_detect(outcome, "AD_", negate = TRUE)) %>%
  filter(., stringr::str_detect(outcome, "PD_", negate = TRUE)) %>%
  mutate(Disease = case_when(
    Disease == "ADnp" ~ "AD",
    Disease == "PDnp" ~ "PD"),
    Tissue = "CSF",
    dummy_names = case_when(
      sex == "Female-specific" & exposure %in% c("CD33","ADGRE2") ~ Disease,
      sex == "Male-specific" ~ "",
      sex == "Sex-combined" ~ ""),
    dummy_genes = case_when(
      sex == "Female-specific" ~ exposure,
      sex == "Male-specific" ~ "",
      sex == "Sex-combined" ~ ""))

mr_results_list <- list(mr_results_1, mr_results_2)

for (i in 1:length(mr_results_list)) {
  
  mr_results_comb <- mr_results_list[[i]]
  
  if (i == 1){
    output_suffix = "AD_PD"
    y_intercept_line = 12.6
  } else if (i == 2){
    output_suffix = "ADnp_PDnp"
    y_intercept_line = 12.6
  }
  
  main <- ggplot(mr_results_comb, aes(x=OR, y=fct_rev(exp_out_ord))) +
    #Add a reference dashed line at 20
    geom_vline(xintercept = 1, linetype = "longdash", colour = "grey") + 
    #Add dot plot and error bars
    geom_errorbar(aes(xmin = low_CI, xmax = upp_CI), width = 0.25) +
    geom_point(size = 4, aes(colour = sex, shape = sex)) + 
    ggtitle("MR results (IVW method)") +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) +
    labs(x="Odds ratio", y = "Exposure") +
    scale_x_continuous(breaks=seq(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI)),), limits=c(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI))), expand=c(0,0) ) +
    scale_shape_manual(values=c(16,15,17)) +
    theme_classic(base_size=14) +
    #Remove legend
    #Also remove y-axis line and ticks
    theme(legend.title = element_blank(),
          legend.position = c(0.83,0.81),
          # legend.box.background = element_rect(color="black", size=1),
          plot.title = element_text(hjust =0.5),
          axis.line.x = element_line(size = 0.6),
          axis.ticks.length=unit(0.3,"cm"),
          axis.text.y  = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y  = element_blank()
    ) + guides(shape = "none") +
    scale_color_manual(values = c("Female-specific" = "#DDAA33", "Male-specific" = "#BB5566", "Sex-combined" = "#004488"),
                       name = "",
                       guide = guide_legend(override.aes = list(shape = c(16,15,17), size = 4, stroke = 1.5, colour = c("#DDAA33", "#BB5566", "#004488"))))
  
  or_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= OR_CI), vjust=0, size = 6) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("OR (95% CI)") +
    xlab("  ") +
    theme_classic(base_size=14) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust =0.5) 
    ) 
  
  # Format FDR values
  format_fdr <- function(fdr_value) {
    if (fdr_value > 0.009) {
      return(sprintf("%.3f", fdr_value))
    } else {
      # Convert to scientific notation with superscript
      exponent <- floor(log10(fdr_value))
      value <- fdr_value / 10^exponent
      return(sprintf("%.2f~x~10^{%d}", value, exponent))
    }
  }
  
  # Apply formatting to create formatted FDR column
  mr_results_comb <- mr_results_comb %>%
    mutate(FDR_formatted = sapply(FDR, format_fdr))
  
  # to insert the pvalues in scientific format: label = formatC(FDR, format = "e", digits = 2)
  fdr_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= FDR_formatted), vjust=0, size = 6, parse = TRUE) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("FDR") +
    xlab("  ") +
    theme_classic(base_size=14) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust =0.5)
    )
  
  coloc_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= PP.H4.abf), vjust=0, size = 6) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("PP H4") +
    xlab("  ") +
    theme_classic(base_size=14) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust =0.5)
    )
  
  nsnps_title <- expression(paste(italic("n"), " SNPs"))
  
  nsnps_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= nsnp), vjust=0, size = 6) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle(nsnps_title) +
    xlab("  ") +
    theme_classic(base_size=14) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust =0.5)
    )
  
  disease_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= dummy_names), size = 6, vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("Disease") +
    xlab("  ") +
    theme_classic(base_size=14) +
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust =0.5)
    )
  
  names_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= dummy_genes), size = 6, vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("Exposure") +
    xlab("  ") + 
    theme_classic(base_size=14) + 
    theme(
      axis.line.y = element_blank(),
      axis.line.x = element_line(color = "white"),
      axis.text.y  = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.ticks.x = element_line(color = "white"),
      axis.ticks.length=unit(0.3,"cm"),
      axis.title.y  = element_blank(),
      axis.text.x = element_text(color="white"),
      plot.title = element_text(hjust = 0.5)
    )
  
  all <- grid.arrange(disease_table, names_table, nsnps_table, main, or_table, fdr_table, coloc_table, widths=c(2,5,3,8,4.5,3,2))
  ggsave(here(project_dir, "results", "cis_MR_CSF", "figures", stringr::str_c("forest_plot_neurodegen_CSF_pqtls_", output_suffix,"_final.jpg")), plot = all, width=15, height=7)
  ggsave(here(project_dir, "results", "cis_MR_CSF", "figures", stringr::str_c("forest_plot_neurodegen_CSF_pqtls_", output_suffix,"_final.pdf")), plot = all, width=17, height=8)
  ggsave(here(project_dir, "results", "cis_MR_CSF", "figures", stringr::str_c("forest_plot_neurodegen_CSF_pqtls_", output_suffix,"_final.tif")), plot = all, width=15, height=7, dpi = 300)

}



