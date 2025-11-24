# Description: generate forest plots for genome-wide MR results

# source: https://rgraphs.com/high-quality-forest-plots-in-r-ggplot2/

# Packages -------------------------------------------------------

library(here)
library(tidyverse)
library(stringr)
library(data.table)
library(grid)
library(gridExtra)
library(forcats)

# Arguments ------------------------------------------------------

project_dir = "/Users/frida/Documents/research-projects/als_immune"

cbbPalette <- c("purple", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# Read files -----------------------------------------------------

mr_results <- read.table(here(project_dir, "results", "genomewide_MR", "main", "MR_main_results_all_pairs_excl_conf_excl_HLA.txt"), sep = "\t", header = T)

# Main -----------------------------------------------------------

# note: function of TwoSampleMR "generate_odds_ratios()" can be used instead of computing them manually, as done below.

mr_results <- mr_results %>%
  mutate(OR = round(exp(b),2), low_CI = round(exp(b-1.96*se),2), upp_CI = round(exp(b+1.96*se),2), pval = round(pval,3),
         exposure_method = 
           case_when(
             method == "Inverse variance weighted (multiplicative random effects)" ~ str_c(exposure,"_IVW"),
             method == "MR Egger" ~ str_c(exposure, "_MREgger"),
             method == "Weighted median" ~ str_c(exposure, "_WM")
           ),
         exposure_name =
           case_when(
             exposure == "eosinoph_meta_sexcombined" ~ "Eosinophil count (SC)",
             exposure == "eosinophper_meta_sexcombined" ~ "Eosinophil % (SC)",
             exposure == "lymph_meta_sexcombined" ~ "Lymphocyte % (SC)",
             exposure == "mono_meta_sexcombined" ~ "Monocyte % (SC)",
             exposure == "neutroph_meta_sexcombined" ~ "Neutrophil count (SC)",
             exposure == "neutrophper_meta_sexcombined" ~ "Neutrophil % (SC)",
             exposure == "basophil_meta_sexcombined" ~ "Basophil count (SC)",
             exposure == "eosinoph_females" ~ "Eosinophil count (F)",
             exposure == "eosinoph_males" ~ "Eosinophil count (M)",
             exposure == "eosinophper_females" ~ "Eosinophil % (F)",
             exposure == "eosinophper_males" ~ "Eosinophil % (M)",
             exposure == "lymph_females" ~ "Lymphocyte % (F)",
             exposure == "lymph_males" ~ "Lymphocyte % (M)",
             exposure == "mono_females" ~ "Monocyte % (F)",
             exposure == "mono_males" ~ "Monocyte % (M)",
             exposure == "neutroph_females" ~ "Neutrophil count (F)",
             exposure == "neutroph_males" ~ "Neutrophil count (M)",
             exposure == "neutrophper_females" ~ "Neutrophil % (F)",
             exposure == "neutrophper_males" ~ "Neutrophil % (M)",
             exposure == "basophil_females" ~ "Basophil count (F)",
             exposure == "basophil_males" ~ "Basophil count (M)"
           ),
         OR_CI = str_c(OR, " (", low_CI, "-", upp_CI, ")"),
         sex = str_extract(exposure, "_[:alpha:]+") %>% str_remove(., "_"),
         sex = case_when(
           sex == "females" ~ "Female-specific",
           sex == "males" ~ "Male-specific",
           sex == "meta" ~ "Sex-combined"
         ),
         dummy_names = case_when(
           sex == "Female-specific" ~ exposure_name %>% stringr::str_remove(., "\\(F\\)"),
           sex == "Male-specific" ~ "",
           sex == "Sex-combined" ~ "")
  ) %>%
  filter(method == "Inverse variance weighted (multiplicative random effects)")

# create levels for plotting:
levels_exp_trait = arrange(mr_results, outcome, exposure) %>% .$exposure_name %>% sort() %>% unique()
mr_results <- mr_results %>%
  mutate(., exp_out_ord = factor(exposure_name, levels = levels_exp_trait))

mr_results_ALS <- mr_results %>%
  filter(outcome %in% c("ALS_females", "ALS_males", "ALS_meta_sexcombined"))

mr_results_PD <- mr_results %>%
  filter(outcome %in% c("PD_females", "PD_males", "PD_meta_sexcombined"))

mr_results_PDnp <- mr_results %>%
  filter(outcome %in% c("PDnp_females", "PDnp_males", "PDnp_meta_sexcombined"))

mr_results_AD <- mr_results %>%
  filter(outcome %in% c("AD_females", "AD_males", "AD_meta_sexcombined"))

mr_results_ADnf <- mr_results %>%
  filter(outcome %in% c("ADnp_females", "ADnp_males", "ADnp_meta_sexcombined"))

# PLOTS -------------------------------

for (outcome in c("ALS", "PD", "PDnp", "AD", "ADnp")) {
  
  if (outcome == "ALS") {
    mr_results = mr_results_ALS
    outcome = "ALS"
  } else if (outcome == "PD") {
    mr_results = mr_results_PD
    outcome = "PD"
  } else if (outcome == "PDnp") {
    mr_results = mr_results_PDnp
    outcome = "PD_exclude_UKB"
  } else if (outcome == "AD") {
    mr_results = mr_results_AD
    outcome = "AD"
  } else if (outcome == "ADnp") {
      mr_results = mr_results_ADnf
      outcome = "AD_exclude_UKB"
  } 
  
  main <- ggplot(mr_results, aes(x=OR, y=fct_rev(exp_out_ord))) +
    #Add a reference dashed line at 20
    geom_vline(xintercept = 1, linetype = "longdash", colour = "grey") + 
    #Add dot plot and error bars
    geom_errorbar(aes(xmin = low_CI, xmax = upp_CI), width = 0.25) +
    geom_point(size = 3, aes(colour = sex, shape = sex)) + 
    ggtitle("MR results (IVW method)") +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + #before = 9.6
    labs(x=stringr::str_c("OR of ", outcome), y = "Exposure") +
    scale_x_continuous(breaks=seq(floor(min(mr_results$low_CI)),ceiling(max(mr_results$upp_CI)),0.4), limits=c(floor(min(mr_results$low_CI)),ceiling(max(mr_results$upp_CI))), expand=c(0,0) ) +
    scale_shape_manual(values=c(16,15,17)) +
    theme_classic(base_size=14) +
    #Remove legend
    #Also remove y-axis line and ticks
    theme(legend.title = element_blank(),
          legend.position = c(0.8, 0.3),
          legend.box.background = element_rect(color="black", size=1),
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
                       guide = guide_legend(override.aes = list(shape = c(16,15,17), size = 3, stroke = 1.5, colour = c("#DDAA33", "#BB5566", "#004488"))))
  
  or_table <- ggplot(data=mr_results) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= OR_CI), vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + 
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
  
  pval_title <- expression(paste(italic("p"), "-value"))
  
  pval_table <- ggplot(data=mr_results) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= pval), vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + 
    ggtitle(pval_title) +
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
  
  nsnps_table <- ggplot(data=mr_results) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= nsnp), vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + 
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
  
  sex_table <- ggplot(data=mr_results) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= sex), vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + 
    ggtitle("Sex") +
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
  
  names_table <- ggplot(data=mr_results) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= dummy_names), vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=21.6, size=1) + 
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
  
  all <- grid.arrange(names_table, nsnps_table, main, or_table, pval_table, widths=c(6,2,8,4,2))
  ggsave(here(project_dir, "results", "genomewide_MR", "figures", str_c("suppl_forest_plot_",outcome,".jpg")), plot = all, width=14, height=7)
  ggsave(here(project_dir, "results", "genomewide_MR", "figures", str_c("suppl_forest_plot_",outcome,".pdf")), plot = all, width=14, height=7)
  
}
 
