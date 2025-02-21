# Description: generate forest plots for cis MR results with plasma pQTLs

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

# 3) MR significant results
# mr_sign_results <- mr_results %>%
#  filter(., FDR < 0.05) %>%
#  dplyr::select(., exposure, Disease, sex)

# 4) COLOC files
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

# 5) MRPRESSO files
mrpresso_files_all <- setNames(list.files(here(project_dir, "results", "cis_MR_plasma", "sensitivity"), pattern = "MR_PRESSO_pqtls_*", full.names = T),
                            nm = list.files(here(project_dir, "results", "cis_MR_plasma", "sensitivity"), pattern = "MR_PRESSO_pqtls", full.names = F) %>%
                              stringr::str_remove_all(., "MR_PRESSO_pqtls_") %>% stringr::str_remove(., ".txt")
)

mrpresso_results_all <- lapply(mrpresso_files_all, function(x) read.table(x, sep = "\t", header = T)) %>%
  vctrs::list_drop_empty() %>%
  dplyr::bind_rows(., .id = "exp_out") %>%
  filter(., Main.MR.results.MR.Analysis == "Raw") %>%
  dplyr::select(exp_out, MRPRESSO_Global_pval = MR.PRESSO.results.Global.Test.Pvalue)

######## Generate tables for manuscript ----------------------------------------

all_table_paper <- mr_results %>%
  dplyr::select(outcome, exposure, method, nsnp, b, se, pval, FDR) %>%
  separate(exposure, c("Exposure_gene_symbol", "Exposure_Olink_ID", "Exposure_Uniprot_ID"), extra = "drop") %>%
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

# write.table(all_table_paper, here(project_dir, "results", "cis_MR_plasma", "main", 
                              "cisMR_results_all_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)

mr_sign_files <- list.files(here(project_dir, "results", "cis_MR_plasma", "main"), pattern = "MR_significant_main_results_immune_pqtls_UKB*", full.names = T)
mr_sign_results <- lapply(mr_sign_files, function(x) read.table(x, sep = "\t", header = T) %>%
                            mutate(exp_out = stringr::str_c(exposure, "_", outcome)) %>%
                            separate(outcome, c("Disease", "sex"), sep = "_", remove = FALSE, extra = "drop") %>%
                            mutate(exp_trait = stringr::str_c(exposure, "_", Disease))) %>%
  rbindlist(., use.names = TRUE, fill = TRUE)

coloc_pp_all <- lapply(coloc_files_all, function(x) read.table(x, sep = "\t", header = T)) %>%
  vctrs::list_drop_empty() %>%
  dplyr::bind_rows(., .id = "outcome") %>%
  mutate(exp_out = stringr::str_c(exposure, "_", outcome), PP.H4.abf = round(PP.H4.abf, 2)) 

sign_table_paper <- mr_sign_results %>%
  dplyr::select(-id.exposure, -id.outcome, -Disease, -sex, -exp_trait) %>%
  separate(exposure, c("Exposure_gene_symbol", "Exposure_Olink_ID", "Exposure_Uniprot_ID"), extra = "drop", remove = F) %>%
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
  ) %>% #left_join(., coloc_pp_all, by = c("exp_out", "exposure")) %>%
  left_join(., mrpresso_results_all, by = "exp_out") %>%
  dplyr::select(-exp_out, -exposure)
  
# write.table(sign_table_paper, here(project_dir, "results", "cis_MR_plasma", "main", 
                                  "cisMR_results_significant_for_manuscript.txt"), sep = "\t", row.names = F, quote = F)

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

# add tiers according to validation status:
tiers <- fread(here(project_dir, "results", "cis_MR_plasma", "validation", "tiers_forest_plot.txt"))
mr_results2 <- left_join(mr_results2, tiers, by = c("Disease", "exposure"))

# create levels for plotting:
levels_exp_trait = arrange(mr_results2, Disease, tier, exposure) %>% .$exp_out
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

sde_test %>%
  filter(SDE == "**") %>%
  dplyr::select(-SDE) %>%
  separate(exposure, c("Exposure_gene_symbol", "Exposure_Olink_ID", "Exposure_Uniprot_ID"), extra = "drop", remove = T) %>%
  mutate(Disease = case_when(
    Disease == "AD" ~ "AD (including UKB)",
    Disease == "ADnp" ~ "AD",
    Disease == "PD" ~ "PD (including UKB)",
    Disease == "PDnp" ~ "PD",
  )) %>%
  write.table(., here(project_dir, "results", "cis_MR_plasma", "SDE_supplementary_table.txt"), sep = "\t", row.names = F, quote = F)

sde_test_join <- sde_test %>% dplyr::select(Disease, exposure, SDE)

mr_results2 <- mr_results2 %>%
  left_join(., sde_test_join, by = c("Disease", "exposure"))

# PLOTS -------------------------------

mr_results_list=mr_results2 %>%
  filter(., stringr::str_detect(outcome, "AD_", negate = TRUE)) %>%
  filter(., stringr::str_detect(outcome, "PD_", negate = TRUE)) %>%
  mutate(Disease = case_when(
    Disease == "ADnp" ~ "AD",
    Disease == "PDnp" ~ "PD",
    TRUE ~ Disease),
    Tissue = "Plasma") %>%
  separate(exposure, c("Gene_symbol", "OlinkID", "UniProtID"), sep = "_") %>%
  group_split(., Disease)

# with UKB in PD and AD:
# mr_results_list=mr_results2 %>%
#  filter(., stringr::str_detect(outcome, "ADnp_", negate = TRUE)) %>%
 # filter(., stringr::str_detect(outcome, "PDnp_", negate = TRUE)) %>%
  #mutate(Disease = case_when(
  #  Disease == "ADnp" ~ "AD",
  #  Disease == "PDnp" ~ "PD",
  #  TRUE ~ Disease),
  #  Tissue = "Plasma") %>%
  # separate(exposure, c("Gene_symbol", "OlinkID", "UniProtID"), sep = "_") %>%
  # group_split(., Disease)

for (i in 1:length(mr_results_list)) {
  
  # tiers_to_label <- array()
  # for (j in 1:3){
  #   tier_label <- filter(mr_results_list[[i]], tier == j) %>% dplyr::select(Gene_symbol) 
  #   tiers_to_label[j] <- sort(tier_label$Gene_symbol) %>% .[1]
  # }

  mr_results_comb <- mr_results_list[[i]] %>%
    mutate(dummy_names = case_when(
      sex == "Female-specific" ~ Disease,
      sex == "Male-specific" ~ "",
      sex == "Sex-combined" ~ ""),
          dummy_genes = case_when(
            sex == "Female-specific" ~ Gene_symbol,
            sex == "Male-specific" ~ "",
            sex == "Sex-combined" ~ ""),
      dummy_tiers = case_when(
        sex == "Female-specific" ~ as.character(tier),
        sex == "Male-specific" ~ "",
        sex == "Sex-combined" ~ ""),
      dummy_sde = case_when(
        sex == "Female-specific" ~ "",
        sex == "Male-specific" ~ as.character(SDE),
        sex == "Sex-combined" ~ "")
      )
  
  output_suffix = mr_results_comb$Disease[1]
  
  # parameters for output: Y intercept line in plot and height size (hgt) for figure file
  if (i == 1){ #AD 39.6 with UKB
    y_intercept_line = 27.6
    hgt=16
    wdt=18
  } else if (i == 2){ #ALS
    y_intercept_line = 3.6
    hgt=5
    wdt=18
  } else if (i == 3) { #PD 9.6 with UKB
    y_intercept_line = 12.6
    hgt=9
    wdt=18
  }
  
  main <- ggplot(mr_results_comb, aes(x=OR, y=fct_rev(exp_out_ord)), color="black") +
    #Add a reference dashed line at 20
    geom_vline(xintercept = 1, linetype = "longdash", colour = "grey") + 
    #Add dot plot and error bars
    geom_errorbar(aes(xmin = low_CI, xmax = upp_CI), width = 0.25) +
    geom_point(size = 4, aes(colour = sex)) + 
    ggtitle("MR results (IVW method)") +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) +
    labs(x=stringr::str_c("Odds ratio of ", output_suffix), y = "Exposure") +
    scale_x_continuous(breaks=seq(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI)),), limits=c(floor(min(mr_results_comb$low_CI)),ceiling(max(mr_results_comb$upp_CI))), expand=c(0,0) ) +
    scale_shape_manual(values=c(15,16,17,18)) +
    theme_classic(base_size=14) +
    #Remove legend
    #Also remove y-axis line and ticks
    theme(legend.title = element_blank(),
          legend.position = c(0.83,0.5),
          # legend.box.background = element_rect(color="black", size=1),
          plot.title = element_text(hjust =0.5),
          axis.line.x = element_line(size = 0.6),
          axis.ticks.length=unit(0.3,"cm"),
          axis.text.y  = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y  = element_blank()
    ) + guides(shape = "none")
  
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
  
  # to insert the pvalues in scientific format: label = formatC(FDR, format = "e", digits = 2)
  fdr_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= formatC(FDR, format = "e", digits = 2)), vjust=0, size = 6) +
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
  
  nsnps_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= nsnp), vjust=0, size = 6) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("N SNPs") +
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
  
  tiers_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= dummy_tiers), size = 6, vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("Tier") +
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
  
  sde_table <- ggplot(data=mr_results_comb) +
    geom_text(aes(y=fct_rev(exp_out_ord), x=1, label= dummy_sde), size = 9, vjust=0) +
    #Add a line above graph
    geom_hline(yintercept=y_intercept_line, size=2) + 
    ggtitle("SDE") +
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
  
  all <- grid.arrange(tiers_table, names_table, nsnps_table, main, or_table, fdr_table, coloc_table, sde_table, widths=c(5,5,3,8,4.5,3,2,3))
  ggsave(here(project_dir, "results", "cis_MR_plasma", "figures", stringr::str_c("forest_plot_neurodegen_plasma_pqtls_", output_suffix,"_sde.jpg")), plot = all, width=wdt, height=hgt)
  ggsave(here(project_dir, "results", "cis_MR_plasma", "figures", stringr::str_c("forest_plot_neurodegen_plasma_pqtls_", output_suffix,"_sde.pdf")), plot = all, width=wdt, height=hgt)
  
}