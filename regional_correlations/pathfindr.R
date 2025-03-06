# PathfindR: https://cran.r-project.org/web/packages/pathfindR/vignettes/intro_vignette.html

# Load libraries
library(here)
library(stringr)
library(tidyverse)
library(data.table)
library(pathfindR)

# Arguments
project_dir = "/Users/frida/Documents/research-projects/als_immune"

# Read results from LAVA
genes_files <- list.files(here(project_dir,"results","lava_neurodegen"), pattern = "genes_for_enrichment_analysis_*", full.names = T)

genes_list <- lapply(genes_files, function(x) fread(x))

geneset  = "GO-All"

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:length(genes_list)) {
  
  genes_df <- genes_list[[i]] %>%
    distinct(., .keep_all = TRUE)
    
    output_df <- run_pathfindR(genes_list[[i]], gene_sets=geneset, min_gset_size = 5, enrichment_threshold = 0.05, adj_method = "fdr")
    write.table(output_df, here(project_dir, "results", "lava_neurodegen", stringr::str_c("pathfindR_", geneset, "_enriched_terms_", i, ".txt")), sep = "\t", row.names = F, quote = F)
    enrichment_chart(result_df = output_df, top_terms = 10)
    # ggsave(here(project_dir, "results", "lava_neurodegen", stringr::str_c("enrichment_analysis_", geneset, "_geneset_", i, ".jpg")), width = 10, height = 10)
    clustered_df <- cluster_enriched_terms(output_df, use_description = TRUE)
    write.table(clustered_df, here(project_dir, "results", "lava_neurodegen", stringr::str_c("pathfindR_clustered_terms_", geneset, "_", i, ".txt")), sep = "\t", row.names = F, quote = F)
    
    output_plot <- clustered_df %>%
      filter(., Status == "Representative") %>%
      mutate(log10p = -log10(lowest_p),
             `p-values` = lowest_p)
    
    ggplot(output_plot, aes(x = log10p, y = fct_reorder(Term_Description, desc(`p-values`)), fill = `p-values`)) +
      geom_dotplot() +
      theme_classic() +
      labs(y = "GO Term Description", x = "-log10(p-value)", title = "Pathway enrichment analysis (GO)")
    
    ggsave(here(project_dir, "results", "lava_neurodegen", stringr::str_c("dotplot_clustering_", geneset, "_", i, ".jpg")), width = 8, height = 8)

}
