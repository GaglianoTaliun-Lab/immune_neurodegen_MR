# Load libraries
library(here)
library(stringr)
library(tidyverse)
library(rtracklayer)
library(data.table)
library(pathfindR)

# Arguments
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
source(here(project_dir, "scripts", "plots_RHR.R"))

# Read results from LAVA
bivar_nom_all <- read.table(here(project_dir,"results","lava","gwas_sexstr","neurodegen_noUKB_sexstr_bivar_all.tsv"), sep = "\t", header = T) %>%
  filter(., p < 0.05)

bivar_ALS_PD_f <- bivar_nom_all %>%
  filter(phen1 %in% c("PDnp_females", "ALS_females")) %>%
  filter(phen2 %in% c("PDnp_females", "ALS_females"))

bivar_AD_PD_m <- bivar_nom_all %>%
  filter(phen1 %in% c("PDnp_males", "ADnp_males")) %>%
  filter(phen2 %in% c("PDnp_males", "ADnp_males"))

bivar_nom_list <- list(bivar_ALS_PD_f, bivar_AD_PD_m)

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Annotation Plots per LD block

ref <- rtracklayer::import(here(project_dir,"reference_data","Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse") 
ref <- ref[ref$type == "gene"]

for (j in 1:length(bivar_nom_list)) {

  if (j == 1) {
    out_name = "ALS_PD_females"
  } else if (j == 2) {
    out_name = "AD_PD_males"
  }

  loci_gr <-
    bivar_nom_list[[j]] %>%
    dplyr::count(locus, chr, start, stop, p, n_snps) %>%
    dplyr::arrange(locus) %>% 
    GenomicRanges::makeGRangesFromDataFrame(
      .,
      keep.extra.columns = TRUE,
      ignore.strand = TRUE,
      seqinfo = NULL,
      seqnames.field = "chr",
      start.field = "start",
      end.field = "stop"
    )

  overlap = list()

  for (i in 1:length(loci_gr)) {
    locus_gr <- loci_gr[i]
    overlap[[i]] <-
        GenomicRanges::findOverlaps(locus_gr, ref) %>%
        tibble::as_tibble() %>%
        mutate(pval = locus_gr$p)
  }

  overlap_df <- rbindlist(overlap)

  pvals <- overlap_df %>% dplyr::select(pval)

  coords <-
        tibble::tibble(
          group_type = as.factor("genes"),
          gene_id = ref[overlap_df$subjectHits]$gene_id,
          gene_name =
            ref[overlap_df$subjectHits]$gene_name %>%
              stringr::str_replace_all("-", "_") %>%
              stringr::str_replace_all("\\.", "_"),
          gene_biotype = ref[overlap_df$subjectHits]$gene_biotype %>% as.factor(),
          chr = ref[overlap_df$subjectHits] %>% GenomeInfoDb::seqnames() %>% as.character(),
          strand = ref[overlap_df$subjectHits] %>% BiocGenerics::strand() %>% as.character(),
          start = ref[overlap_df$subjectHits] %>% BiocGenerics::start(),
          end = ref[overlap_df$subjectHits] %>% BiocGenerics::end(),
          width = ref[overlap_df$subjectHits] %>% BiocGenerics::width()
        ) %>%
        cbind(., pvals) %>%
        dplyr::filter(
          gene_biotype %in% c("protein_coding")
        ) %>%
        dplyr::select(gene_name, lava_bivar_pval = pval)

  output_df <- run_pathfindR(coords, gene_sets="GO-All")
  clustered_df <- cluster_enriched_terms(output_df)

write.table(coords, here(project_dir, "results", "lava", "gwas_sexstr", stringr::str_c("genes_for_enrichment_analysis_", out_name, ".txt")), sep = "\t", row.names = F, quote = F)
write.table(output_df, here(project_dir, "results", "lava", "gwas_sexstr", stringr::str_c("pathfindR_enriched_terms_", out_name, ".txt")), sep = "\t", row.names = F, quote = F)
write.table(clustered_df, here(project_dir, "results", "lava", "gwas_sexstr", stringr::str_c("pathfindR_clustered_terms_", out_name, ".txt")), sep = "\t", row.names = F, quote = F)

}

