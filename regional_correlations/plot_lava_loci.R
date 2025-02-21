# Load libraries
library(here)
library(stringr)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(rtracklayer)
library(gghighlight)

# Arguments
project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
source(here(project_dir, "scripts", "plots_RHR.R"))

# Read results from LAVA
bivar_nom1 <- read.table(here(project_dir,"results","lava","gwas_sexstr","neurodegen_sexstr_bivar_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05)
  bivar_nom2 <- read.table(here(project_dir,"results","lava","gwas_sexstr","neurodegen_PDnp_ADnp_sexstr_bivar_all.tsv"), sep = "\t", header = T) %>%
    filter(., p < 0.05)

bivar_nom <- rbind(bivar_nom1, bivar_nom2)

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Annotation Plots per LD block

ref <- rtracklayer::import(here(project_dir,"reference_data","Homo_sapiens.GRCh37.87.gtf"))
ref <- ref %>% keepSeqlevels(c(1:22), pruning.mode = "coarse") 
ref <- ref[ref$type == "gene"]

loci_gr <-
  bivar_nom %>%
  dplyr::count(locus, chr, start, stop, n_snps) %>%
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

fig_list = vector(mode = "list", length = length(loci_gr))

for(i in 1:length(loci_gr)){
  
  fig_list[[i]] <- 
    plot_locus(
    locus_gr = loci_gr[i], 
    ref = ref
    )
  
  names(fig_list)[i] <- str_c("locus_", loci_gr[i]$locus)
  
}

pdf(here(project_dir, "results", "lava", "gwas_sexstr", "figures", "phenotypes_LDblock_annotations.pdf"))
fig_list
dev.off()

