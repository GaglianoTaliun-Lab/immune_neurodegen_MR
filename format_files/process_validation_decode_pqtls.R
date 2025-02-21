# Description: formatting of deCODE pQTLs summary statistics for MR validation
# Note: there are two aptamers available for GPNMB (see https://menu.somalogic.com/)

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(devtools)
library(stringr)
library(BiocManager)
library(dplyr)
library(tidyr)
library(data.table)
library(biomaRt)
library(BSgenome)
library(colochelpR)
library(vctrs)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# Read file with effect allele frequency ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

eaf_decode <- fread(here(project_dir, "gwas_sumstats", "pQTLs_deCODE", "assocvariants.annotated.txt.gz")) %>%
  dplyr::select(Name, effectAllele, otherAllele, effectAlleleFreq)
cat("Loaded deCODE EAFs.\n")

pqtl_files <- list.files(here(project_dir, "gwas_sumstats", "pQTLs_deCODE"), pattern = "*.txt.gz", full.names = T) %>%
  stringr::str_subset(., "/home/fridald4/projects/def-gsarah/fridald4/als_immune/gwas_sumstats/pQTLs_deCODE/assocvariants.annotated.txt.gz", negate = T)

genes_to_test <- basename(pqtl_files) %>% 
  stringr::str_remove(., ".txt.gz") %>% 
  stringr::str_extract(., "_[:alnum:]+_[:alnum:]+") %>% 
  stringr::str_remove(., "_[:digit:]+_")

# Load gtf file
cat("Loading GTF file...\n")
gtf_file <- read.table(here(project_dir, "reference_data", "Homo_sapiens.GRCh38.v46.gtf"),
                       sep = "\t", header = F, skip = 5) %>% .[,c(1:5,9)]
colnames(gtf_file) = c("chr","source","type","start","end","info") 

gtf_file <- gtf_file %>%
  mutate(ensembl_id = stringr::str_extract(info, "ENS[:alpha:][:alnum:]+"),
    gene_name = stringr::str_extract(info, "gene_name [:graph:]+") %>%
      stringr::str_remove_all(., "gene_name ") %>%
      stringr::str_remove_all(., ";") %>%
      stringr::str_remove_all(., '"')
  ) %>%
  dplyr::select(.,
         source, type, chr, start, end, ensembl_id, gene_name) %>%
  filter(., type == "gene") %>%
  filter(., gene_name %in% genes_to_test)
cat("Done.\n")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                             Format pQTLs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (i in 1:length(genes_to_test)){

  gene=genes_to_test[i]

  if (!(gene %in% gtf_file$gene_name)) {
    cat("The gene to test is not present in the GTF file. Skipping without output.\n")
    next
  }
  cat("Processing pQTLs for ", genes_to_test[i], "...\n")

  gtf_file_subset <- gtf_file %>% filter(., gene_name == gene)

  chr = as.numeric(gtf_file_subset[1, "chr"] %>% stringr::str_remove(., "chr"))
  start = gtf_file_subset[1, "start"]
  end = gtf_file_subset[1, "end"]

  # Read file for pQTLs, filter FDR < 0.01 and extract cis-pQTLs:
  decode_gene <- fread(pqtl_files[i]) %>%
      mutate(gene_name = gene, Chrom = stringr::str_remove(Chrom, "chr"), P = Pval) %>%
      filter(., !is.na(Pval)) %>%
      mutate(FDR = p.adjust(P, method = "BH")) %>%
      filter(., FDR < 0.01) %>%
      filter(., Chrom == chr & Pos >= start-500000 & Pos <= end+500000) %>%
      filter(., !is.na(ImpMAF) | !is.na(rsids)) %>%
      dplyr::select(-effectAllele, -otherAllele) %>%
      right_join(., eaf_decode, by = "Name")

  if (nrow(decode_gene) > 0) {

    # Format for TwoSample MR package
    cat("Formatting for TwoSampleMR...\n")
    gene_MR <- decode_gene %>%
      dplyr::select(SNP = rsids,
                    CHR = Chrom,
                    BP = Pos,
                    Phenotype = gene_name,
                    effect_allele = effectAllele,
                    other_allele = otherAllele,
                    eaf = effectAlleleFreq,
                    beta = Beta,
                    se = SE,
                    pval = Pval,
                    samplesize = N) %>%
      filter(., !is.na(SNP)) %>%
      filter(., other_allele != "!") %>%
      distinct(., SNP, .keep_all= TRUE)
    
    if (nrow(gene_MR) > 0){
      fwrite(gene_MR, file = here(project_dir, "gwas_sumstats","pQTLs_deCODE","cis_pqtls", stringr::str_c("deCODE_pqtls_" ,gene, ".tsv")), sep = "\t")
      cat("Done.\n")
    } else {
      cat("There are no cis-pQTLs, not writing any output for ", gene, ".\n")
    }

  } else {
    cat("There are no cis-pQTLs, not writing any output for ", gene, ".\n")
  }
}