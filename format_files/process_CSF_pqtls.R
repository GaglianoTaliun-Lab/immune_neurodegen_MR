# Description: formatting of pQTLs summary statistics for MR.

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
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                             Format pQTLs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Read names of files for pQTLs (SOMAseqIDs):
files <- list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Yang2021","2021_pQTL_CSF_results_gz"), 
                             full.names = F, pattern = "*.glm.linear.gz") %>%
  stringr::str_remove(., ".glm.linear.gz") %>% stringr::str_remove(., "reducepqtls.") %>% as.data.frame()
colnames(files) <- "SOMAseqID"

# Read table with uniprot IDs linked to SOMAseqIDs:
seq_ids <- read.table(here(project_dir,"gwas_sumstats","pQTLs_Yang2021", "CSF_SOMAscan1.3k_analyte_info.csv"),
                      sep = ",", header = T) %>%
  dplyr::select(gene_name = EntrezGeneSymbol,
                SOMAseqID,
                UniProt)

# Read table with only immune proteins (from Immport genes):
immune_proteins <- read.table(here(project_dir, "gwas_sumstats","pQTLs_Zhang2022", "Immport_gene_list.txt"), sep = "\t", header = F) 
colnames(immune_proteins) <- c("ensembl_id", "gene_name")
# immune_proteins <- immune_proteins %>%
#  mutate(., uniprot_id = str_remove(uniprot_id, "UniProtKB:"))

# read GTF used to get gene CHR and BP and keep only cis-pQTLs:
gtf_file <- read.table(here(project_dir, "reference_data", "Homo_sapiens.GRCh38.v46.gtf"),
                       sep = "\t", header = F, skip = 5) %>% .[,c(1:5,9)]
colnames(gtf_file) = c("chr","source","type","start","end","info")

# Filter to keep only proteins in pQTLs within the immune proteins and get full path for each file:
files_to_keep <- left_join(files, seq_ids, by = "SOMAseqID") %>%
  filter(., gene_name %in% immune_proteins$gene_name) %>%
  mutate(full_name = stringr::str_c(here(project_dir,"gwas_sumstats","pQTLs_Yang2021","2021_pQTL_CSF_results_gz"),"/reducepqtls.",SOMAseqID,".glm.linear.gz"))

# Create a list of the paths for immune proteins only and name each element with the gene name:
files <- setNames(files_to_keep$full_name, nm = files_to_keep$gene_name)

# Allele frequencies from 1000G EUR to add into final data. According to PLINK, Allele 1 (REF) is the minor allele
af_1kgp <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_eur.frq")) %>%
 dplyr::select(SNP, A1, A2, MAF)

# Calculate FDR, create a dataframe for all genes and join with gene names
cat("Calculating FDR, filtering to keep FDR < 0.05, and creating a dataframe for all genes...\n")
pqtls_all <- lapply(files, function(x) fread(x) %>%
                    filter(., !is.na(P)) %>%
                    dplyr::mutate(P = as.numeric(P), FDR = p.adjust(P, method = "BH")) %>%
                    filter(., FDR < 0.05)
)
pqtls_all <- list_drop_empty(pqtls_all)
cat("Done.\n")

cat("List to dataframe and back to add column with gene name, and to remove empty elements in list...\n")
for (n in names(pqtls_all)) {
  pqtls_all[[n]][['gene_name']] = n
}
pqtls_significant_df = rbindlist(pqtls_all, fill = TRUE)

cat("Across all proteins, there are a total of ", nrow(pqtls_significant_df), " significant QTLs.\n")
pqtls_sig_list <- split(pqtls_significant_df, pqtls_significant_df$gene_name)
cat("Done.\n")

# Format for TwoSample MR package (one dataframe per gene)
# no data for allele frequencies - add from 1KGP-eur
cat("Formatting for TwoSampleMR...\n")
pqtls_list_MR <- lapply(pqtls_sig_list, function(x)
  separate(x, ID, c("CHR", "BP", "REF", "ALT"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(CHR = as.factor(stringr::str_remove(CHR, "chr"))) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144) %>%
  dplyr::select(SNP,
                CHR,
                BP,
                Phenotype = gene_name,
                effect_allele = ALT,
                other_allele = REF,
                beta = BETA,
                se = SE,
                pval = P,
                samplesize = OBS_CT) %>%
  left_join(., af_1kgp, by = "SNP") %>%
  mutate(eaf = case_when(
    A1 == effect_allele ~ MAF,
    A2 == effect_allele ~ 1-MAF    
  )) %>%
  dplyr::select(-A1, -A2, -MAF) %>%
  distinct(., SNP, .keep_all= TRUE)
)
cat("Done.\n")

# Keep only cis-pQTLs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cat("Filter to keep only cis-pQTLs, and save each file...\n")
gtf_file <- gtf_file %>%
  mutate(
    ensembl_id = stringr::str_extract(info, "ENS[:alpha:][:alnum:]+"),
    gene_name = stringr::str_extract(info, "gene_name [:graph:]+") %>%
      stringr::str_remove_all(., "gene_name ") %>%
      stringr::str_remove_all(., ";") %>%
      stringr::str_remove_all(., '"')
  ) %>%
  dplyr::select(.,
         source, type, chr, start, end, ensembl_id, gene_name)

# keep only genes in the pQTL list for MR:
genes_gtf <- filter(gtf_file, type == "gene" & source == "HAVANA")
genes_pqtls_gtf <- filter(genes_gtf, gene_name %in% names(pqtls_list_MR)) %>%
  mutate(chr = str_remove(chr, "chr")) %>%
  dplyr::select(chr, start, end, gene_name)

# Extract cis-pQTLs for each gene and write a table with cis-pQTLs:
for (i in 1:length(pqtls_list_MR)) {
  
  chr = filter(genes_pqtls_gtf, gene_name == names(pqtls_list_MR[i])) %>% dplyr::select(chr) %>% as.numeric(.[1,1])
  start = filter(genes_pqtls_gtf, gene_name == names(pqtls_list_MR[i])) %>% dplyr::select(start) %>% as.numeric(.[1,1])
  end = filter(genes_pqtls_gtf, gene_name == names(pqtls_list_MR[i])) %>% dplyr::select(end) %>% as.numeric(.[1,1])
  
  cis_pqtls_out <- filter(pqtls_list_MR[[i]], CHR == chr & BP >= start-500000 & BP <= end+500000) %>%
    filter(., !is.na(eaf) | !is.na(SNP))
  
  if (nrow(cis_pqtls_out) > 0) {
    
    fwrite(cis_pqtls_out, file = here(project_dir, "gwas_sumstats","pQTLs_Yang2021","cis_pqtls",
                                    stringr::str_c(names(pqtls_list_MR)[i],".tsv")), sep = "\t")
  }
}
cat("Done.\n")

#################################

files_to_clump <- list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Yang2021","cis_pqtls"), full.names = F, pattern = "*.tsv") %>%
  stringr::str_remove_all(., ".tsv") %>% as.data.frame()

write.table(files_to_clump, here(project_dir, "gwas_sumstats", "pQTLs_Yang2021", "list_sign_genes.txt"), sep = "\t", row.names = F, quote = F, col.names = F)

cat("There are",nrow(files_to_clump),"significant genes across all CSF pQTLs.\n")
