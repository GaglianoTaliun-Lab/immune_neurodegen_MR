# Description: formatting of Zhang et al. (2022) pQTLs summary statistics for MR validation.

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(devtools)
library(stringr)
library(BiocManager)
library(dplyr)
library(data.table)
library(biomaRt)
library(BSgenome)
library(colochelpR)
library(rutils)
library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh37

# Choose subset of proteins for validation:
coloc_sign_files <- list.files(here(project_dir, "results", "coloc", "summary_results"), pattern = "significant_PPH4_0.7_", full.names = T) %>%
    stringr::str_subset(., "CSF_immune_exposures.tsv", negate = TRUE)

sign_proteins <- lapply(coloc_sign_files, function(x) read.table(x, sep = "\t", header = TRUE)) %>%
  rbindlist(.) %>%
  separate(., exposure, c("gene_name", "olink_id", "uniprot_id"))
  dplyr::select(., gene_symbol) %>% unique(.)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                             pQTLs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

files <- setNames(list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Zhang2022","EA"), 
                             full.names = T, pattern = "*PHENO1.glm.linear"),
                  nm = list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Zhang2022","EA"), 
                                  full.names = F, pattern = "*PHENO1.glm.linear") %>%
                    stringr::str_remove_all(., ".PHENO1.glm.linear"))

seq_ids <- read.table(here(project_dir,"gwas_sumstats","pQTLs_Zhang2022", "seqid.txt"),
                      sep = "\t", header = T) %>%
  dplyr::rename(gene_name = entrezgenesymbol)

# Calculate FDR, create a dataframe for all genes and join with gene names
pqtls_all <- lapply(files, function(x) read.table(x, sep = "\t", header= F, skip = 1, 
                    col.names = c("CHR", "BP", "ID", "REF", "ALT", "A1", "A1_FREQ", "TEST", "OBS_CT", "BETA", "SE", "T_STAT", "P", "ERRCODE")) %>%
                    dplyr::mutate(P = as.numeric(P),
                                  FDR = p.adjust(P, method = "BH")))

pqtls_all <- pqtls_all %>%
  dplyr::bind_rows(., .id = "seqid_in_sample") %>%
  left_join(., seq_ids, by = "seqid_in_sample") %>% 
  filter(., gene_name %in% sign_proteins$gene_symbol)
colnames(pqtls_all) <- c("seqid_in_sample","CHR","BP","ID","REF","ALT","A1","A1_FREQ","TEST","OBS_CT",
                         "BETA","SE","T_STAT","P","ERRCODE","FDR","uniprot_id","gene_name", "chromosome_name", "transcription_start_site")

# Keep only pQTLs with FDR < 0.05
pqtls_significant <- pqtls_all %>%
  dplyr::filter(FDR < 0.05)
  
# Format for TwoSample MR package (one dataframe per gene)
pqtls_significant <- pqtls_significant %>%
  dplyr::mutate(samplesize = 7213,
                CHR = as.factor(CHR),
                other_allele = case_when(A1 != REF ~ REF, TRUE ~ ALT)) %>%
  filter(., str_detect(ID, "rs")) %>%
  dplyr::select(SNP = ID,
                Phenotype = gene_name,
                effect_allele = A1,
                other_allele,
                eaf = A1_FREQ,
                beta = BETA,
                se = SE,
                pval = P,
                samplesize) %>%
  colochelpR::convert_rs_to_loc(df = ., "SNP", dbsnp_144) %>%
  tidyr::separate(., loc, c("CHR", "BP"), sep = ":")

pqtls_sig_list <- split(pqtls_significant, pqtls_significant$Phenotype)

for(i in 1:length(pqtls_sig_list)){
  if (nrow(pqtls_sig_list[[i]] > 0)) {
    distinct(pqtls_sig_list[[i]]$SNP, .keep_all= TRUE) %>%
    fwrite(.,
      file =
        here(project_dir, "gwas_sumstats","pQTLs_Zhang2022","pqtls_for_MR",
             stringr::str_c(names(pqtls_sig_list)[i],".tsv")), sep = "\t")
  }
}

files_to_clump <- list.files(path = here(project_dir,"gwas_sumstats","pQTLs_Zhang2022","pqtls_for_MR"), full.names = F, pattern = "*.tsv") %>%
  stringr::str_remove_all(., ".tsv") %>% as.data.frame()

write.table(files_to_clump, here(project_dir, "gwas_sumstats", "pQTLs_Zhang2022", "list_sign_genes_for_validation.txt"), sep = "\t", row.names = F, quote = F, col.names = F)

cat("There are",nrow(files_to_clump),"significant genes across all pQTLs.\n")
