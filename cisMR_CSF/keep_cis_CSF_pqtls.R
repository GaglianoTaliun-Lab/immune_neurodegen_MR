# Description: script to keep cis-pQTLs from CSF data (Yang et al. 2021) for colocalisation

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

project_dir <- "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
dbsnp_144 <- SNPlocs.Hsapiens.dbSNP144.GRCh38

# Read files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

# read MR results and keep only MR significant files:
mr_sign_files <- list.files(here(project_dir, "MR", "results", "main"), pattern = "MR_main_results_immune_CSF_pqtls_*", full.names = T)
mr_sign_res <- lapply(mr_sign_files, function(x) fread(x))
mr_sign_df <- rbindlist(mr_sign_res, use.names = TRUE, fill = TRUE) %>%
  filter(., FDR < 0.05)

# Allele frequencies from 1000G EUR to add into final data. According to PLINK, Allele 1 (REF) is the minor allele
af_1kgp <- fread(here(project_dir, "reference_data", "g1000_eur", "g1000_eur.frq")) %>%
 dplyr::select(SNP, A1, A2, MAF)

# read GTF file and format it:
gtf_file <- read.table(here(project_dir, "reference_data", "Homo_sapiens.GRCh38.v46.gtf"),
                       sep = "\t", header = F, skip = 5) %>% .[,c(1:5,9)]

colnames(gtf_file) = c("chr","source","type","start","end","info")

# Main ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Get full path for each protein file:
files_to_keep <- left_join(files, seq_ids, by = "SOMAseqID") %>%
  mutate(full_name = stringr::str_c(here(project_dir,"gwas_sumstats","pQTLs_Yang2021","2021_pQTL_CSF_results_gz"),"/reducepqtls.",SOMAseqID,".glm.linear.gz"))

# Create a list of the paths and name each element with the gene name:
files_pqtls <- setNames(files_to_keep$full_name, nm = files_to_keep$gene_name)

# filter list of pQTLs to keep only if present in MR results:
sign_res <- unique(mr_sign_df$exposure)
files_pqtls <- files_pqtls[sign_res]

# read pQTLs data for MR significant proteins and format dataframes for coloc:
cis_pqtls <- lapply(files_pqtls, function(x)
  fread(x) %>%
  separate(., ID, c("CHR", "BP", "REF", "ALT"), sep = ":", remove = FALSE) %>%
  dplyr::mutate(CHR = as.factor(stringr::str_remove(CHR, "chr")), varbeta = SE^2) %>%
  colochelpR::convert_loc_to_rs(df = ., dbsnp_144) %>%
  dplyr::select(SNP,
                CHR,
                BP,
                BETA,
                SE,
                pvalues = P,
                N = OBS_CT,
                varbeta) %>%
  left_join(., af_1kgp, by = "SNP") %>%
  dplyr::select(-A1, -A2) %>%
  distinct(., SNP, .keep_all= TRUE)
  )

# format GTF file to get MR significant proteins (genes) coordinates:
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

# keep only genes within MR results:
genes_gtf <- filter(gtf_file, type == "gene" & source == "HAVANA")
genes_pqtls_gtf <- filter(genes_gtf, gene_name %in% mr_sign_df$exposure) %>%
  mutate(chr = str_remove(chr, "chr")) %>%
  dplyr::select(chr, start, end, gene_name)

# loop to keep only cis-pQTLs based on gene positions +/- 500kb and write outputs:
for (i in 1:length(cis_pqtls)) {
  
  chr = filter(genes_pqtls_gtf, gene_name == names(cis_pqtls[i])) %>% dplyr::select(chr) %>% as.numeric(.[1,1])
  start = filter(genes_pqtls_gtf, gene_name == names(cis_pqtls[i])) %>% dplyr::select(start) %>% as.numeric(.[1,1])
  end = filter(genes_pqtls_gtf, gene_name == names(cis_pqtls[i])) %>% dplyr::select(end) %>% as.numeric(.[1,1])
  
  cis_pqtls_out <- filter(cis_pqtls[[i]], CHR == chr & BP >= start-500000 & BP <= end+500000)
  
  if (nrow(cis_pqtls_out) > 0) {
    
    fwrite(cis_pqtls_out, file = here(project_dir, "gwas_sumstats","pQTLs_Yang2021","cis_pqtls_for_coloc",
                                    stringr::str_c(names(cis_pqtls)[i],".tsv")), sep = "\t")
  }
  
}



