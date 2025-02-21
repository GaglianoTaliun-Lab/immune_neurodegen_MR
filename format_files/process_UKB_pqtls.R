# Description: formatting of UKB pQTLs summary statistics across 1,248 immune-related proteins.
# Note: gene positions in olink map are in GRCh38, and therefore the summary statistics are output in GRCh38, but this does not
# affect the MR analysis since the rsids are used to match to the outcome rsids.

# Libraries ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

library(here)
library(devtools)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)

# Arguments ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#                                             Format pQTLs
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Olink map with 1,248 immune-related proteins
olink_map <- as.data.frame(fread(here(project_dir, "gwas_sumstats", "pQTLs_UKB", "olink_map_with_immune_proteins.tsv")))

for (i in 1:nrow(olink_map)){

  # get info about the protein being processed (positions are based on GRCh38):
  directory_name = olink_map[i, "file_name"] %>% stringr::str_remove(., ".tar")
  name_id = olink_map[i, "UKBPPP_ProteinID"]
  Panel = olink_map[i, "Panel"]
  chr = olink_map[i, "chr"]
  gene_start = olink_map[i, "gene_start"]
  gene_end = olink_map[i, "gene_end"]
  gene_name = olink_map[i, "HGNC.symbol"]
  olink_id = olink_map[i, "OlinkID"]
  protein_id = olink_map[i, "UniProt"]

  pheno_id = stringr::str_c(gene_name , "_", olink_id,"_",protein_id)
  file_name = stringr::str_c("discovery_chr", chr, "_", name_id, ":", Panel, ".gz")

  cat("Processing ", pheno_id, ", located in chromosome", chr, "...\n")

  # read rsid map from olink for that chromosome:
  rsid_map <- fread(stringr::str_c("/scratch/fridald4/als_immune/ukb_pqtls/olink_rsid_map_mac5_info03_b0_7_chr", chr, "_patched_v2.tsv.gz")) %>%
    dplyr::select(rsid, GENPOS = POS38)

  # Main: keep genome-wide significant cis-pQTLs, and map to rsids; format columns for TwoSample MR.
  cis_pqtls_out <- fread(stringr::str_c("/scratch/fridald4/als_immune/ukb_pqtls/", directory_name, "/", file_name)) %>%
    mutate(P = 10^-LOG10P, Phenotype = pheno_id, FDR = p.adjust(P, method = "BH")) %>%
    filter(., !is.na(P)) %>%
    dplyr::mutate(P = as.numeric(P)) %>%
    filter(., FDR < 0.01) %>%
    filter(CHROM == chr & GENPOS >= gene_start - 500000 & GENPOS <= gene_end + 500000) %>%
    left_join(., rsid_map, by = "GENPOS") %>%
    dplyr::select(SNP = rsid,
                  CHR = CHROM,
                  BP = GENPOS,
                  Phenotype,
                  effect_allele = ALLELE1,
                  other_allele = ALLELE0,
                  eaf = A1FREQ,
                  beta = BETA,
                  se = SE,
                  pval = P,
                  samplesize = N) %>%
  distinct(., SNP, .keep_all = TRUE) %>%
  filter(., !is.na(eaf)) %>%
  filter(., !is.na(SNP))

  if (nrow(cis_pqtls_out) > 0) {
      
    fwrite(cis_pqtls_out, file = here(project_dir, "gwas_sumstats","pQTLs_UKB","cis_pqtls", stringr::str_c(pheno_id, ".tsv")), sep = "\t")
    cat("Done.\n")

  } else {
    cat("There are no cis-pQTLs, not writing any output for this protein.\n")
  }

}
