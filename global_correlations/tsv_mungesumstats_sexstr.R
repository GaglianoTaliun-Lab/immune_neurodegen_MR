# create tsv file for munge_sumstats from .lava.gz summary statistics
# it removes the N column (which causes an error in munge_sumstats.py script)
# it removes p-values==NA
# it removes variants without rsid

#------------------------------------------------------- Packages ------------------------------------------------------------

library(here)
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)

#------------------------------------------------------- Arguments -----------------------------------------------------------

project_dir="/home/fridald4/projects/def-gsarah/fridald4/als_immune"

lava_sumstats_files <- list.files(here(project_dir, "gwas_sumstats", "lava_sumstats", "gwas_sexstratified"), 
                                 pattern = "*.lava.gz", 
                                 recursive = T,
                                 full.names = T
)

#------------------------------------------------------- Main ----------------------------------------------------------------

# loop across all files in GWAS LAVA sumstats to output a tsv file with or without the N column
for (pheno in 1:length(lava_sumstats_files)){
  
  # cut full sumstats path to get only the GWAS name
  lava_sumstats_name <- basename(lava_sumstats_files[pheno]) %>%
    #stringr::str_remove(., here(project_dir,"gwas_sumstats","lava_sumstats","gwas","")) %>%
    stringr::str_remove(., ".lava.gz")
  
  # read lava sumstats
  pheno_sumstats <- fread(lava_sumstats_files[pheno])

  if ("A1_AF" %in% names(pheno_sumstats)) {
  
    # remove p-values = NA and top p-values to 1E-300; remove N and variant_ID columns; remove SNPs without rsid
    pheno_sumstats <- pheno_sumstats %>%
      filter(., !is.na(P)) %>%
      filter(., SNP != "") %>%
      mutate(P = as.numeric(P)) %>%
      mutate(
        P =
          case_when(
            P <= 1E-300 ~ 1E-300,
      TRUE ~ P)) %>%
      select(.,
            -N,
            -variant_ID,
            -A1_AF
      )
    
    # save sumstats as tsv files
    write.table(pheno_sumstats,
                here(project_dir,"gwas_sumstats","munge_sumstats_tsv",stringr::str_c("tmp_",lava_sumstats_name,".tsv")),
                sep = "\t",
                row.names = F,
                quote = F
    )
  }
}
