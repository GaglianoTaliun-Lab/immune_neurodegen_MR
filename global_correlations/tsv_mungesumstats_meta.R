# create tsv file for munge_sumstats from METAL output (with CHR and BP)
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

sumstats_files <- list.files(here(project_dir,"MR","meta_analysis","metal_output"), 
                                 pattern = "*_wCHR_BP.tbl", 
                                 recursive = T,
                                 full.names = T
)

#------------------------------------------------------- Main ----------------------------------------------------------------

# loop across all files in METAL sumstats to output a tsv file with or without the N column
for (pheno in 1:length(sumstats_files)){
  
  # cut full sumstats path to get only the GWAS name
  sumstats_name <- basename(sumstats_files[pheno]) %>%
    #stringr::str_remove(., here(project_dir,"gwas_sumstats","lava_sumstats","gwas","")) %>%
    stringr::str_remove(., "_wCHR_BP.tbl")
  
  # read METAL sumstats
  pheno_sumstats <- fread(sumstats_files[pheno])

  if ("A1_AF" %in% names(pheno_sumstats)) {
  
    # remove p-values = NA and top p-values to 1E-300; remove Direction column; remove SNPs without rsid; rename columns
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
            SNP,
            CHR,
            BP,
            A1,
            A2,
            BETA,
            SE,
            P,
            -A1_AF,
            -N
      )
    
    # save sumstats as tsv files
    write.table(pheno_sumstats,
                here(project_dir,"gwas_sumstats","munge_sumstats_tsv",stringr::str_c("tmp_",sumstats_name,".tsv")),
                sep = "\t",
                row.names = F,
                quote = F
    )
  }
}
