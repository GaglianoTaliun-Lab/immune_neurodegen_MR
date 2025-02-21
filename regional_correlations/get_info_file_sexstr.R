# Get file info file for LAVA  - only for neurodegenerative traits sex-stratified (ALS, AD, PD)

# Load packages -----------------------------------------------------------

library(dplyr)
library(here)
library(stringr)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"

# Read file ---------------------------------------------------------------

N_samples <- read.table(here(project_dir,"gwas_sample_sizes_sexstr_neurodegen.tsv"), sep="\t", header = T) %>% 
  arrange(., trait, by_group=F) %>%
  dplyr::select(
    phenotype = trait,
    cases = n_cases,
    controls = n_controls) %>%
  filter(., stringr::str_detect(phenotype, "sexcombined", negate = TRUE))

# Main -------------------------------------------------------------------

file_path <- data.frame(filename = list.files(
  path = here(project_dir,"gwas_sumstats","lava_sumstats","gwas_sexstratified"),
  full.names = T,
  recursive = T,
  pattern = ".lava.gz"
))

file_name <- data.frame(basename = list.files(
  path = here(project_dir,"gwas_sumstats","lava_sumstats","gwas_sexstratified"),
  full.names = F,
  recursive = T,
  pattern = ".lava.gz"
)) %>% mutate(basename = stringr::str_remove(basename, ".lava.gz"))

files <- cbind(file_path, file_name)

input_info <- left_join(N_samples, files, by = c("phenotype" = "basename"))

write.table(
  input_info,
  here(project_dir,"info_files","input_info_gwas_sexstr.txt"),
  sep = "\t", row.names = F, quote = F
)