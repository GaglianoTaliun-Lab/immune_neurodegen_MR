# Description: run univariate and bivariate tests for sex-stratified neurodegenerative traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
out_dir <- here(project_dir,"results","lava", "gwas_sexstr")
    
phenotypes <- c("ALSX_males","ALSX_females","PDX_males","PDX_females", "ADX_males", "ADX_females")
loc_file = here(project_dir, "reference_data", "blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile")

args <- list(
    ref_prefix = here(project_dir,"reference_data","g1000_eur","g1000_eur"),
    loc_file = loc_file,
    info_file = here(project_dir,"info_files","input_info_gwas_sexstr.txt"),
    sample_overlap_file = NULL,
    phenotypes = phenotypes,
    output_filename = str_c(phenotypes, collapse = ":")
)

print(args)

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)
input <-
  LAVA::process.input(
    input.info.file = args$info_file,
    sample.overlap.file = args$sample_overlap_file,
    ref.prefix = args$ref_prefix,
    phenos = args$phenotypes
  )

# Main --------------------------------------------------------------------

# Print progress
print(str_c("Starting LAVA analysis for ", n_loci, " loci"))
progress <-
  quantile(
    x = 1:n_loci,
    probs = seq(.05,1,.05)
  ) %>%
  ceiling()

# Set univariate threshold to 0.05/n_loci
univar_threshold <-
  0.05/n_loci

univar = bivar = list()

for (i in 2496:2589) {  # start and end loci number for chromosome X
  
  if (i %in% progress) print(str_c("..", names(progress[which(progress==i)])))     # (printing progress)
  
  # Process locus
  locus <-
    LAVA::process.locus(
      loci[i,],
      input
    )
  
  # It is possible that the locus cannot be defined for various reasons (e.g. too few SNPs),
  # The !is.null(locus) check is necessary before calling the analysis functions.
  if (!is.null(locus)) {
    
    # extract some general locus info for the output
    loc_info <-
      data.frame(
        locus = locus$id,
        chr = locus$chr,
        start = locus$start,
        stop = locus$stop,
        n_snps = locus$n.snps,
        n_pcs = locus$K
      )
    
    # Run the univariate and bivariate tests
    loc_out <-
      LAVA::run.univ.bivar(
        locus,
        univ.thresh = univar_threshold
      )
    
    # Bind
    univar[[i]] <-
      loc_info %>%
      dplyr::bind_cols(loc_out$univ)
    
    if(!is.null(loc_out$bivar)){
      
      bivar[[i]] <-
        loc_info %>%
        dplyr::bind_cols(loc_out$bivar)
      
    }
    
  }
  
}

# Save data ---------------------------------------------------------------

saveRDS(
  univar,
  file = file.path(out_dir, str_c(args$output_filename, ".chrx.univ.lava.rds"))
)
saveRDS(
  bivar,
  file = file.path(out_dir, str_c(args$output_filename, ".chrx.bivar.lava.rds"))
)