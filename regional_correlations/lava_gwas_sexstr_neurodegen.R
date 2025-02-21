# Description: run univariate and bivariate tests for sex-stratified neurodegenerative traits

# Load packages -----------------------------------------------------------

library(here)
library(LAVA)
library(tidyverse)
library(stringr)

# Set arguments -----------------------------------------------------------

# Argument must be 1 or 2, where 1 means to use the PD dataset that includes proxies, whereas 2 means to use the PD dataset without proxies.
args <- commandArgs(TRUE)
arg_PD_proxies <- as.numeric(args[1])

project_dir = "/home/fridald4/projects/def-gsarah/fridald4/als_immune"
out_dir <- here(project_dir,"results","lava", "gwas_sexstr")

if (arg_PD_proxies == 1) {
    
    phenotypes <- c("ALS_males","ALS_females","PD_males","PD_females", "AD_males", "AD_females")
    loc_file = here(project_dir, "reference_data", "blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile")

} else if (arg_PD_proxies == 2) {

    phenotypes <- c("ALS_males","ALS_females","PDnp_males","PDnp_females", "AD_males", "AD_females")
    loc_file = here(project_dir,	"reference_data", "blocks_s2500_EUR_chr1-23full.GRCh37_hg19.locfile")

} else {
  
  stop("Error: need to provide argument 1 (i.e., include PD dataset with proxies) or 2 (i.e., include PD dataset without proxies).")
  
}

args <- list(
    ref_prefix = here(project_dir,"reference_data","g1000_eur","g1000_eur"),
    loc_file = loc_file,
    info_file = here(project_dir,"info_files","input_info_gwas_sexstr.txt"),
    sample_overlap_file = here(project_dir, "results", "sample_overlap","sample_overlap_sexstr_neurodegen.txt"),
    phenotypes = phenotypes,
    output_filename = str_c(phenotypes, collapse = ":")
)

print(args)

# Load data ---------------------------------------------------------------

loci <- LAVA::read.loci(args$loc_file)
n_loci <- nrow(loci)
nloci_autosomes <- filter(loci, CHR != 23) %>% nrow(.)
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

for (i in 1:nloci_autosomes) {
  
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
  file = file.path(out_dir, str_c(args$output_filename, ".autosomes.univ.lava.rds"))
)
saveRDS(
  bivar,
  file = file.path(out_dir, str_c(args$output_filename, ".autosomes.bivar.lava.rds"))
)

