# PURPOSE
# Prepare data needed for imputation
# INPUT
# 	chr: 	string or numeric
#	start: 	numeric value for start of the region
#	end: 	numeric value for end of the region
# OUTPUT
#	if it doesn't exist, creates the phylogenetic distance matrix D
#	{indir}/distance/distance_chr{chr}_{start}-{end}.RDS
#	inputs:	list of length(number of strains) --> stores the strain-specific data needed to run haploQA HMM

prepare_for_imputation <- function(chr, start, end, strain_idx=5) {
  options(stringsAsFactors = FALSE)
  library(ape) # for phylogenetic distance
  library(data.table)
  #
  source("which_remove.R")
  source("get_D_dist.R")
  source("get_votes.R")
  source("get_missing_indices.R")
  source("get_best_strains.R")
  #
  # read in the merged datafile on the specified region
  indir <- paste0("../data/out/", chr, "/")
  infile <- paste0(indir, "merged_chr", chr, "_", start, "-", end, ".csv")
  df <- fread(infile, sep=",", skip=1, header=F, data.table=FALSE)
  cnames <- read.csv(infile, skip=0, header=F, nrows = 1)
  colnames(df) <- cnames
  #
  ddir <- paste0(indir, "distance/")
  if (!dir.exists(ddir)) {
    dir.create(ddir)
  }
  dfile <- paste0(ddir, "distance_chr", chr, "_", start, "-", end, ".RDS")
  # if phylogenetic distances on this interval are not saved, calculate it
  if (file.exists(dfile)) {
    D <- readRDS(dfile)
  } else {
    print("calculating distance matrix...")
    D <- get_D_dist(df[, strain_idx:ncol(df)])
    colnames(D) <- rownames(D) <- colnames(df)[strain_idx:ncol(df)]
    saveRDS(D, file=dfile)
    print(paste0("distance matrix save as: ", dfile))
  }
  #
  # remember to filter SNPs on which we can impute. Require at least data from 4 strains for a SNP to be imputable
  scols <- strain_idx:ncol(df)
  snp_remove <- which_remove(df, scols, nstrain = 5)
  print(paste0("There are less that 4 known genotypes for ", length(snp_remove), " SNPs. No imputation will be done for these SNPs."))
  #
  df <- df[-snp_remove, ]
  #
  is.nucleo <- function(x) {x %in% c("A","C","T","G","H")}
  # inputs is where we will store the strain-specific data needed to run haploQA HMM
  inputs <- list()
  for (i in 1:length(scols)) {
    is.sample <- scols[i]
    inputs[[i]] <- list()
    inputs[[i]]$strain <- colnames(df)[is.sample]
  }
  #
  # For each strain, determine which strains are the most phylogenetically similar (using D) and do not have too much missing data. 
  #
  # For each strain, store in d the distance between the strain and every other strain where the first entry in d is the distance to itself (0)
  # the names attached to d are the strain names in df
  impdir <- paste0(indir, "imputation_data/")
  if (!dir.exists(impdir)) {
    dir.create(impdir)
  }
  outdir <- paste0(impdir, "chr", chr, "_", start, "-", end, "/")
  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }
  subdir <- paste0(outdir, "submit/")
  if (!dir.exists(subdir)) {
     dir.create(subdir)
     system(paste0('cp submit_jobs.sh ', outdir))
  }
  shdir <- paste0(outdir, "sh/")
  if (!dir.exists(shdir)) {
     dir.create(shdir)
  }
  #
  for (i in 1:length(inputs)) {
    is.sample <- which(colnames(D) == inputs[[i]]$strain)
    others <- setdiff(which(colnames(D) %in% colnames(df)[scols]), is.sample)
    inputs[[i]]$d <- c(D[is.sample,is.sample], D[is.sample, others])
  } 
  #
  # Determine which strains are most similar (based on d) and do not have too much missing information.
  #
  # Calculate this once. 
  # For all strains, determine the indices in df (SNPs) that have missing values
  #
  # for each target strain, determine the best predictive strains based on phylogenetic distance and missingness
  # if the number of strains > 50, recommend using a distributed version
print(length(inputs))
  included_stats <- get_best_strains(df, inputs, cols=scols, outdir = outdir, strain_idx=5)
  saveRDS(included_stats, paste0(outdir, "included_stats.RDS"))
  print(paste0("summary statistics of best strain choices saved in: ", outdir, "included_stats.RDS"))
}

