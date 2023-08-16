# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  collects the haploqa results for each strain on a given region (chr, start, end)
# INPUT
#	chr:		integer or character representing the chromosome
#	start:		integer or character representing the start position of the region
#	end:		integer or character representing the end position of the region
#	strain_idx:	integer representing hte first column in the merged dataset that contains strain data
#				assumed that all subsequent columns contain strain data
#	testset:	logical representing if a test set was held out (default=FALSE)
# OUTPUT
#		writes the merged_and_imputed dataframe as a .csv for the region
#		writes the imputed_prob dataframe as a .csv for the region
#		writes the included_stats dataframe as a .csv of summary statistics for each strain
# NOTES
## For each strain/sample in the region, this script takes as input 
##  1. inputs.RData: the input file already created 
##  1. max_likelihood_states: the output from haplohmm.py, 
##  2. snps: the SNPs in bp38 that were included in the analysis,
##  3. snp_remove: the indices of the SNPs in df that were not used in the analysis
## REQUIREMENTS
##  R version >= 3.5
##  data.table package: this can be changed to read.csv but it is faster to use fread and fwrite
## functions
##  get_numeric_prediction.R
##  get_nucleo.R
##  match_bps.R
##  get_strain_name.R
##  is_compliment.R
##  get_ref_alt.R
##
## imputed_prob_chr{chr}_{start}-{end}.csv: a csv that contains 
##  the estimated accuracy (0, 0.999) in the region for each imputed genotype (1=observed, 0=missing)

get_imputed <- function(chr, start, end, strain_idx, testset=FALSE) {
  library(data.table)
  options(stringsAsFactors = FALSE)
  source("get_numeric_prediction.R")
  source("get_nucleo.R")
  source("match_bps.R")
  source("get_strain_name.R")
  source("is_compliment.R")
  source("get_ref_alt.R")
  is.nucleo <- function(x) {x %in% c("A","C","T","G","H")}
  #
  region <- paste0("chr", chr, "_", start, "-", end)
  indir <- paste0("../data/out/", chr, "/imputation_data/",region, "/")
  dfile <- paste0("../data/out/", chr, "/merged_", region, ".csv")
  #
  df <- fread(dfile, sep=",", data.table=FALSE)
  #
  # if we have a test set, make a imputed df to store the accuracy measures
  if (testset) {
    dfile <- paste0("../data/out/", chr, "/imputed_prob_", region, ".csv")
    #
    imp <- fread(dfile, sep=",", data.table=FALSE)
  }
  #
  # read in summary stats
  fnames <- list.files(indir, pattern = "included_stats")
  included_stats <- readRDS(paste0(indir, fnames[1]))
  if (length(fnames) > 1) {
    for (ff in 2:length(fnames)) {
      included_stats <- rbind(included_stats, readRDS(paste0(indir, fnames[ff])))
    }
    nn <- nrow(included_stats)
    included_stats <- cbind(included_stats, n_obs=rep(NA, nn),
                          n_imputed=rep(NA, nn), p_imputed=rep(NA, nn),
                         n_test=rep(NA, nn), acc=rep(NA, nn), acc_noCompliments=rep(NA, nn))
  }
  #
  nn <- nrow(included_stats)
  for (i in 1:nn) {
    print(paste0("Processing ", included_stats$strain[i]))
    dname <- paste0(indir, get_strain_name(included_stats$strain[i]), "/")
    #
    if (included_stats$strain[i]=="") { 
      dname="ERROR" 
    }
    if (!dir.exists(dname)) {
      print("ERROR: No directory exists for this strain.")
    } else {
      # read in inputs and output
      #
      snp_remove <- read.csv(paste0(dname, "snp_remove.csv"), head=T)
      snp_remove <- snp_remove[,1]
      snps <- read.csv(paste0(dname, "snps.csv"), head=T)
      snps <- snps[,1]
      obs_numeric <- read.csv(paste0(dname, "observation_ab_codes.csv"), head=T)
      obs_numeric <- obs_numeric[,1]
      input_numeric <- read.csv(paste0(dname, "haplotype_ab_codes.csv"), head=T)
      out_numeric <- read.csv(paste0(dname, "max_likelihood_states.csv"), head=F)
      dd <- readRDS(paste0(dname, "inputs.RDS"))
      #
      #
      if (length(snp_remove > 0)) {
        alleles <- dd$alleles[-snp_remove, ]
      } else {
        alleles <- dd$alleles
      }
      #
      pred_numeric <- get_numeric_prediction(out_numeric, input_numeric)
      pred_nuc <- get_nucleo(pred_numeric, alleles[, 2:5])
      # which of the imputed do we keep
      is.set <- which(obs_numeric == 0)
      to_merge <- data.frame(bp38=alleles[is.set, 1], imputed=pred_nuc[is.set] )
      is_N <- which(to_merge$imputed=="N")
      if (length(is_N) > 0 ) { 
        to_merge <- to_merge[-which(to_merge$imputed=="N"),] 
      }
      #
      #
      included_stats$n_impute[i] <- length(is.set)
      included_stats$n_snps_removed[i] <- length(snp_remove)
      included_stats$n_obs[i] <- length(obs_numeric) - length(is.set)
      included_stats$n_imputed[i] <- nrow(to_merge)
      included_stats$p_imputed[i] <- round(included_stats$n_imputed[i]/included_stats$n_impute[i]*100, 2)
      
      if (testset) {
        test <- read.csv(paste0(dname, "test.csv"))
        test2 <- merge(test, to_merge, by="bp38", all.x=F, all.y=F, sort=F)
        if (nrow(test2) > 0) {
          included_stats$n_test[i] <- nrow(test2)
          included_stats$acc[i] <- round(mean(test2$test == test2$imputed | is_compliment(test2$test, test2$imputed)), 3)
          included_stats$acc_noCompliments[i] <- round(mean(test2$test == test2$imputed), 3)
        }
      }
      # 
      # merge into the df
      is.sample <- which(colnames(df) == dd$strain)
      to_merge$df_idx <- match_bps(to_merge$bp38, df$bp38)
      to_merge$obs <- df[to_merge$df_idx, is.sample]
      to_merge <- to_merge[which(!is.nucleo(to_merge$obs)),]
      if (nrow(to_merge) > 0) {
        df[to_merge$df_idx, is.sample] <- to_merge$imputed
        imp[to_merge$df_idx, is.sample] <- min(0.999, included_stats$acc[i])  
      }
    }
  }
  #
  # save and export results
  if (testset) {
    saveRDS(included_stats, file=paste0(indir, "included_stats_testset.RDS"))
    print(paste0("Accuracy summary: ", summary(included_stats$acc, na.rm=T) ))
  } else {
    saveRDS(included_stats, file=paste0(indir, "included_stats_notestset.RDS"))
  }

  # replace any chr with just the number
  df$chr <- gsub("chr", "", df$chr)
  # fix the observed column so it matches the imputed data
  df$observed <- apply(df[, strain_idx:ncol(df)], 1, get_ref_alt)
  # fix  and save
  imp$chr <- gsub("chr", "",  imp$chr)
  imp$observed <- df$observed
  
  # remove rows that have no nucleotides
  remove <- which(rowSums(imp[, strain_idx:ncol(imp)]) == 0)
  if (length(remove) > 0) {
    df <- df[-remove, ]
    imp <- imp[-remove, ]
  }
  print(paste0(length(remove), " loci were removed due to no data. ", nrow(df), " loci remain."))
  # save it
  fwrite(df, paste0("../data/out/", chr, "/merged_and_imputed_chr", chr, "_", start, "-", end, ".csv"), 
            row.names = FALSE, sep=",")
  fwrite(imp, paste0("../data/out/", chr, "/imputed_prob_chr", chr, "_", start, "-", end, ".csv"),
            row.names = FALSE, sep=",")
}

 
