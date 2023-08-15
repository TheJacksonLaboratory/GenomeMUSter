# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  for each target strain in inputs, identifies the optimal predictor strains 
#	based on co-occuring missingness and phylogenetic distance 
# INPUT
#       df:           	dataframe with SNPs as rows, strains as columns
#	inputs:		list of lists for each strain
#	cols:		integer vector of strain columns in df
#	outdir:		character, output directory
#	strain_idx:	integer, first column that is a strain in df, assumed that all columns after this column are strains
# OUTPUT
#	included_stats	dataframe containing summary statistics for each strain
#
#
get_best_strains <- function(df, inputs, cols, outdir, strain_idx=5) {
  source("num_missing.R")
  source("get_votes.R")
  source("get_strain_name.R")
  source("get_missing_indices.R")
  source("generate_haploqa_input.R")
  #
  is.nucleo <- function(x) {x %in% c("A","C","T","G","H")}
  #
  mfile <- paste0(outdir, "missing_indices.RDS")
  if (file.exists(mfile)) {
    missing_indices <- readRDS(mfile)
  } else { # compute for all of df and save it
    missing_indices <- get_missing_indices(df[, strain_idx:ncol(df)])
    saveRDS(missing_indices, file=mfile)
  }
  cnames <- colnames(df)[cols]
  n <- length(inputs)
  included_stats <-  data.frame(strain=rep("", n) , n_best=rep(NA,n), d_median=rep(NA,n), 
                                missing_median=rep(NA,n), n_impute=rep(NA,n), n_snps_removed=rep(NA, n))
  difficult <- NULL
  for (i in 1:n) {
    print(paste0("Processing ",i," of ",n))
    to_impute <- which(!is.nucleo(df[, inputs[[i]]$strain]))
    if (length(to_impute) > 0 & nrow(df) - length(to_impute) > 3) {
      # first, calculate the number missing for each strain being considered as input
      d <- inputs[[i]]$d[-1]
      # calculate the 10th percentile of the distances to other strains
      q10 <- quantile(d, probs = 0.1, na.rm = T)
      # for each potential predictor strain, calculate the number missing in the same region missing for the target strain
      missing <- num_missing(index1 = to_impute, strain_names = names(d), missing_indices)
      dd <- data.frame(strain=names(d), d, missing)
      # order first by the number missing and second by distance (see METHODS)
      dd <- dd[order(dd$missing, dd$d),]
      # we want at least 2 strains with the fewest missing
      # these must have d's that are not NA
      use <- which(!is.na(dd$d[1:20]))
      # if most of the first 20 entries have NA as the distance, this needs special handling
      if (length(use) < 2) {
        dd <- dd[which(!is.na(dd$d)),]
        use <- 1:min(5, nrow(dd))
        best <- dd[use,]
      } else {
        # best are the first 2 entries (fewest missing and non-NA d) as well as any other of the 
        #    top 20 strains that are reasonably close in phylogenetic distance
        best <- dd[union(use[1:2], which(dd$d[1:20] <= q10)),]
        # if we don't have at least 4 strains, see if relaxing the distance cutoff gives us more
        if (nrow(best) < 4) {
          q20 <- quantile(d, probs = 0.2, na.rm=T)
          best <- dd[union(use[1:2], which(dd$d[1:20] <= q20)),]
        } 
        # if we still have < 4 strains, try for one more using the 3rd lowest distance
        if (nrow(best) < 4) {
          ds <- sort(dd$d[1:20])
          best <- dd[union(use[1:2],which(dd$d[1:20] <= ds[3])),]
        }
      }
        included_stats[i, 1:4] <- c(inputs[[i]]$strain, nrow(best), round(median(best$d),3), round(median(best$missing)/length(to_impute)*100, 1))
        print(paste0("For ", inputs[[i]]$strain,": ", nrow(best), " strains were included, with d_median = ",round(median(best$d),3),"; missing_median = ",round(median(best$missing)/length(to_impute)*100, 1),"% of ",length(to_impute)))
        # store the strain names to be used for imputation 
        inputs[[i]]$best_strains <- best$strain
        # store the major and minor alleles for these best strains
        inputs[[i]]$alleles <- get_votes(df[, best$strain])
        inputs[[i]]$alleles <- cbind(bp38=df$bp38, inputs[[i]]$alleles)
        #
        # create test set and save in the correct folder
        sname <- get_strain_name(inputs[[i]]$strain)
        sdir <- paste0(outdir, sname, "/")
        if (!dir.exists(sdir)) {
          dir.create(sdir)
        }
        # export it
        saveRDS(inputs[[i]], file=paste0(sdir,"inputs.RDS"))
        print(paste0("input file saved in: ", sdir, "inputs.RDS"))
        #
        included_stats$n_snps_removed[i] <- generate_haploqa_input(df, inputs[[i]], outdir)
    } else {
      print(paste0("**************Nothing to impute for ",inputs[[i]]$strain, " or too few observed SNPs in this region."))
    }
  }
  included_stats
}
