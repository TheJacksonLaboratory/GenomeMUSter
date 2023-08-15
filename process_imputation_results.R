# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  cleans the final dataset and calculates summary statistics
#	after imputed genotypes are merged with observed data
# INPUT
#	chr		integer or character denoting the chromosome
#	strain_idx	integer, first column in the data that has strain data
#				assumed all subsequent columns contain strain data
#	testset		logical, default=FALSE, set TRUE if there was a held out test set
# OUTPUT
#	clean_merged_and_imputed .csv and clean_imputed_prob .csv for import into database
#	summary statistics over the chromosome and regions within chromosome
#
process_imputation_results <- function(chr, strain_idx, testset=FALSE) {
  source("get_ref_alt.R")
  library(data.table)
  #library(dplyr)
  
  indir <- paste0("../data/out/", chr, "/imputation_data/")
  dirs <- list.dirs(indir, recursive=F, full.names=F)

  if (testset) {
    fname <- paste0(indir, dirs[1], "/included_stats_testset.RDS")
  } else { fname <- paste0(indir, dirs[1], "/included_stats_notestset.RDS") }
  d1 <- readRDS(fname)
  d1 <- cbind(region=rep(dirs[1], nrow(d1)), d1)

  if (length(dirs) > 1) {
    for (i in 2:length(dirs)) {
      if (testset) {
        fname <- paste0(indir, dirs[i], "/included_stats_testset.RDS")
      } else { fname <- paste0(indir, dirs[i], "/included_stats_notestset.RDS") }
      d2 <- readRDS(fname)
      d2 <- cbind(region=rep(dirs[i], nrow(d2)), d2)
      d1 <- rbind(d1, d2)
    }
  }

  # check to see if we didn't have enough test cases to compute accuracy. If so, the imputed values will be removed.
  change <- which(d1$p_imputed > 0 & is.na(d1$acc))
  if (length(change) > 0) { d1[change, 9:13] <- NA }

  strains <- unique(d1$strain)
  strains <- strains[strains!=""]
  n <- length(strains)
  dd <- data.frame(strains, n_impute=rep(NA, n), n_imputed=rep(NA, n), n_not_imputed=rep(NA, n), p_imputed=rep(NA, n), completeness=rep(NA,n),
                   n_test=rep(NA, n), n_correct=rep(NA, n), n_correct_noCompliments=rep(NA, n), acc=rep(NA, n), acc_noCompliments=rep(NA, n))
  CI <- array(NA, dim=c(n,2))

  for (i in 1:n) {
    is.strain <- which(d1$strain==strains[i])
    mm <- d1[is.strain, ]
    dd[i, 2:ncol(dd)] <- c(sum(mm$n_impute, na.rm = T), # n to impute
                         sum(mm$n_imputed, na.rm = T), # n imputed
                         sum(mm$n_impute, na.rm = T) - sum(mm$n_imputed, na.rm = T), # n not imputed
                         round(sum(mm$n_imputed, na.rm = T)/ sum(mm$n_impute, na.rm = T), 3), # proportion imputed
                         round(1 -  (sum(mm$n_impute, na.rm = T) - sum(mm$n_imputed, na.rm = T))/(sum(mm$n_obs, na.rm=T) + sum(mm$n_test, na.rm=T) + sum(mm$n_imputed) ), 3), # proportion complete (have genotype)
                         sum(mm$n_test, na.rm=T), # number in test set
                         round(sum(mm$acc*mm$n_test, na.rm=T)), # number correct in test set
                         round(sum(mm$acc_noCompliments*mm$n_test, na.rm=T)), # number correct ignoring compliments
                         round(sum(mm$acc*mm$n_test, na.rm=T) / sum(mm$n_test, na.rm=T), 3), # accuracy in test set
                         round(sum(mm$acc_noCompliments*mm$n_test, na.rm=T) / sum(mm$n_test, na.rm=T), 3)) # acc ignoring compliments
    CI[i,] <- binom.test(x=dd$n_correct[i], n=dd$n_test[i])$conf.int[1:2]
  }

  dd$acc_L95 <- round(CI[,1], 3)
  dd$acc_U95 <- round(CI[,2], 3)
  write.csv(dd, paste0(indir, "imputation_summary_statistics_chr", chr, ".csv"), row.names = F)

  print("Accuracy")
  print(summary(dd$acc))

  print("Completeness")
  print(summary(dd$completeness))

  print("Accuracy < 0.5")
  print(dd$strain[which(dd$acc < 0.5)])

  print("Accuracy < 0.7")
  print(dd$strain[which(dd$acc < 0.7)])


########################################################
# generate more summary statistics and clean up the data
########################################################
  which_imputed <- function(x) {
    which(imp[, x] < 1 & imp[, x] > 0)
  }

  regions <- NULL
  counts <- NULL
  fns <- list.files(paste0("../data/out/", chr, "/"), "^merged_and_imputed_chr")
  for (fn in fns) {
    region <- gsub(".*_", "", fn)
    region <- gsub("\\.csv", "", region)
    region <- as.numeric(unlist(strsplit(region, "-")))
    regions <- rbind(regions, data.frame(start=region[1], end=region[2], fn=fn))
  }

  minstart <- min(regions$start)
  maxend <- max(regions$end)

  regions <- cbind(regions, data.frame(minstart=0, maxend=0, nSNPs=0))
  sdata <- NULL

  removed <- NULL
  for (i in 1:length(dirs)) {
    print(paste0("Processing ", dirs[i], "..."))
    fname <- paste0("../data/out/", chr, "/merged_and_imputed_", dirs[i], ".csv")
    df <- fread(fname, skip=1, header = FALSE, sep=",", data.table=F)
    cnames <- read.csv(fname, nrows=1, header=FALSE)
    fname2 <- paste0("../data/out/", chr, "/imputed_prob_", dirs[i], ".csv")
    imp <- fread(fname2, skip=1, header = FALSE, sep=",", data.table=F)
    colnames(df) <- colnames(imp) <- cnames
   
    # replace chr1 with 1
    df$chr <- gsub("chr", "", df$chr)
    imp$chr <- gsub("chr", "", imp$chr)

    obs <- as.vector(df$observed, mode="character")
    # remove rows with all "N"
    remove <- which(obs=="N")
    if (length(remove) > 0) {
      df <- df[-remove, ]
      imp <- imp[-remove, ]
      print(paste0(length(remove), " Ns removed from ", fname))
    }

    bp <- as.vector(df$bp38, mode="character")
    # remove duplicated bps
    remove <- which(duplicated(bp))
    if (length(remove) > 0) {
      df <- df[-remove, ]
      imp <- imp[-remove, ]
      print(paste0(length(remove), "dups removed from ", fname))
    }

    # for first region, remove and bps > end
    if (regions$start[i] == minstart) { # first region
      remove <- which(df$bp38 > regions$end[i])
      if (length(remove) > 0) {
        df <- df[-remove, ]
        imp <- imp[-remove, ]
        print(paste0("Removed ", length(remove), " bps not in region."))
      }
  # for last region, remove bps < start
    } else if (regions$end[i] == maxend) { # last region
      remove <- which(df$bp38 < regions$start[i])
      if (length(remove) > 0) {
        df <- df[-remove, ]
        imp <- imp[-remove, ]
        print(paste0("Removed ", length(remove), " bps not in region."))
      }
   } else {
    remove <- which(df$bp38 > regions$end[i] | df$bp38 < regions$start[i])
    if (length(remove) > 0) {
      df <- df[-remove, ]
      imp <- imp[-remove, ]
      print(paste0("Removed ", length(remove), " bps not in region."))
    }
  }
  #
  counts <- rbind(counts, data.frame(chr=chr, fn=fname, nSNPs=nrow(df)))
  #
  regions$minstart[i] <- min(df$bp38)
  regions$maxend[i] <- max(df$bp38)
  regions$nSNPs[i] <- nrow(df)
  #
  scols <- strain_idx:ncol(imp)
  for (scol in scols) {
    # check if there weren't enough test cases to compute accuracy. If so, remove these imputed values.
    if (any(is.na(imp[, scol]))) {
      replace <- which(is.na(imp[, scol]))
      imp[replace, scol] <- 0
      df[replace, scol] <- "N"
    }
    sdata <- rbind(sdata,
      data.frame(chr=chr, start=regions$start[i], end=regions$end[i], nSNPs=regions$nSNPs[i],
                 strain=colnames(imp)[scol], n_known=sum(imp[, scol]==1), n_blank=sum(imp[, scol]==0), n_imputed=sum(imp[, scol] > 0 & imp[, scol] < 1)))

  }
  # get ref and alt alleles 
  obs <- apply(df[, scols], 1, get_ref_alt)
  df$observed <- imp$observed <- obs
  #
  if (chr=="M") df$chr <- imp$chr <- "M"
  if (chr=="Y") df$chr <- imp$chr <- "Y"
  fn <- paste0("../data/out/", chr, "/clean_merged_and_imputed_", dirs[i], ".csv")
  fwrite(df, fn, row.names = FALSE)
  #
  fn <- paste0("../data/out/", chr, "/clean_imputed_prob_", dirs[i], ".csv")
  fwrite(imp, fn, row.names = FALSE)


  }


# collect across the chromosome
  strains <- unique(sdata$strain)
  cdata <- NULL
  nsnps <- sum(regions$nSNPs)

  for (strain in strains) {
    cc <- sdata[which(sdata$strain==strain), ]
    cdata <- rbind(cdata, data.frame(chr=chr, nSNPs=nsnps, strain=strain, n_known=sum(cc$n_known), n_blank=sum(cc$n_blank), n_imputed=sum(cc$n_imputed)))
  }

  write.csv(regions, paste0("../data/out/", chr, "/chr", chr, "_region_counts.csv"), row.names=FALSE)
  write.csv(sdata, paste0("../data/out/", chr, "/chr", chr, "_strain-by-region_counts.csv"), row.names=FALSE)
  write.csv(cdata, paste0("../data/out/", chr, "/chr", chr, "_strain-by-chr_counts.csv"), row.names=FALSE)

  #print(regions)
  print(paste0("Number SNPs in MUSter for chr ", chr, ": ", sum(regions$nSNP)))

}

