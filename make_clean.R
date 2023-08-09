# fills in consensus SNP calls and returns a cleaned dataset

make_clean <- function(df=df, called=called, mismatches=mismatches) {
  cleaned <- df
  mstrains <- unique(called$strain)
  mis_strains <- unique(mismatches$strain)
  for (i in 1:length(mstrains)) {
    yy <- NULL
    xx <- called[which(called$strain == mstrains[i]), c('bp38', 'consensus')]
    is.strain <- which(colnames(cleaned) == mstrains[i])
    if (mstrains[i] %in% mis_strains) {
      yy <- mismatches[which(mismatches$strain==mstrains[i]), 'bp38']
    }
    if (length(is.strain) == 1) {
      idx <- match_bps(xx$bp38, cleaned$bp38)
      cleaned[idx, is.strain] <- xx$consensus
      if (length(yy) > 0) {
        idx <- match_bps(yy, cleaned$bp38)
        cleaned[idx, is.strain] <- ""
      }
    } else {
      print(paste0('ERROR: check strain matching in make_clean(): ', mstrains[i], ' & ', colnames(cleaned)[is.strain]))
    }
  }
  # check if we have mismatches to still remove
  missed <- setdiff(mis_strains, mstrains)
  if (length(missed) > 0) {
    for (i in 1:length(missed)) {
      is.strain <- which(colnames(cleaned) == missed[i])
      yy <- mismatches[which(mismatches$strain==missed[i]), 'bp38']
      if (length(is.strain) == 1) {
        idx <- match_bps(yy, cleaned$bp38)
        cleaned[idx, is.strain] <- ""
      } else {
        print(paste0('ERROR: check strain matching in make_clean(): ', mstrains[i], ' & ', colnames(cleaned)[is.strain]))
      }
    }
  }
  replace_N <- function(x) {
    x <- as.character(x)
    replace <- which(is.na(x) | x=="")
    if (length(replace) > 0 ) {
      x[replace] <- "N"
    }
    x
  }
  cleaned[, 7:ncol(cleaned)] <- apply(cleaned[, 7:ncol(cleaned)], 2, replace_N)
  return(cleaned)
}
