# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  calculates the consensus call across all datasets
#   	ignores blanks and NAs
# INPUT
#       mm:            	dataframe of all unique genomic locations with mismatches
#	dataset:	list for each dataset used in the merging process
# OUTPUT
#	list(	called,	        dataframe of consensus calls (clear majority)
#		mismatches,	dataframe of mismatches that could not be called (ties)
#		compliments	dataframe of mismatches that were compliment calls	
# NOTES
# filenames is a vector of filenames of the datasets
# pnames are the dataset names
#
#
get_consensus <- function(mm, dataset) {
  # source functions
  source("consensus.R")
  # clean mm
  mm <- mm[,1:3]
  dups <- which(duplicated(mm))
  if (length(dups) > 0) mm <- mm[-which(duplicated(mm)),]
  mm <- mm[order(mm$strain, mm$chr, mm$bp38),]

  nfiles <- length(dataset$name)
  out <- mm
  mstrains <- unique(mm$strain)
  # for each strain, then for each dataset, grab all SNP calls
  #   identify (and return) if there is consensus SNP call across all datasets
  for (j in 1:length(mstrains)) {
    is.strain <- which(mm$strain == mstrains[j])
    x <- data.frame(bp38 = mm$bp38[is.strain])
    for (k in 1:nfiles) {
      if (mstrains[j] %in% dataset$strains[[k]]) {
        grab <-
          dataset$data[[k]][which(dataset$data[[k]]$bp38 %in% mm$bp38[is.strain]), c('bp38', mstrains[j])]
        if (nrow(grab) == 0) {
          grab <- cbind(x, rep(NA, nrow(x)))
        }
      } else {
        grab <- cbind(x, rep(NA, nrow(x)))
      }
      colnames(grab)[2] <- dataset$name[k]
      grab$strain <- mstrains[j]
      if (dataset$name[k] %in% colnames(out)) {
        grab <- grab[order(grab$bp38), ]
        if (any(duplicated(grab[,1:2]))) {
          grab <- grab[-which(duplicated(grab[,1:2])), ]
        }
        jj <-
          which(out$strain == mstrains[j] & out$bp38 %in% grab$bp38)
        if (all(out$bp38[jj] == grab$bp38)) {
          out[jj, dataset$name[k] ] <- grab[, 2]
        } else {
          print("ERROR: problem merging in consensus!")
          stop()
        }
      } else {
        out <- merge(out, grab, all.x = T, sort = F)
      }
      out <- out[order(out$strain, out$chr, out$bp38), ]
    }
  }
  out$consensus <- unlist(apply(out[, 4:ncol(out)], 1, consensus))
  mismatches <- out[which(is.na(out$consensus)), ]
  is.comp <- which(apply(mismatches, 1, check_compliments))
  comp <- NULL
  if (length(is.comp) > 0) {
    comp <- mismatches[is.comp, ]
    mismatches <- mismatches[-is.comp, ]
  }
  called <- out[which(!is.na(out$consensus)), ]
  return(list(called = called, mismatches = mismatches, compliments=comp))
}
