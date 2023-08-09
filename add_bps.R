# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
# Add genomic locations (bps) to the dataframe during the merging process
# INPUT
#       addfrom:    	dataframe adding from
#       addto:  	dataframe adding to
# OUTPUT
#	addto:		dataframe with new genomic locations (bps) added
# assumes add_strains has been called so that addto has all the strains contained in addfrom
add_bps <- function(addfrom, addto) {
  add <- which(!(addfrom$bp38 %in% addto$bp38))
  dd <- addfrom[add,]
  addto.col.idx <- match_bps(colnames(dd), colnames(addto))
  dummy <- data.frame(array(NA, dim=c(nrow(dd), ncol(addto))))
  colnames(dummy) <- colnames(addto)
  dummy[, addto.col.idx] <- dd
  if (all(colnames(addto)[addto.col.idx] == colnames(dummy)[addto.col.idx])) {
    addto <- rbind(addto, dummy)
  } else {
    print(paste0('ERROR: COLUMNS DO NOT MATCH! CHECK add_bps!'))
  }
  return(addto)
}
