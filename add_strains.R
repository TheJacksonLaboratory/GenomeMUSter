# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
# Add strains (columns) to the dataframe during the merging process
# INPUT
#       addfrom:        dataframe adding from
#       addto:          dataframe adding to
#	strains:	character vector, strain names in addfrom
#	dset_strains:	character vector, strain names in addto
#	pname:		string, name of addfrom dataset (saves mismatch information)
#	outdir:		string, output directory 
#	chr:		string or numeric, chromosome
#	start:		integer, start of region
#	end:		integer, end of region
#	checkcalls	logical, default=FALSE. Set to TRUE to output disagreements of genotype between addfrom and addto
# OUTPUT
#       addto:          dataframe with new strains (columns)  added
#
#
# This function only adds strains to the dataset for locations already present.
# This is done prior to adding bps not already in the full data.
#
# If checkcalls=TRUE, this function will call check_calls() and print a file with the
# location/strain of any mismatch calls between the full dataset and the new dataset.
# It only returns mismatches that are not complementary (i.e., A/T amd C/G are not mismatches but A/G is.)
#

add_strains <- function(addfrom=input,addto=df, strains=strains, dset_strains, pname, outdir=outdir, chr=chr, start=start, end=end, checkcalls=FALSE) {
  if (checkcalls) {check_calls(addfrom, addto, strains, dset_strains, pname, outdir=outdir, chr=chr, start=start, end=end)}
  add <- setdiff(dset_strains, strains)
  dd <- addfrom[which(addfrom$bp38 %in% addto$bp38),c('bp38',add)]
  addto.row.idx <- match_bps(dd$bp38, addto$bp38)
  addto.col.idx <- ncol(addto) + 1
  dummy <- data.frame(array(NA,dim=c(nrow(addto), ncol(dd))))
  colnames(dummy) <- c('xbp38',colnames(dd)[2:ncol(dd)])
  dummy[addto.row.idx,] <- dd
  if (all(dummy[addto.row.idx,1]==addto$bp38[addto.row.idx], na.rm = TRUE)) {
    dummy <- dummy[,-1]
    addto <- cbind(addto, dummy)
    #if there's only one strain to add, treat it differently
    if (length(add) == 1) {
      colnames(addto)[addto.col.idx] <- add
    }
  } else {
    print(paste0('ERROR: BP ROWS DO NOT MATCH! CHECK IT!'))
  }
  return(addto)
}
