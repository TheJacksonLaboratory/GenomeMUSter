# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
# Check for mismatches in genotype calls  during the merging process
# INPUT
#       addfrom:        dataframe adding from
#       addto:          dataframe adding to
#       strains:        character vector, strain names in addfrom
#       dset_strains:   character vector, strain names in addto
#       pname:          string, name of addfrom dataset (saves mismatch information)
#       outdir:         string, output directory 
#       chr:            string or numeric, chromosome
#       start:          integer, start of region
#       end:            integer, end of region
# OUTPUT
#       csv:          	writes csv file to {outdir}/mismatches/ with genotype mismatch information between addto and addfrom
#
#
# assumes add_strains has been called so that addto has all the strains contained in addfrom
# If checkcalls=TRUE in add_strains(),  check_calls() outputs a file with the location/strain
# of any mismatch calls between the full dataset and the new dataset.
# It only returns mismatches that are not complementary (i.e., A/T amd C/G are not mismatches but A/G is.)

check_calls <- function(addfrom, addto,  strains, dset_strains, pname,
                        outdir=outdir, chr=chr, start=start, end=end) {
  fn <- paste0(outdir, 'mismatches/mismatches_chr', chr, '_', start, '-', end, '_', pname, '.csv')
  print(paste0('          Mismatch information in: ',fn))
  out <- out2 <- data.frame(strain="", chr=NA, bp38=NA, call="")
  dd <- addfrom
  scheck <- intersect(strains, dset_strains)
  addto.row.idx <- match_bps(dd$bp38, addto$bp38)
  ignore_addfrom <- c(1,which(is.na(addto.row.idx)))
  addto.row.idx <- addto.row.idx[-ignore_addfrom]
  addfrom.row.idx <- c(1:nrow(dd))[-ignore_addfrom]
  nmismatches <- 0
  for (i in 1:length(scheck)) {
    x <- addto[addto.row.idx, scheck[i]]
    y <- dd[addfrom.row.idx, scheck[i]]
    # include compliments as a mismatch to start so we can see if we can find consensus
    #  if we can't find consensus, ignore compliments and do not call them mismatches
    is.mismatch <- which(!((x == y) | y=="") )
    nmis <- length(is.mismatch)
    nmismatches <- nmismatches + nmis
    if (nmis > 0) {
      dd.idx <- addfrom.row.idx[is.mismatch]
      mm <- data.frame(strain=rep(scheck[i],nmis), chr=dd$chr[dd.idx],
                       bp38=dd$bp38[dd.idx], call=dd[dd.idx, scheck[i]])
      out <- rbind(out,mm)
    }
  }
  out <- out[-1,];
  write.csv(out, file=fn, quote = FALSE, row.names = FALSE)
}
