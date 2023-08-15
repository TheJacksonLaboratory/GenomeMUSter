# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  checks that haplohmm.py siccessfully completed for each strain and each region in the chr
# INPUT
#       chr     integer or character denoting the chromosome
#
check_haplohmm_output <- function(chr) {
  options(stringsAsFactors = FALSE)
  indir <- paste0("../data/out/",chr, "/imputation_data/") 
  fnames <- list.dirs(indir, recursive=FALSE, full.names=FALSE)
  check <- data.frame("region", "strain")
  for (i in 1:length(fnames)) {
    print(paste0("Processing chr ", chr, " ",  fnames[i]))
    strains <- list.dirs(paste0(indir, fnames[i], "/"), recursive=FALSE, full.names=FALSE)
    strains <- strains[-which(strains %in% c("logs", "input", "sh", "submit"))]
   if (length(strains) > 0) {
   print(paste0("checking data for ", length(strains), " strains..."))
   for (j in 1:length(strains)) {
     sname <- paste0(indir, fnames[i], "/", strains[j], "/max_likelihood_states.csv")
     if (!file.exists(sname) | file.size(sname)==0) {
        write.table(paste0(strain, "/"), file=paste0(indir, fnames[i], "/config2.txt"), append=T)
        check <- rbind(check, c(fnames[i], strains[j]))
      } 
    }
   }   
  }
  if (nrow(check) > 1) {
    print(check)
    print(paste0("submit strains in ", unique(check[, 1])))
  } else {
    print("No errors. Done!")
  }
}

