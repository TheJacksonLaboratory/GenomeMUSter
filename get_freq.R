# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  calculates row frequencies of allelic state calls
# INPUT
#       x:      	dataframe or array of character entries 
# OUTPUT
#       freq_out:      array of 6 columns that contains the frequencies
#			of each genotype call for a given row (SNP)
#			A, T, C, G, H, N, where H = het and N = either NA or ""
#
get_freq <- function(x) {
  n <- ncol(x)
  input <- apply(x, 1, as.character)
  nuc <- c('A', 'T', 'C', 'G', 'H')
  freq_out <- array(NA, dim=c(nrow(x), 6), dimnames = list(NULL, c(nuc, 'N')))
  for (i in 1:length(nuc)) {
    freq_out[,i] <- colSums(input == nuc[i], na.rm = T)/n
  }
  freq_out[, 6] <- 1- rowSums(freq_out[,1:5])
  return(freq_out)
}
