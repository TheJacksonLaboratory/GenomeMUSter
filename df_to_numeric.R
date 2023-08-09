
# Takes a dataframe, X, of snp calls (snps are rows, strains are columns)
# and a dataframe of allele calls (major in column 1, minor in column 2, third call in column 3)
# X and alleles have snp calls equal to A, C, G, T, H, "", NA, N
#
# Returns a dataframe of 0, 1, 2, 3 where the mjor allele maps to 1, 
# minor allele maps to 2, third call maps to 3, and no call maps to 0
#
#
df_to_numeric <- function(X, alleles) {
  #
  sample_to_numeric <- function(x, alleles) {
    y <- rep(0, length(x))
    y[x==alleles[,1]] <- 1
    notnas <- which(!is.na(alleles[,2]))
    if (length(notnas) > 0) {
      twos <- which(x[notnas] == alleles[notnas,2])
      y[notnas[twos]] <- 2
    }
    notnas <- which(!is.na(alleles[,3]))
    if (length(notnas) > 0) {
      threes <- which(x[notnas] == alleles[notnas, 3])
      y[notnas[threes]] <- 3
    }
    y
  }
  if (!is.null(dim(X))) {
    Y <- apply(X, 2, sample_to_numeric, alleles)
  } else { # it's a vector
    Y <- sample_to_numeric(X, alleles)
  }
}