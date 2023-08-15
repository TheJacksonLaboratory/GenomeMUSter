# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
# Prepares genotype data for HaploQA by converting genotypes to numbers
# INPUT
#       X:        	dataframe with genomic location (SNP) as rows and strains as columns
#       alleles:        dataframe of allele calls for each SNP in X
#				major in column 1, minor in column 2, third call in column 3
# OUTPUT
#       Y:            	dataframe of 0, 1, 2, 3 where the major allele maps to 1, 
# 				minor allele maps to 2, third call maps to 3, and no call maps to 0
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
