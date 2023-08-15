# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  converts numeric allelic state calls to genotype calls
# INPUT
#       x:              vector of numeric allele calls
#	alleles:	dataframe of major, minor, third allele genotypes 
# OUTPUT
#	y:		vector of genotype calls
# NOTES
# Takes the numeric input 0, 1, 2, 3 and returns the nucleotide SNP call
# 1 = major allele, 2 = minor allele, 3 = third allele, 0 = no call
# any of these can be "H"
#
get_nucleo <- function(x, alleles) {
  y <- alleles[,1]
  y[which(x==2)] <- alleles[which(x==2),2]
  y[which(x==3)] <- alleles[which(x==3),3]
  y[which(x==0)] <- "N"
  y
}
