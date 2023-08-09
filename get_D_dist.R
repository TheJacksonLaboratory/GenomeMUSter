
# caclulates the pairwise phylogenetic distance across all strains in the dataset
# takes a dataframe X with rows as SNP and columns as strain
# returns an square array d of the phylogenetic distances between strains
# the distance is calculated using dist.dna() from the 'ape' R package
# all defaults left in place except the model for dist.dna changed to 'F84'
#
# from the documentation, 
# F84: This model generalizes K80 by relaxing the assumption of equal base frequencies. 
# It was first introduced by Felsenstein in 1984 in Phylip, 
# and is fully described by Felsenstein and Churchill (1996). 
# The formulae used in this function were taken from McGuire et al. (1999).
#
# note: this calculation seems to ignore 'H'
#
# note: If the sequences are very different, 
# most evolutionary distances are undefined 
# and a non-finite value (Inf or NaN) is returned. 
# You may do dist.dna(, model = "raw") to check whether some values are higher than 0.75.
# 

get_D_dist <- function(X) {
  require(ape)
  not_empty <- function(x) {
    which(x %in% c("A", "G", "C", "T", "H"))
  }
  n <- ncol(X)
  d <- array(NA, dim=c(n,n))
  diag(d) <- 0
  for (i in 1:n) {
    a <- X[,i]
    for (j in 1:n) {
      if (is.na(d[i,j])) {
        b <- X[,j]
        is.observed <- intersect(not_empty(a), not_empty(b))
        if (length(is.observed) > 0) {
          ab <- as.DNAbin( rbind(a[is.observed], b[is.observed]))
          d[i,j] <- d[j,i] <- dist.dna(ab, model = 'F84') 
        }
      }
    }
  }
  d
}