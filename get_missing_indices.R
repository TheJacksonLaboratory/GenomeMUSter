# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  identifies co-occurring missingness when identifying optimal predictor strains
# INPUT
#       X:              dataframe of strain data in the region (strains as columns)
# OUTPUT
#	C:		list, a vector of indices that are missing allelic state calls
#
#
get_missing_indices <- function(X) {
  is.missing <- function(x) {
    which(!(x %in% c("A", "T", "C", "G", "H")))
  }
  if (is.null(dim(X))) {
    C <- list(is.missing(X))
  } else {
    C <- apply(X, 2, is.missing) 
  }
  C
}
