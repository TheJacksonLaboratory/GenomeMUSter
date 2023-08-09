
# takes a dataframe, X, as input and returns a list where each element of 
# the list is a vector of indices that are missing values

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