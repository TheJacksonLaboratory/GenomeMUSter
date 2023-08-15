# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  calculates the proportion of agreement between a dataset allelic state calls and the consenus calls
# INPUT
#       x:            numeric vector of dataset calls
#	consensus:	numeric vector of consensus calls
# OUTPUT
#			numeric mean agreement
#
get_agreement <- function(x, consensus) {
  compare <- which(!is.na(x) & x != "")
  if (length(compare) > 0) {
    c(round(mean(x[compare] == consensus[compare]), 3), round(length(compare)))
  } else {
    c(NA, 0)
  }
}
