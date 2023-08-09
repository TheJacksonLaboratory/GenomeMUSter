# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
# Check mismatches for compliments during the merging process
# INPUT
#       xx:        vector, individual dataset votes for genotype for a given strain at a given location
# OUTPUT
#       logical, TRUE if disagreements are compliments (e.g., A, T). FALSE if not
#
check_compliments <- function(xx) {
  source("is_compliment.R")
  votes <- as.character(xx[4:length(xx)])
  votes <- votes[votes!="" & !is.na(votes) & votes!="NA"]
  votes <- unique(votes)
  n <- length(votes)
  y <- NULL
  for (i in 1:n) {
    x <- votes[i]
    for (j in 2:n) {
      y <- as.numeric(c(y, is_compliment(x, votes[j]) | x==votes[j]))
    }
  }
  all(y==1)
}
