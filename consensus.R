# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  identifies the majority vote (consensus) across datasets for a strain and site.
# INPUT
#       votes   character vector of allelic state calls across datasets
# OUTPUT
#       c       consensus allelic state, or NA if no consensus was reached
consensus <- function(votes) {
  c <- prop <- NA
  tb <- table(unlist(votes[which(!is.na(votes) & votes!="")]))
  if (length(tb) == 0) {
    print(paste0("ERROR"))
    stop()
  }
  wins <- which(tb == max(tb))
  if (length(wins) == 1) {
    c <- names(tb)[wins]
    #prop <- round( length(which(tb==tb[wins]))/length(tb))
  }
  c
}
