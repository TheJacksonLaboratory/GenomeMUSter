# votes is a vector of the votes across all datasets
# returns the majority vote or returns NA if no consensus could be reached
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
