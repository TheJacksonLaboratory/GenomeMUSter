

# function that calculates the proportion agreement among a dataset call and the consensus call
get_agreement <- function(x, consensus) {
  compare <- which(!is.na(x) & x != "")
  if (length(compare) > 0) {
    c(round(mean(x[compare] == consensus[compare]), 3), round(length(compare)))
  } else {
    c(NA, 0)
  }
}
