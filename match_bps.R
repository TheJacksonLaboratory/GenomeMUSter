# Function to find the indices that match for merging data efficiently (w/o calling merge)
# 
# It's just the match function. I put it here because I get confused as to which is which
# 
# returns the indices of x which match the data in newdf
# e.g., if newdf are the bps and x are the bps in df, 
#   match_bps returns the indices in df that match these bps in the correct order
#   if the bps in newdf are not in df, match_bps returns NA in that spot
match_bps <- function(newdf, x) {
  match(newdf, x)
}