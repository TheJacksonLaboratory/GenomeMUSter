# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#   uses match instead of merge to merge in new data
# INPUT
#       newdf:		vector of bps in the new dataset to be merged with x
#	x:		vector of bps
# OUTPUT
#		indices in x that match newdf
# NOTES
# returns the indices of x which match the data in newdf
# e.g., if newdf are the bps and x are the bps in df, 
#   match_bps returns the indices in df that match these bps in the correct order
#   if the bps in newdf are not in df, match_bps returns NA in that spot
match_bps <- function(newdf, x) {
  match(newdf, x)
}
