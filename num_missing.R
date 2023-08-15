# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#	calculated co-occurring missingness from missing_indices
# INPUT
#	index1		indices that are missing allelic state calls for the target strain
#	strain_names	character vector, potential optimal predictor strain names
#	missing_idx	list od missing indices for each strain
# OUTPUT
#	miss		vector of number of co-occurring missing between the target strain and potential predictor strains
#
# NOTES
# num_missing takes the indices of the imputed region (index1) and a character vector of strains of interest (strain_names), 
#   and the missing indices across all SNPs for all strains
# it returns the number missing in the imputed region for each strain in strain_names
num_missing <- function(index1, strain_names, missing_idx=missing_indices) {
  miss <- rep(NA, length(strain_names))
  for (i in 1:length(strain_names)) {
    miss[i] <- length(intersect(index1, missing_idx[[strain_names[i]]]))
  }
  miss
}
