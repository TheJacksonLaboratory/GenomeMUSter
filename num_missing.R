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