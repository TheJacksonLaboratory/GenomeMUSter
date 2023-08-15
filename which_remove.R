# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  checks for indices for which there are too few observed calls (< 5 observed allelic states)
#	these we cannot reliable impute.
# INPUT
#       dd:              dataframe with SNPs as rows and strains as columns
#	strain_cols	numeric vector if strain columns in dd
#	nstrain		integer, default=5 (require observed allelic states for at least nstrain
# OUTPUT
#		indices in dd to be removed
# NOTES
# call this once
# returns indices of SNPs in dataframe dd that should be removed 
# from genotyping imputation due to too few observations
# too few observations are any observations < nstrain
which_remove <- function(dd, strain_cols, nstrain=5) {
  sc <- strain_cols
  not_empty <- function(x) {
    sum(as.numeric(!(x == "" | is.na(x) | x=="N")))
  }
  num_values <- apply(dd[, sc], 1, not_empty)
  which(num_values < nstrain)
}
