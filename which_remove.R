
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
