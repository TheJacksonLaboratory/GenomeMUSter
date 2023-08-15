# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  for naming directories and files by strain names, 
#	replaces all special characters with underscores
# INPUT
#       x:              character string of strain name
# OUTPUT
#	sname:		character strin of strain name with underscores instead of special characters
#
get_strain_name <- function(x) {
  sname <- gsub("\\/", "_", x)
  sname <- gsub("<", "_", sname)
  sname <- gsub(">", "_", sname)
  sname <- gsub(" ", "_", sname)
  sname <- gsub("\\.", "_", sname)
  sname <- gsub("-", "_", sname)
  sname <- gsub("\\(", "_", sname)
  sname <- gsub("\\)", "_", sname)
  sname
}
