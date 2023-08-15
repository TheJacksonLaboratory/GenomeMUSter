# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  returns the reference, alternate, and third allelic state calls
# INPUT
#       x:              vector of allelic state calls with names(x)=strains
# OUTPUT
# 	out		character, reference/alternate/third
# NOTES
# if reference is unknown, reference = "N"
#
get_ref_alt <- function(x) {
  if ("C57BL/6J" %in% names(x)) {
    out <- as.character(x[which(names(x)=="C57BL/6J")])
    if (is.na(out) | out %in% c("", "N") ) {
      out <- "N"
    }
  } else {
    out <- "N"
  }
  x <- as.character(x)
  x[is.na(x) | x==""] <- "N"
  others <- setdiff(unique(x), out)
  others <- setdiff(others, "N")
  if (length(others) > 0) {
    if ("H" %in% others) { 
      others <- setdiff(others, "H")
      if (length(others) > 0) { out <- paste(c(out, others), collapse = "/") }
      out <- paste(c(out, "H"), collapse = "/")
    } else {
      out <- paste(c(out, others), collapse = "/")
    }
  }
  out 
}
