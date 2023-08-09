
# replaces all special characters with underscores
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
