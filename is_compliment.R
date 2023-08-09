# function to check if 2 bp calls are compliments

is_compliment <- function(x,y) {
  (x=='A' & y=='T') |  (x=='T' & y=='A') |(x=='C' & y=='G') | (x=='G' & y=='C')
}
