# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  checks if two mismatched calls are nucleotide compliments
# INPUT
#	x, y	characters, allelic state calls
# OUTPUT
#		logical (TRUE if nucleotide compliments, FALSE otherwise)
#
is_compliment <- function(x,y) {
  (x=='A' & y=='T') |  (x=='T' & y=='A') |(x=='C' & y=='G') | (x=='G' & y=='C')
}
