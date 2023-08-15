# AUTHOR
#  Robyn L Ball, PhD (robyn dot ball at jax dot org)
# PURPOSE
#  converts haplotype predictions to numeric allelic state calls
# INPUT
#       predicted:	dataframe od predictions
#	X:		matrix of haplotype ab calls             
# OUTPUT
#	y:		numeric prediction of allelic states
# NOTES
# Takes the predicted column and the input matrix X 
# returns the  numeric value of the prediction (0, 1, 2, 3)
#
get_numeric_prediction <- function(predicted, X){
  column_to_prediction <- function(x) {
    x[2:length(x)][x[1]]
  }
  X1 <- cbind(predicted[,1] + 1, X)
  y1 <- apply(X1, 1, column_to_prediction)
  X2 <- cbind(predicted[,2] + 1, X)
  y2 <- apply(X2, 1, column_to_prediction)
  y <- y1
  y[which(y1==0 & y2!=0)] <- y2[which(y1==0 & y2!=0)]
  #y[which(y1!=0 & y1!=y2)] <- 3
  y
}
