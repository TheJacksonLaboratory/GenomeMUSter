# gets major, minor, and third alleles as well as probability of minor allele

get_votes <- function(X) {
  xx <- data.frame(major=rep("NA",nrow(X)), minor=rep("NA", nrow(X)), third=rep("NA", nrow(X)))
  votes <- function(x) {
    major <- minor <- third <- prob <- "NA"
    x <- as.character(x)
    tb <- table(x[x %in% c("A", "C", "G", "T", "H")])
    if (length(tb) > 0) {
      if (length(tb) == 1) {
        major <- rownames(tb)
      } else if (length(tb) > 1) {
        major <- names(tb)[which.max(tb)]
        tb <- tb[-which.max(tb)]
        minor <- names(tb)[which.max(tb)]
        tb <- tb[-which.max(tb)[1]]
        if (length(tb) > 0) {
          third <- names(tb)[1]
        }
        prob <- mean(x[x %in% c("A", "C", "G", "T", "H")] == minor) 
      }
    }
    c(major, minor, third, prob)
  }
  xx <- t(apply(X, 1, votes))
  colnames(xx) <- c("major", "minor", "third", "prob")
  return(xx)
}
