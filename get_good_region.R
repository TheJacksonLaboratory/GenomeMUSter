

is.nucleo <- function(x) {x %in% c("A","C","T","G","H")}
n_is.nucleo <- function(x) mean(is.nucleo(x))
n_nucleos <- apply(df, 1, n_is.nucleo)
summary(n_nucleos)
quantile(n_nucleos, probs=seq(0, 1, 0.1))
