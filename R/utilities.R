#' @importFrom Rcpp evalCpp
#' @useDynLib adjVUS, .registration = TRUE
#' @import utils
#' @import np

fx <- function(s, mu, sig, error, ...){
  res <- numeric(length(s))
  res[s == 0] <- res[s == 1] <- 0
  ss <- s[s != 0 & s != 1]
  t <- log(ss/(1 - ss))
  tem1 <- npudist(tdat = error[[1]], edat = (t - mu[1])/sig[1], ...)$dist
  tem3 <- 1 - npudist(tdat = error[[3]], edat = (t - mu[3])/sig[3], ...)$dist
  tem2 <- npudens(tdat = error[[2]], edat = (t - mu[2])/sig[2], ...)$dens/sig[2]
  res[s != 0 & s != 1] <- tem1*tem2*tem3/(ss*(1 - ss))
  return(res)
}

integ_simpson <- function(y, h){
  n <- length(y) - 1
  s <- y[1] + y[n + 1] + 2*sum(y[seq(2, n, by=2)]) + 4 *sum(y[seq(3, n - 1, by = 2)])
  return(s*h/3)
}
