#' @importFrom Rcpp evalCpp
#' @useDynLib adjVUS, .registration = TRUE
#' @import utils
#' @import np

fx <- function(s, mu, sig, error, bw1, bw2, bw3, ...){
  res <- numeric(length(s))
  res[s == 0] <- res[s == 1] <- 0
  ss <- s[s != 0 & s != 1]
  t <- log(ss/(1 - ss))
  if(missing(bw1)) u1 <- npudist(tdat = error[[1]], edat = (t - mu[1])/sig[1], ...)
  else u1 <- npudist(bws = bw1, tdat = error[[1]], edat = (t - mu[1])/sig[1], ...)
  tem1 <- u1$dist
  if(missing(bw3)) u3 <- npudist(tdat = error[[3]], edat = (t - mu[3])/sig[3], ...)
  else u3 <- npudist(bws = bw3, tdat = error[[3]], edat = (t - mu[3])/sig[3], ...)
  tem3 <- 1 - u3$dist
  if(missing(bw2)) u2 <- npudens(tdat = error[[2]], edat = (t - mu[2])/sig[2], ...)
  else u2 <- npudens(bws = bw2, tdat = error[[2]], edat = (t - mu[2])/sig[2], ...)
  tem2 <- u2$dens/sig[2]
  res[s != 0 & s != 1] <- tem1*tem2*tem3/(ss*(1 - ss))
  out <- list()
  out$res <- res
  out$bw <- c(u1$bw, u2$bw, u3$bw)
  return(out)
}

fx_y <- function(s, Y1, Y2, Y3, bw1, bw2, bw3, ...){
  res <- numeric(length(s))
  res[s == 0] <- res[s == 1] <- 0
  ss <- s[s != 0 & s != 1]
  t <- log(ss/(1 - ss))
  if(missing(bw1)) u1 <- npudist(tdat = Y1, edat = t, ...)
  else u1 <- npudist(bws = bw1, tdat = Y1, edat = t, ...)
  tem1 <- u1$dist
  if(missing(bw3)) u3 <- npudist(tdat = Y3, edat = t, ...)
  else u3 <- npudist(bws = bw3, tdat = Y3, edat = t, ...)
  tem3 <- 1 - u3$dist
  if(missing(bw2)) u2 <- npudens(tdat = Y2, edat = t, ...)
  else u2 <- npudens(bws = bw2, tdat = Y2, edat = t, ...)
  tem2 <- u2$dens
  res[s != 0 & s != 1] <- tem1*tem2*tem3/(ss*(1 - ss))
  out <- list()
  out$res <- res
  out$bw <- c(u1$bw, u2$bw, u3$bw)
  return(out)
}


integ_simpson <- function(y, h){
  n <- length(y) - 1
  s <- y[1] + y[n + 1] + 2*sum(y[seq(2, n, by = 2)]) + 4 *sum(y[seq(3, n - 1, by = 2)])
  return(s*h/3)
}
