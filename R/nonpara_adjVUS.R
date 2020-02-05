#' @importFrom Rcpp evalCpp
#' @useDynLib adjVUS, .registration = TRUE
#' @import utils
#' @import np


#' @title First nonparametric estimation method of VUS
#' @description \code{vus_np1} computes covariate-specific estimates of the VUS by using integrate with approximated distribution.
#' @param mu  a vector of three estimated conditional means for given covariate's value.
#' @param sig a vector of three estimated conditional variances for given covariate's value.
#' @param errors a list of three estimated errors in regression models.
#' @param ... additional arguments supplied to the \code{npudist()} and \code{npudens()} functions.
#'
#' @export
vus_np1 <- function(mu, sig, errors, ...){
  s <- seq(0, 1, by = 0.002)
  zs <- fx(s, mu = mu, sig = sqrt(sig), error = errors, ...)
  res <- integ_simpson(zs, 0.002)
  return(res)
}

#' @export
vus_np1_bst <- function(mu_B, sig_B, errors_B, ...){
  s <- seq(0, 1, by = 0.002)
  res <- list()
  res$t <- numeric(length(errors_B))
  for(i in 1:length(errors_B)){
    zs <- fx(s, mu = mu_B[,i], sig = sqrt(sig_B[,i]), error = errors_B[[i]], ...)
    res$t[i] <- integ_simpson(zs, 0.002)
  }
  res$t_bar <- mean(res$t)
  res$t_sd <- sd(res$t)
  return(res)
}


###-------------------------------------------------------------------------------------

#' @title Second nonparametric estimation method of VUS
#' @description \code{vus_np2} computes covariate-specific estimates of the VUS by using working samples.
#' @param mu  a vector of three estimated conditional means for given covariate's value.
#' @param sig a vector of three estimated conditional variances for given covariate's value.
#' @param errors  a list of three estimated errors in regression models.
#'
#' @export
vus_np2 <- function(mu, sig, errors){
  Y1 <- mu[1] + sqrt(sig[1])*errors[[1]]
  Y2 <- mu[2] + sqrt(sig[2])*errors[[2]]
  Y3 <- mu[3] + sqrt(sig[3])*errors[[3]]
  res <- vusC_full(Y1, Y2, Y3)
  return(res)
}

#' @export
vus_np2_bst <- function(mu_B, sig_B, errors_B){
  res <- list()
  res$t <- numeric(length(errors_B))
  for(i in 1:length(errors_B)){
    Y1 <- mu_B[1,i] + sqrt(sig_B[1,i])*errors_B[[i]][[1]]
    Y2 <- mu_B[2,i] + sqrt(sig_B[2,i])*errors_B[[i]][[2]]
    Y3 <- mu_B[3,i] + sqrt(sig_B[3,i])*errors_B[[i]][[3]]
    res$t[i] <- vusC_full(Y1, Y2, Y3)
  }
  res$t_bar <- mean(res$t)
  res$t_sd <- sd(res$t)
  return(res)
}
