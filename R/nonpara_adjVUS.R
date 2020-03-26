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
  u <- fx(s, mu = mu, sig = sqrt(sig), error = errors, ...)
  zs <- u$res
  res <- integ_simpson(zs, 0.002)
  out <- list()
  out$vus <- res
  out$bws <- u$bw
  return(out)
}

#' @export
vus_np1_fy <- function(Y1, Y2, Y3, ...){
  s <- seq(0, 1, by = 0.002)
  res <- numeric(ncol(Y1))
  bws <- matrix(0, nrow = 3, ncol = ncol(Y1))
  for(i in 1:ncol(Y1)){
    u <- fx_y(s, Y1[,i], Y2[,i], Y3[,i], ...)
    zs <- u$res
    res[i] <- integ_simpson(zs, 0.002)
    bws[,i] <- u$bw
  }
  out <- list()
  out$vus <- res
  out$bws <- bws
  return(out)
}


#' @export
vus_np1_bst <- function(mu_B, sig_B, errors_B, bw1, bw2, bw3, ...){
  s <- seq(0, 1, by = 0.002)
  res <- list()
  res$t <- numeric(length(errors_B))
  for(i in 1:length(errors_B)){
    temp <- try(fx(s, mu = mu_B[,i], sig = sqrt(sig_B[,i]), error = errors_B[[i]], bw1, bw2, bw3, ...), silent = TRUE)
    if(inherits(temp, "try-error")) res$t[i] <- NA
    else {
      zs <- temp$res
      res$t[i] <- integ_simpson(zs, 0.002)
    }
  }
  res$t_bar <- mean(res$t, na.rm = TRUE)
  res$t_sd <- sd(res$t, na.rm = TRUE)
  return(res)
}

#' @export
vus_np1_fy_bst <- function(Y1_B, Y2_B, Y3_B, bw1, bw2, bw3, ci.level = 0.95, ...){
  s <- seq(0, 1, by = 0.002)
  R <- length(Y1_B)
  nx <- ncol(Y1_B[[1]])
  res <- list()
  res$t <- matrix(0, nrow = nx, ncol = R)
  for(i in 1:R){
    for(j in 1:nx){
      temp <- try(fx_y(s, Y1_B[[i]][,j], Y2_B[[i]][,j], Y3_B[[i]][,j], bw1[j], bw2[j], bw3[j], ...), silent = TRUE)
      if(inherits(temp, "try-error")) res$t[j,i] <- NA
      else {
        zs <- temp$res
        res$t[j,i] <- integ_simpson(zs, 0.002)
      }
    }
  }
  res$t_bar <- rowMeans(res$t, na.rm = TRUE)
  res$t_sd <- apply(res$t, 1, sd, na.rm = TRUE)
  res$pci <- apply(res$t, 1, function(x){
    temp <- sort(x)
    res <- temp[c(floor(length(temp)*(1 - ci.level)/2), floor(length(temp)*(1 - (1 - ci.level)/2)))]
  })
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
vus_np2 <- function(Y1, Y2, Y3){
  res <- numeric(ncol(Y1))
  for(i in 1:ncol(Y1)){
    res[i] <- vus(Y1[,i], Y2[,i], Y3[,i])
  }
  return(res)
}

#' @export
vus_np2_bst <- function(Y1_B, Y2_B, Y3_B, ci.level = 0.95){
  R <- length(Y1_B)
  nx <- ncol(Y1_B[[1]])
  res <- list()
  res$t <- matrix(0, nrow = nx, ncol = R)
  for(i in 1:R){
    for(j in 1:nx){
      res$t[j,i] <- vus(Y1_B[[i]][,j], Y2_B[[i]][,j], Y3_B[[i]][,j])
    }
  }
  res$t_bar <- rowMeans(res$t, na.rm = TRUE)
  res$t_sd <- apply(res$t, 1, sd, na.rm = TRUE)
  res$pci <- apply(res$t, 1, function(x){
    temp <- sort(x)
    res <- temp[c(floor(length(temp)*(1 - ci.level)/2), floor(length(temp)*(1 - (1 - ci.level)/2)))]
  })
  return(res)
}
