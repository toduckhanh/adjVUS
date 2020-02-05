#' @import utils
#' @import np

#' @export
reg_locl <- function(X, Y, x, ckertype = "epanechnikov", bwmethod = "cv.ls", bwscaling = FALSE, bws = NA, ...){
  if(is.na(bws)){
    out_locl <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod, ckertype = ckertype, bwscaling = bwscaling, ...)
  } else{
    out_locl <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod, ckertype = ckertype, bwscaling = bwscaling,
                      bws = bws, ...)
  }
  resid <- Y - out_locl$mean
  mu_x <- predict(out_locl, newdata = data.frame(X = x))
  res <- list(h = out_locl$bw, mu_x = mu_x, resid = resid, mu_X = out_locl$mean)
  class(res) <- "reg_locl"
  return(res)
}

#' @export
reg_var_locl <- function(X, resid, x, ckertype = "epanechnikov", bwmethod = "cv.ls", bwscaling = FALSE, bws = NA, ...){
  Y <- log(resid^2 + 1/length(resid))
  if(is.na(bws)){
    out_locl_var <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod,
                          ckertype = ckertype, bwscaling = bwscaling, ...)
  } else{
    out_locl_var <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod,
                          ckertype = ckertype, bwscaling = bwscaling, bws = bws, ...)
  }
  d_hat <- 1/mean(resid^2*exp(-out_locl_var$mean))
  var_x <- exp(predict(out_locl_var, newdata = data.frame(X = x)))/d_hat
  var_X <- exp(out_locl_var$mean)/d_hat
  error <- resid/sqrt(var_X)
  res <- list(h = out_locl_var$bw, var_x = var_x, error = error, var_X = var_X)
  class(res) <- "reg_var_locl"
  return(res)
}

#' @export
reg_locl_bts <- function(X, Y, resid, x, h1, h2, ckertype = c("epanechnikov", "epanechnikov"),
                         bwmethod = c("cv.ls", "cv.ls"), bwscaling = c(FALSE, FALSE), R = 250, ...){
  n <- length(Y)
  ind <- replicate(R, sample(1:n, size = n, replace = TRUE))
  fun_boot <- function(ind, X, Y, resid, x){
    X.b <- X[ind]
    Y.b <- Y[ind]
    resid.b <- resid[ind]
    out_mu.b <- reg_locl(X.b, Y.b, x, ckertype = ckertype[1], bwmethod = bwmethod[1], bwscaling = bwscaling[1],
                         bws = h1) #, ...
    out_var.b <- reg_var_locl(X.b, resid.b, x, ckertype = ckertype[2], bwmethod = bwmethod[2],
                              bwscaling = bwscaling[2], bws = h2) # , ...
    res.b <- list()
    res.b$mu_x <- out_mu.b$mu_x
    res.b$var_x <- out_var.b$var_x
    res.b$error <- out_var.b$error
    return(res.b)
  }
  res_B <- apply(ind, 2, function(i) fun_boot(i, X = X, Y = Y, resid = resid, x = x))
  mu_x_B <- var_x_B <- matrix(0, nrow = length(x), ncol = R)
  error_B <- matrix(0, nrow = n, ncol = R)
  for(i in 1:R){
    mu_x_B[,i] <- res_B[[i]]$mu_x
    var_x_B[,i] <- res_B[[i]]$var_x
    error_B[,i] <- res_B[[i]]$error
  }
  res <- list()
  res$mu_x_B <- mu_x_B
  res$var_x_B <- var_x_B
  res$sd.mu_x <- apply(mu_x_B, 1, sd)
  res$sd.var_x <- apply(var_x_B, 1, sd)
  res$error <- error_B
  class(res) <- "reg_locl_bts"
  return(res)
}

