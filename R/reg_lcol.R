#' @import utils
#' @import np
#' @import edfun

#' @export
reg_locl <- function(X, Y, x_eval, ckertype = "epanechnikov", bwmethod = "cv.ls", bwscaling = FALSE, bws = NA, ...){
  if(is.na(bws)){
    out_locl <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod, ckertype = ckertype, bwscaling = bwscaling, ...)
  } else{
    out_locl <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod, ckertype = ckertype, bwscaling = bwscaling,
                      bws = bws, ...)
  }
  res <- list()
  res$resid <- Y - out_locl$mean
  res$mu_x <- predict(out_locl, newdata = data.frame(X = x_eval))
  res$h <- out_locl$bw
  res$mu_X <- out_locl$mean
  class(res) <- "reg_locl"
  return(res)
}

#' @export
reg_var_locl <- function(X, out_reg_locl, x_eval, ckertype = "epanechnikov", bwmethod = "cv.ls", bwscaling = FALSE,
                         bws = NA, ...){
  Y <- log(out_reg_locl$resid^2 + 1/length(out_reg_locl$resid))
  if(is.na(bws)){
    out_locl_var <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod,
                          ckertype = ckertype, bwscaling = bwscaling, ...)
  } else{
    out_locl_var <- npreg(Y ~ X, regtype = "ll", bwmethod = bwmethod,
                          ckertype = ckertype, bwscaling = bwscaling, bws = bws, ...)
  }
  d_hat <- 1/mean(out_reg_locl$resid^2*exp(-out_locl_var$mean))
  res <- list()
  res$h <- out_locl_var$bw
  res$var_x <- exp(predict(out_locl_var, newdata = data.frame(X = x_eval)))/d_hat
  res$var_X <- exp(out_locl_var$mean)/d_hat
  res$error <- out_reg_locl$resid/sqrt(res$var_X)
  res$y_hat <- matrix(0, nrow = length(res$error), ncol = length(x_eval))
  for(i in 1:length(x_eval)) {
    res$y_hat[,i] <- out_reg_locl$mu_x[i] + sqrt(res$var_x[i])*res$error
  }
  class(res) <- "reg_var_locl"
  return(res)
}

#' @export
reg_locl_bts_1 <- function(X, Y, x_eval, h1, h2, ckertype = c("epanechnikov", "epanechnikov"),
                           bwmethod = c("cv.ls", "cv.ls"), bwscaling = c(FALSE, FALSE), R = 250, ...){
  n <- length(Y)
  ind <- replicate(R, sample(1:n, size = n, replace = TRUE))
  fun_boot <- function(ind, X, Y, x_eval){
    X.b <- X[ind]
    Y.b <- Y[ind]
    out_mu.b <- reg_locl(X.b, Y.b, x_eval, ckertype = ckertype[1], bwmethod = bwmethod[1], bwscaling = bwscaling[1],
                         bws = h1) #, ...
    out_var.b <- reg_var_locl(X.b, out_mu.b, x_eval, ckertype = ckertype[2], bwmethod = bwmethod[2],
                              bwscaling = bwscaling[2], bws = h2) # , ...
    res.b <- list()
    res.b$mu_x <- out_mu.b$mu_x
    res.b$var_x <- out_var.b$var_x
    res.b$error <- out_var.b$error
    res.b$y_hat <- out_var.b$y_hat
    return(res.b)
  }
  res_B <- apply(ind, 2, function(i) fun_boot(i, X = X, Y = Y, x_eval = x_eval))
  mu_x_B <- var_x_B <- matrix(0, nrow = length(x_eval), ncol = R)
  error_B <- matrix(0, nrow = n, ncol = R)
  Y_B <- replicate(R, matrix(NA, nrow = n, ncol = length(x_eval)), simplify = FALSE)
  for(i in 1:R){
    mu_x_B[,i] <- res_B[[i]]$mu_x
    var_x_B[,i] <- res_B[[i]]$var_x
    error_B[,i] <- res_B[[i]]$error
    Y_B[[i]] <- res_B[[i]]$y_hat
  }
  res <- list()
  res$mu_x_B <- mu_x_B
  res$var_x_B <- var_x_B
  res$sd.mu_x <- apply(mu_x_B, 1, sd)
  res$sd.var_x <- apply(var_x_B, 1, sd)
  res$error <- error_B
  res$y_hat <- Y_B
  class(res) <- "reg_locl_bts"
  return(res)
}


### Boostrap procedure uses the simulation of errors from the estimated distributions
#' @export
reg_locl_bts_2 <- function(out_reg_locl, out_reg_var_locl, X, x, h1, h2, ckertype = c("epanechnikov", "epanechnikov"),
                           bwmethod = c("cv.ls", "cv.ls"), bwscaling = c(FALSE, FALSE), R = 250, ...){
  n <- length(out_reg_var_locl$error)
  mu_x_B <- var_x_B <- matrix(0, nrow = length(x), ncol = R)
  error_B <- matrix(0, nrow = n, ncol = R)
  err_ll_funs <- edfun(x = out_reg_var_locl$error)
  err_ll_rfun <- err_ll_funs$rfun
  for(i in 1:R){
    err_ll_b <- err_ll_rfun(n)
    Y1_b <- out_reg_locl$mu_X + sqrt(out_reg_var_locl$var_X)*err_ll_b
    out_ll_b <- reg_locl(X = X, Y = Y1_b, x = x, bws = out_reg_locl$h)
    out_ll_var_b <- reg_var_locl(X = X, resid = out_ll_b$resid, x = x, bws = out_reg_var_locl$h)
    mu_x_B[,i] <- out_ll_b$mu_x
    var_x_B[,i] <- out_ll_var_b$var_x
    error_B[,i] <- out_ll_var_b$error
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
