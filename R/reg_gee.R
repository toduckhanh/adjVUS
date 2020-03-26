#' @import utils
#' @import stats
#' @import numDeriv
#' @import MASS
#' @import edfun

#### Using own routine
#' @export
gee_fiscor <- function(bet_0, gam_0, fun_mu, fun_sig2, z_mu, z_sig2, y){
  bet_j <- bet_0
  gam_j <- gam_0
  gam_err <- 1
  bet_err <- 1
  iter <- 0
  while((gam_err > 1e-9 | bet_err > 1e-9) & iter < 200){
    mu_j <- fun_mu(bet_j, z = z_mu)
    sig2_j <- fun_sig2(gam_j, z = z_sig2)
    jac_mu_j <- jacobian(func = fun_mu, x = bet_j, z = z_mu)
    jac_sig2_j <- jacobian(func = fun_sig2, x = gam_j, z = z_sig2)
    deta_sig2 <- ginv(crossprod(jac_sig2_j/sig2_j)) %*% colSums(jac_sig2_j*((y - mu_j)^2 - sig2_j)/sig2_j^2)
    gam_j <- as.numeric(gam_j + deta_sig2)
    gam_err <- sqrt(sum(deta_sig2^2))
    deta_mu <- ginv(crossprod(jac_mu_j/sqrt(sig2_j))) %*% colSums(jac_mu_j*(y - mu_j)/sig2_j)
    bet_j <- as.numeric(bet_j + deta_mu)
    bet_err <- sqrt(sum(deta_mu^2))
    iter <- iter + 1
  }
  res <- list()
  res$beta_mu <- bet_j
  res$gam_sig2 <- gam_j
  res$iter <- iter
  res$convergence_mu <- res$convergence_sig2 <- 0
  if(bet_err > 1e-9) res$convergence_mu <- 1
  if(gam_err > 1e-9) res$convergence_sig2 <- 1
  return(res)
}

#' @export
gee_const <- function(beta_0, sig2_0, z, y, fun_mu){
  ## target function
  ff <- function(par, z, y, fun_mu){
    n_par <- length(par)
    par[n_par] <- exp(par[n_par])
    mu <- fun_mu(par = par[-n_par], z = z)
    jac_mu <- jacobian(fun_mu, x = par[-n_par], z = z)
    if(n_par > 2) u1 <- colSums(jac_mu * (y - mu)/par[n_par])
    else u1 <- sum(jac_mu * (y - mu)/par[n_par])
    u2 <- sum(((y - mu)^2 - par[n_par]))
    return(sum(u1^2) + u2^2)
  }
  out <- optim(par = c(beta_0, sig2_0), fn = ff, method = "L-BFGS-B", z = z, y = y, fun_mu = fun_mu,
               control = list(maxit = 500))
  res <- list()
  res$par <- c(out$par[1:length(beta_0)], exp(out$par[length(out$par)]))
  res$convergence <- out$convergence
  return(res)
}

#' @export
reg_gee <- function(X_mu, x_eval, X_sig2 = NULL, Y, fun_mu, fun_sig2 = NULL, bet_start = NULL, gam_start = NULL){
  if(is.null(bet_start) & is.null(gam_start)){
    temp <- lm(Y ~ X_mu - 1)
    bet_0 <- temp$coeff
    sig2_0 <- var(temp$resid)
  }
  else{
    bet_0 <- bet_start
    sig2_0 <- gam_start
  }
  res <- list()
  if(is.null(X_sig2)){ # constant variance
    out_const <- gee_const(beta_0 = bet_0, sig2_0 = sig2_0, z = X_mu, y = Y, fun_mu = fun_mu)
    par_const <- out_const$par
    mu_X <- fun_mu(par_const[-length(par_const)], X_mu)
    mu_x <- fun_mu(par_const[-length(par_const)], x_eval)
    var_x <- par_const[length(par_const)]
    res$mess <- "Constant variance"
    res$mu_x <- mu_x
    res$mu_X <- mu_X
    res$var_x <- res$var_X <- var_x
    res$error <- (Y - res$mu_X)/sqrt(res$var_X)
    res$y_hat <- matrix(0, nrow = length(res$error), ncol = length(x_eval))
    for(i in 1:length(x_eval)) {
      res$y_hat[,i] <- res$mu_x[i] + sqrt(res$var_x)*res$error
    }
    res$bet_mu <- par_const[-length(par_const)]
    res$gam_sig2 <- log(par_const[length(par_const)])
    res$bet_start <- bet_0
    res$gam_start <- sig2_0
  } else{ # non-constant variance
    if(is.null(gam_start)) gam_0 <- c(sqrt(sig2_0), rep(0, ncol(X_sig2) - 1))
    else gam_0 <- gam_start
    out_nconst <- gee_fiscor(bet_0 = bet_0, gam_0 = gam_0, fun_mu = fun_mu, fun_sig2 = fun_sig2,
                             z_mu = X_mu, z_sig2 = X_sig2, y = Y)
    mu_X <- fun_mu(out_nconst$beta_mu, X_mu)
    mu_x <- fun_mu(out_nconst$beta_mu, x_eval)
    var_X <- fun_sig2(out_nconst$gam_sig2, X_sig2)
    var_x <- fun_sig2(out_nconst$gam_sig2, x_eval)
    res$mess <- "non-constant variance"
    res$mu_x <- mu_x
    res$mu_X <- mu_X
    res$var_x <- var_x
    res$var_X <- var_X
    res$error <- (Y - res$mu_X)/sqrt(res$var_X)
    res$y_hat <- matrix(0, nrow = length(res$error), ncol = length(x_eval))
    for(i in 1:length(x_eval)) {
      res$y_hat[,i] <- res$mu_x[i] + sqrt(res$var_x[i])*res$error
    }
    res$bet_mu <- out_nconst$beta_mu
    res$gam_sig2 <- out_nconst$gam_sig2
    res$conveg <- c(out_nconst$convergence_mu, out_nconst$convergence_sig2)
    res$bet_start <- bet_0
    res$gam_start <- gam_0
  }
  return(res)
}

#' @export
reg_gee_bts_1 <- function(X_mu, x_eval, X_sig2 = NULL, Y, fun_mu, fun_sig2 = NULL, bet_start, gam_start, R = 250){
  n <- length(Y)
  ind <- replicate(R, sample(1:n, size = n, replace = TRUE))
  fun_boot <- function(ind, X_mu, X_sig2 = NULL, Y, x_eval){
    if(inherits(X_mu, "matrix")) X_mu.b <- X_mu[ind,]
    else X_mu.b <- X_mu[ind]
    Y.b <- Y[ind]
    if(is.null(X_sig2)) X_sig2.b <- NULL
    else{
      if(inherits(X_sig2, "matrix")) X_sig2.b <- X_sig2[ind,]
      else X_sig2.b <- X_sig2[ind]
    }
    out_gee.b <- reg_gee(X_mu.b, x_eval, X_sig2.b, Y.b, fun_mu, fun_sig2, bet_start = bet_start, gam_start = gam_start)
    res.b <- list()
    res.b$mu_x <- out_gee.b$mu_x
    res.b$var_x <- out_gee.b$var_x
    res.b$error <- out_gee.b$error
    res.b$y_hat <- out_gee.b$y_hat
    return(res.b)
  }
  res_B <- apply(ind, 2, function(i) {
    out <- try(fun_boot(i, X_mu = X_mu, X_sig2 = X_sig2, Y = Y, x_eval = x_eval), silent = TRUE)
    if(!inherits(out, "try-error")) return(out)
    else return(NA)
    })
  if(inherits(x_eval, "matrix")) nx <- nrow(x_eval)
  else nx <- length(x_eval)
  mu_x_B <- var_x_B <- matrix(NA, nrow = nx, ncol = R)
  error_B <- matrix(NA, nrow = n, ncol = R)
  Y_B <- replicate(R, matrix(NA, nrow = n, ncol = nx), simplify = FALSE)
  for(i in 1:R){
    if(all(!is.na(res_B[[i]]))){
      mu_x_B[,i] <- res_B[[i]]$mu_x
      var_x_B[,i] <- res_B[[i]]$var_x
      error_B[,i] <- res_B[[i]]$error
      Y_B[[i]] <- res_B[[i]]$y_hat
    }
  }
  res <- list()
  res$mu_x_B <- mu_x_B
  res$var_x_B <- var_x_B
  res$sd.mu_x <- apply(mu_x_B, 1, sd, na.rm = TRUE)
  res$sd.var_x <- apply(var_x_B, 1, sd, na.rm = TRUE)
  res$error <- error_B
  res$y_hat <- Y_B
  class(res) <- "reg_gee_bts"
  return(res)
}

### Boostrap procedure uses the simulation of errors from the estimated distribution
#' @export
reg_gee_bts_2 <- function(out_reg_gee, X_mu, X_sig2, x_eval, fun_mu, fun_sig2, bet_start, gam_start, R = 250){
  err_gee_funs <- edfun(x = out_reg_gee$error)
  err_gee_rfun <- err_gee_funs$rfun
  n <- length(out_reg_gee$error)
  if(inherits(x_eval, "matrix")) nx <- nrow(x_eval)
  else nx <- length(x_eval)
  mu_x_B <- var_x_B <- matrix(NA, nrow = nx, ncol = R)
  error_B <- matrix(NA, nrow = n, ncol = R)
  for(i in 1:R){
    err_gee_b <- err_gee_rfun(n)
    Y_gee_b <- out_reg_gee$mu_X + sqrt(out_reg_gee$var_X)*err_gee_b
    if(!missing(fun_sig2)){
      out_gee_b <- try(reg_gee(X_mu = X_mu, x_eval = x_eval, X_sig2 = X_sig2, Y = Y_gee_b, fun_mu = fun_mu,
                               fun_sig2 = fun_sig2, bet_start = bet_start, gam_start = gam_start), silent = TRUE)
    } else{
      out_gee_b <- try(reg_gee(X_mu = X_mu, x_eval = x_eval, Y = Y_gee_b, fun_mu = fun_mu, bet_start = bet_start,
                               gam_start = gam_start), silent = TRUE)
    }
    if(!inherits(out_gee_b, "try-error")){
      mu_x_B[,i] <- out_gee_b$mu_x
      var_x_B[,i] <- out_gee_b$var_x
      error_B[,i] <- out_gee_b$error
    }
  }
  res <- list()
  res$mu_x_B <- mu_x_B
  res$var_x_B <- var_x_B
  res$sd.mu_x <- apply(mu_x_B, 1, sd, na.rm = TRUE)
  res$sd.var_x <- apply(var_x_B, 1, sd, na.rm = TRUE)
  res$error <- error_B
  class(res) <- "reg_gee_bts"
  return(res)
}



