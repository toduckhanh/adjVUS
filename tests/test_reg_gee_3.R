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
    res$error <- (Y - mu_X)/sqrt(res$var_X)
    res$beta <- par_const[-length(par_const)]
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
    res$error <- (Y - mu_X)/sqrt(res$var_X)
    res$beta_mu <- out_nconst$beta_mu
    res$gam_sig2 <- out_nconst$gam_sig2
    res$conveg <- c(out_nconst$convergence_mu, out_nconst$convergence_sig2)
    res$bet_start <- bet_0
    res$gam_start <- gam_0
  }
  return(res)
}

reg_gee_bts <- function(X_mu, x_eval, X_sig2 = NULL, Y, fun_mu, fun_sig2 = NULL, bet_start, gam_start, R = 250){
  n <- length(Y)
  ind <- replicate(R, sample(1:n, size = n, replace = TRUE))
  fun_boot <- function(ind, X_mu, X_sig2 = NULL, Y, x_eval){
    X_mu.b <- X_mu[ind,]
    Y.b <- Y[ind]
    if(is.null(X_sig2)) X_sig2.b <- NULL
    else X_sig2.b <- X_sig2[ind,]
    out_gee.b <- reg_gee(X_mu.b, x_eval, X_sig2.b, Y.b, fun_mu, fun_sig2, bet_start = bet_start, gam_start = gam_start)
    res.b <- list()
    res.b$mu_x <- out_gee.b$mu_x
    res.b$var_x <- out_gee.b$var_x
    res.b$error <- out_gee.b$error
    return(res.b)
  }
  res_B <- apply(ind, 2, function(i) fun_boot(i, X_mu = X_mu, X_sig2 = X_sig2, Y = Y, x_eval = x_eval))
  mu_x_B <- var_x_B <- matrix(0, nrow = nrow(x_eval), ncol = R)
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
  class(res) <- "reg_gee_bts"
  return(res)
}


