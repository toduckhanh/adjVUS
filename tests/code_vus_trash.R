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




# #' @export
vus_np2_jack <- function(vus_est, mu, sig, errors){
  Y1 <- mu[1] + sqrt(sig[1])*errors[[1]]
  Y2 <- mu[2] + sqrt(sig[2])*errors[[2]]
  Y3 <- mu[3] + sqrt(sig[3])*errors[[3]]
  n1 <- length(Y1)
  n2 <- length(Y2)
  n3 <- length(Y3)
  vus_est_mult <- vus_est * n1 * n2 * n3
  var_term <- vusC_full_core(Y1, Y2, Y3)
  var_term_i <- sapply(1:n1, function(x) sum(var_term$ind1[-x]))
  var_term_j <- sapply(1:n2, function(x) sum(var_term$ind2[-x]))
  var_term_k <- sapply(1:n3, function(x) sum(var_term$ind3[-x]))
  ans_var <- var((vus_est_mult - var_term_i)/(n2*n3))/n1 +
    var((vus_est_mult - var_term_j)/(n1*n3))/n2 +
    var((vus_est_mult - var_term_k)/(n1*n2))/n3
  return(ans_var)
}

// [[Rcpp::export]]
List vusC_full_core(NumericVector tt1, NumericVector tt2, NumericVector tt3){
  int nn1 = tt1.size(), nn2 = tt2.size(), nn3 = tt3.size();
  NumericVector out1(nn1), out2(nn2), out3(nn3);
  NumericMatrix I_ij(nn2, nn1);
  NumericMatrix I_ik(nn3, nn1);
  for(int i = 0; i < nn1; i++){
    NumericMatrix M(nn3, nn2);
    for(int j = 0; j < nn2; j++){
      for(int k = 0; k < nn3; k++){
        M(k,j) = indvus(tt1[i], tt2[j], tt3[k]);
      }
    }
    NumericVector temp = colSums(M);
    I_ij(_, i) = temp;
    I_ik(_, i) = rowSums(M);
    out1[i] = sum(temp);
  }
  out2 = rowSums(I_ij);
  out3 = rowSums(I_ik);
  List out;
  out["ind1"] = out1;
  out["ind2"] = out2;
  out["ind3"] = out3;
  return out;
}

# #' @import BB
# #' @export
gee_BB <- function(par, fun_mu, fun_sig2, z_mu, z_sig2, y, np1, ...){
  ff_BB <- function(par, fun_mu, fun_sig2, z_mu, z_sig2, y, np1){
    bet_j <- par[1:np1]
    gam_j <- par[(np1 + 1):length(par)]
    mu_j <- fun_mu(bet_j, z = z_mu)
    res <- numeric(length(par))
    jac_mu_j <- jacobian(func = fun_mu, x = bet_j, z = z_mu)
    if(missing(fun_sig2)){
      sig2_j <- gam_j <- exp(gam_j)
      res[(np1 + 1):length(par)] <- sum(((y - mu_j)^2 - sig2_j))
    } else{
      sig2_j <- fun_sig2(gam_j, z = z_sig2)
      jac_sig2_j <- jacobian(func = fun_sig2, x = gam_j, z = z_sig2)
      res[(np1 + 1):length(par)] <- colSums(jac_sig2_j*((y - mu_j)^2 - sig2_j)/sig2_j^2)
    }
    res[1:np1] <- colSums(jac_mu_j*(y - mu_j)/sig2_j)
    return(res)
  }
  res <- list()
  if(missing(fun_sig2)){
    out_BB <- BBsolve(par, fn = ff_BB, fun_mu = fun_mu, z_mu = X1, y = Y1, np1 = 2, ...)
    res$beta_mu <- out_BB$par[1:np1]
    res$gam_sig2 <- exp(out_BB$par[(np1 + 1):length(par)])
  } else {
    out_BB <- BBsolve(par, fn = ff_BB, fun_mu = fun_mu, fun_sig2 = fun_sig2, z_mu = X1,
                      z_sig2 = X1, y = Y1, np1 = 2, ...)
    res$beta_mu <- out_BB$par[1:np1]
    res$gam_sig2 <- out_BB$par[(np1 + 1):length(par)]
  }
  res$convergence <- out_BB$convergence
  res$error <- out_BB$residual
  return(res)
}

#' @export
reg_gee_BB <- function(X_mu, x_eval, X_sig2, Y, fun_mu, fun_sig2, par_start, np1, ...){
  res <- list()
  if(missing(fun_sig2)){ # constant variance
    out_const <- gee_BB(par = par_start, z_mu = X_mu, fun_mu = fun_mu, y = Y, np1 = np1, ...)
    par_const <- out_const$par
    mu_X <- fun_mu(out_const$beta_mu, X_mu)
    mu_x <- fun_mu(out_const$beta_mu, x_eval)
    var_x <- out_const$gam_sig2
    res$mess <- "Constant variance"
    res$mu_x <- mu_x
    res$mu_X <- mu_X
    res$var_x <- res$var_X <- var_x
    res$error <- (Y - mu_X)/sqrt(res$var_X)
    res$beta <- par_const[-length(par_const)]
    res$par_start <- par_start
    res$convergence <- out_const$convergence
  } else{ # non-constant variance
    out_nconst <- gee_BB(par = par_start, z_mu = X_mu, z_sig2 = X_sig2, fun_mu = fun_mu, fun_sig2 = fun_sig2,
                         y = Y, np1 = np1, ...)
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
    res$convergence <- out_nconst$convergence
    res$par_start <- par_start
  }
  return(res)
}
