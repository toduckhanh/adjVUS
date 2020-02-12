#' @import utils
#' @import numDeriv

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
    deta_sig2 <- solve(crossprod(jac_sig2_j/sig2_j)) %*% colSums(jac_sig2_j*((y - mu_j)^2 - sig2_j)/sig2_j^2)
    gam_j <- as.numeric(gam_j + deta_sig2)
    gam_err <- sqrt(sum(deta_sig2^2))
    deta_mu <- solve(crossprod(jac_mu_j/sqrt(sig2_j))) %*% colSums(jac_mu_j*(y - mu_j)/sig2_j)
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

