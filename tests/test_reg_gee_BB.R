require(BB)
require(numDeriv)

gee_BB <- function(par, fun_mu, fun_sig2, z_mu, z_sig2, y, np1){
  bet_j <- par[1:np1]
  gam_j <- par[(np1 + 1):length(par)]
  mu_j <- fun_mu(bet_j, z = z_mu)
  sig2_j <- fun_sig2(gam_j, z = z_sig2)
  jac_mu_j <- jacobian(func = fun_mu, x = bet_j, z = z_mu)
  jac_sig2_j <- jacobian(func = fun_sig2, x = gam_j, z = z_sig2)
  res <- numeric(length(par))
  res[1:np1] <- colSums(jac_mu_j*(y - mu_j)/sig2_j)
  res[(np1 + 1):length(par)] <- colSums(jac_sig2_j*((y - mu_j)^2 - sig2_j)/sig2_j^2)
  return(res)
}

gee_BB(par = c(0.5, 1, 0.25, 2), fun_mu = fun_mu, fun_sig2 = fun_sig2, z_mu = X1, z_sig2 = X1, y = Y1, np1 = 2)

out_gee_1_BB <- BBsolve(c(0.5, 1, 0.25, 2), fn = gee_BB, fun_mu = fun_mu, fun_sig2 = fun_sig2, z_mu = X1,
                        z_sig2 = X1, y = Y1, np1 = 2)

mu_gee_1 <- fun_mu(par = out_gee_1_BB$par[1:2], z = X1)
var_gee_1 <- fun_sig2(par = out_gee_1_BB$par[3:4], z = X1)
error_gee_1 <- (Y1 - mu_gee_1)/sqrt(var_gee_1)

err_gee_1_funs <- edfun(x = error_gee_1)
err_gee_1_rfun <- err_gee_1_funs$rfun

err_gee_1_b <- err_gee_1_rfun(sum(D_f == 1))
Y1_gee_b <- mu_gee_1 + sqrt(var_gee_1)*err_gee_1_b

BBsolve(c(0.5, 1, 0.25, 2), fn = gee_BB, fun_mu = fun_mu, fun_sig2 = fun_sig2, z_mu = X1,
        z_sig2 = X1, y = Y1_gee_b, np1 = 2)

gee_fiscor(bet_0 = c(0.5, 1), gam_0 = c(0.25, 2), fun_mu = fun_mu, fun_sig2 = fun_sig2, z_mu = X1, z_sig2 = X1,
           y = Y1_gee_b)



