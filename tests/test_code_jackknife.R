fun_mu <- function(par, z) as.numeric(par[1] + par[2]*z)
fun_sig2 <- function(par, z) as.numeric((par[1] + par[2]*z)^2)

prob_dise <- c(1/3, 1/3, 1/3)
vus_np2_est <- vus_np2_est_jack <- numeric(1000)
n <- 150

x_eval <- 0.5
for(i in 1:1000){
  D <- t(rmultinom(n, 1, prob_dise))
  D_f <- apply(D, 1, function(x) which(x == 1))
  X <- runif(n, 0, 1)
  X1 <- X[D_f == 1]
  X2 <- X[D_f == 2]
  X3 <- X[D_f == 3]
  #
  sig_Y1 <- sqrt(fun_sig2(par = c(0.25, 2), z = X1))
  sig_Y2 <- sqrt(fun_sig2(par = c(0.25, 2), z = X2))
  sig_Y3 <- sqrt(fun_sig2(par = c(0.25, 2), z = X3))
  #
  mu_Y1 <- fun_mu(par = c(0.5, 1), z = X1)
  mu_Y2 <- fun_mu(par = c(1, 1.5), z = X2)
  mu_Y3 <- fun_mu(par = c(2.5, 1.5), z = X3)
  #
  Y1 <- mu_Y1 + sig_Y1*rnorm(sum(D_f == 1), 0, 1)
  Y2 <- mu_Y2 + sig_Y2*rnorm(sum(D_f == 2), 0, 1)
  Y3 <- mu_Y3 + sig_Y3*rnorm(sum(D_f == 3), 0, 1)
  #
  out_ll_1 <- reg_locl(X = X1, Y = Y1, x = x_eval)
  out_ll_var_1 <- reg_var_locl(X = X1, resid = out_ll_1$resid, x = x_eval)
  #
  out_ll_2 <- reg_locl(X = X2, Y = Y2, x = x_eval)
  out_ll_var_2 <- reg_var_locl(X = X2, resid = out_ll_2$resid, x = x_eval)
  #
  out_ll_3 <- reg_locl(X = X3, Y = Y3, x = x_eval)
  out_ll_var_3 <- reg_var_locl(X = X3, resid = out_ll_3$resid, x = x_eval)
  #
  out_ll_1_bst <- reg_locl_bts(X = X1, Y = Y1, resid = out_ll_1$resid, x = x_eval, h1 = out_ll_1$h,
                               h2 = out_ll_var_1$h, R = 250)
  out_ll_2_bst <- reg_locl_bts(X = X2, Y = Y2, resid = out_ll_2$resid, x = x_eval, h1 = out_ll_2$h,
                               h2 = out_ll_var_2$h, R = 250)
  out_ll_3_bst <- reg_locl_bts(X = X3, Y = Y3, resid = out_ll_3$resid, x = x_eval, h1 = out_ll_3$h,
                               h2 = out_ll_var_3$h, R = 250)
  #
  vus_np2_est[i] <- vus_np2(mu = c(out_ll_1$mu_x, out_ll_2$mu_x, out_ll_3$mu_x),
                            sig = c(out_ll_var_1$var_x, out_ll_var_2$var_x, out_ll_var_3$var_x),
                            errors = list(out_ll_var_1$error, out_ll_var_2$error, out_ll_var_3$error))
  vus_np2_est_jack[i] <- vus_np2_jack(vus_np2_est[i], mu = c(out_ll_1$mu_x, out_ll_2$mu_x, out_ll_3$mu_x),
                                      sig = c(out_ll_var_1$var_x, out_ll_var_2$var_x, out_ll_var_3$var_x),
                                      errors = list(out_ll_var_1$error, out_ll_var_2$error, out_ll_var_3$error))
}

mean(vus_np2_est)
sd(vus_np2_est)
mean(sqrt(vus_np2_est_jack))


