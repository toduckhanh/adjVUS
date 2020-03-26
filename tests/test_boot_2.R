### test bootstrap procedure by simulating the error terms

fun_mu <- function(par, z) as.numeric(par[1] + par[2]*z)
fun_sig2 <- function(par, z) as.numeric((par[1] + par[2]*z)^2)

prob_dise <- c(1/3, 1/3, 1/3)
n <- 150
D <- t(rmultinom(n, 1, prob_dise))
D_f <- apply(D, 1, function(x) which(x == 1))
X <- runif(n, 0, 1)
X1 <- X[D_f == 1]
X2 <- X[D_f == 2]
X3 <- X[D_f == 3]

sig_Y1 <- sqrt(fun_sig2(par = c(0.25, 2), z = X1))
sig_Y2 <- sqrt(fun_sig2(par = c(0.25, 2), z = X2))
sig_Y3 <- sqrt(fun_sig2(par = c(0.25, 2), z = X3))

mu_Y1 <- fun_mu(par = c(0.5, 1), z = X1)
mu_Y2 <- fun_mu(par = c(1, 1.5), z = X2)
mu_Y3 <- fun_mu(par = c(2.5, 1.5), z = X3)

Y1 <- mu_Y1 + sig_Y1*rnorm(sum(D_f == 1), 0, 1)
Y2 <- mu_Y2 + sig_Y2*rnorm(sum(D_f == 2), 0, 1)
Y3 <- mu_Y3 + sig_Y3*rnorm(sum(D_f == 3), 0, 1)

x_eval <- 0.5
out_ll_1 <- reg_locl(X = X1, Y = Y1, x = x_eval)
out_ll_var_1 <- reg_var_locl(X = X1, resid = out_ll_1$resid, x = x_eval)

out_ll_2 <- reg_locl(X = X2, Y = Y2, x = x_eval)
out_ll_var_2 <- reg_var_locl(X = X2, resid = out_ll_2$resid, x = x_eval)

out_ll_3 <- reg_locl(X = X3, Y = Y3, x = x_eval)
out_ll_var_3 <- reg_var_locl(X = X3, resid = out_ll_3$resid, x = x_eval)

vus_np2_ll <- vus_np2(mu = c(out_ll_1$mu_x, out_ll_2$mu_x, out_ll_3$mu_x),
                      sig = c(out_ll_var_1$var_x, out_ll_var_2$var_x, out_ll_var_3$var_x),
                      errors = list(out_ll_var_1$error, out_ll_var_2$error, out_ll_var_3$error))

require(edfun)
system.time({
err_ll_1_funs <- edfun(x = out_ll_var_1$error)
err_ll_1_rfun <- err_ll_1_funs$rfun

err_ll_2_funs <- edfun(x = out_ll_var_2$error)
err_ll_2_rfun <- err_ll_2_funs$rfun

err_ll_3_funs <- edfun(x = out_ll_var_3$error)
err_ll_3_rfun <- err_ll_3_funs$rfun

R <- 250
mu_ll_1_x_b <- var_ll_1_x_b <- mu_ll_2_x_b <- var_ll_2_x_b <- mu_ll_3_x_b <- var_ll_3_x_b <- numeric(R)
error_ll_1_b <- matrix(0, nrow = sum(D_f == 1), ncol = R)
error_ll_2_b <- matrix(0, nrow = sum(D_f == 2), ncol = R)
error_ll_3_b <- matrix(0, nrow = sum(D_f == 3), ncol = R)

for(i in 1:R){
  err_ll_1_b <- err_ll_1_rfun(sum(D_f == 1))
  Y1_b <- out_ll_1$mu_X + sqrt(out_ll_var_1$var_X)*err_ll_1_b
  #
  out_ll_1_b <- reg_locl(X = X1, Y = Y1_b, x = x_eval, bws = out_ll_1$h)
  out_ll_var_1_b <- reg_var_locl(X = X1, resid = out_ll_1_b$resid, x = x_eval, bws = out_ll_var_1$h)
  #
  mu_ll_1_x_b[i] <- out_ll_1_b$mu_x
  var_ll_1_x_b[i] <- out_ll_var_1_b$var_x
  error_ll_1_b[,i] <- out_ll_var_1_b$error
  #
  err_ll_2_b <- err_ll_2_rfun(sum(D_f == 2))
  Y2_b <- out_ll_2$mu_X + sqrt(out_ll_var_2$var_X)*err_ll_2_b
  #
  out_ll_2_b <- reg_locl(X = X2, Y = Y2_b, x = x_eval, bws = out_ll_2$h)
  out_ll_var_2_b <- reg_var_locl(X = X2, resid = out_ll_2_b$resid, x = x_eval, bws = out_ll_var_2$h)
  #
  mu_ll_2_x_b[i] <- out_ll_2_b$mu_x
  var_ll_2_x_b[i] <- out_ll_var_2_b$var_x
  error_ll_2_b[,i] <- out_ll_var_2_b$error
  #
  err_ll_3_b <- err_ll_3_rfun(sum(D_f == 3))
  Y3_b <- out_ll_3$mu_X + sqrt(out_ll_var_3$var_X)*err_ll_3_b
  #
  out_ll_3_b <- reg_locl(X = X3, Y = Y3_b, x = x_eval, bws = out_ll_3$h)
  out_ll_var_3_b <- reg_var_locl(X = X3, resid = out_ll_3_b$resid, x = x_eval, bws = out_ll_var_3$h)
  #
  mu_ll_3_x_b[i] <- out_ll_3_b$mu_x
  var_ll_3_x_b[i] <- out_ll_var_3_b$var_x
  error_ll_3_b[,i] <- out_ll_var_3_b$error
}
})
mean(mu_ll_1_x_b)
mean(var_ll_1_x_b)

mu_ll_x_B_1 <- rbind(mu_ll_1_x_b, mu_ll_2_x_b, mu_ll_3_x_b)
var_ll_x_B_1 <- rbind(var_ll_1_x_b, var_ll_2_x_b, var_ll_3_x_b)
error_ll_B_1 <- list()
for(i in 1:250){
  error_ll_B_1[[i]] <- list(error_ll_1_b[,i], error_ll_2_b[,i], error_ll_3_b[,i])
}

vus_np2_B_1 <- vus_np2_bst(mu_ll_x_B_1, var_ll_x_B_1, error_ll_B_1)

system.time({
out_ll_1_bst <- reg_locl_bts(X = X1, Y = Y1, resid = out_ll_1$resid, x = x_eval, h1 = out_ll_1$h,
                             h2 = out_ll_var_1$h, R = 250)
out_ll_2_bst <- reg_locl_bts(X = X2, Y = Y2, resid = out_ll_2$resid, x = x_eval, h1 = out_ll_2$h,
                             h2 = out_ll_var_2$h, R = 250)
out_ll_3_bst <- reg_locl_bts(X = X3, Y = Y3, resid = out_ll_3$resid, x = x_eval, h1 = out_ll_3$h,
                             h2 = out_ll_var_3$h, R = 250)
})
mu_ll_x_B_2 <- rbind(out_ll_1_bst$mu_x_B, out_ll_2_bst$mu_x_B, out_ll_3_bst$mu_x_B)
var_ll_x_B_2 <- rbind(out_ll_1_bst$var_x_B, out_ll_2_bst$var_x_B, out_ll_3_bst$var_x_B)
error_ll_B_2 <- list()
for(i in 1:250){
  error_ll_B_2[[i]] <- list(out_ll_1_bst$error[,i], out_ll_2_bst$error[,i], out_ll_3_bst$error[,i])
}

vus_np2_B_2 <- vus_np2_bst(mu_ll_x_B_2, var_ll_x_B_2, error_ll_B_2)





