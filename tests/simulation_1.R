
fun_mu_1 <- function(par, z) as.numeric(par[1] + sin(par[2]*z))
fun_mu_2 <- function(par, z) as.numeric(par[1] + par[2]*sin(par[3]*z))
fun_mu_3 <- function(par, z) as.numeric(par[1] + par[2]*sin(par[3]*z))

fun_sig2 <- function(par, z) as.numeric((z %*% par)^2)

prob_dise <- c(1/3, 1/3, 1/3)

n <- 150
D <- t(rmultinom(n, 1, prob_dise))
D_f <- apply(D, 1, function(x) which(x == 1))
X <- runif(n, 0, 1)
X1 <- X[D_f == 1]
X2 <- X[D_f == 2]
X3 <- X[D_f == 3]

Y1 <- fun_mu_1(par = c(0.2, 2), z = X1) + rnorm(sum(D_f == 1), 0, 1)*fun_sig2(par = c(0.5, 1), z = cbind(1, X1))
Y2 <- fun_mu_2(par = c(0.6, 1.2, 2), z = X2) + rnorm(sum(D_f == 2), 0, 1)*fun_sig2(par = c(0.5, 1), z = cbind(1, X2))
Y3 <- fun_mu_3(par = c(1, 1.5, 2), z = X3) + rnorm(sum(D_f == 3), 0, 1)*fun_sig2(par = c(0.5, 1), z = cbind(1, X3))

out_gee_1 <- reg_gee(X_mu = X1, x_eval = 0.8, Y = Y1, fun_mu = fun_mu_1, bet_start = c(0.2, 2), gam_start = 1.3)
out_gee_2 <- reg_gee(X_mu = X2, x_eval = 0.8, Y = Y2, fun_mu = fun_mu_2, bet_start = c(0.4, 1.2, 2), gam_start = 1.3)
out_gee_3 <- reg_gee(X_mu = X3, x_eval = 0.8, Y = Y3, fun_mu = fun_mu_3, bet_start = c(0.5, 1.5, 2), gam_start = 1.3)

vus_np1(mu = c(out_gee_1$mu_x, out_gee_2$mu_x, out_gee_3$mu_x),
        sig = c(out_gee_1$var_x, out_gee_2$var_x, out_gee_3$var_x),
        errors = list(out_gee_1$error, out_gee_2$error, out_gee_3$error))

vus_np2(mu = c(out_gee_1$mu_x, out_gee_2$mu_x, out_gee_3$mu_x),
        sig = c(out_gee_1$var_x, out_gee_2$var_x, out_gee_3$var_x),
        errors = list(out_gee_1$error, out_gee_2$error, out_gee_3$error))

