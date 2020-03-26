library(adjVUS)

fun_mu_1 <- function(par, z) as.numeric(par[1] + sin(par[2]*z*pi))
fun_mu_2 <- function(par, z) as.numeric(par[1] + par[2]*sin(par[3]*z*pi))
fun_mu_3 <- function(par, z) as.numeric(par[1] + par[2]*sin(par[3]*z))
fun_sig2 <- function(par, z) as.numeric((par[1] + par[2]*z)^2)

x_eval <- seq(0.5, 1.5, length.out = 21)

n1 <- n2 <- n3 <- 50
bet1 <- c(-0.3, 2)
bet2 <- c(1.5, 1, 2)
bet3 <- c(2, 1, 1.5)
gam <- c(0.5, 1.2)

R <- 10

set.seed(1)

X1 <- runif(n1, 0.5, 1.5)
X2 <- runif(n2, 0.5, 1.5)
X3 <- runif(n3, 0.5, 1.5)
sig_Y1 <- sqrt(fun_sig2(par = gam, z = X1))
sig_Y2 <- sqrt(fun_sig2(par = gam, z = X2))
sig_Y3 <- sqrt(fun_sig2(par = gam, z = X3))
mu_Y1 <- fun_mu_1(par = bet1, z = X1)
mu_Y2 <- fun_mu_2(par = bet2, z = X2)
mu_Y3 <- fun_mu_3(par = bet3, z = X3)
Y1 <- mu_Y1 + sig_Y1*rnorm(n1, 0, 1)
Y2 <- mu_Y2 + sig_Y2*rnorm(n2, 0, 1)
Y3 <- mu_Y3 + sig_Y3*rnorm(n3, 0, 1)
#
fun_mu_gee <- function(par, z) as.numeric(par[1] + par[2]*z)

out_gee_1 <- reg_gee(X_mu = X1, x_eval = x_eval, Y = Y1, fun_mu = fun_mu_gee,
                     bet_start = c(0, 1), gam_start = 1)
out_gee_2 <- reg_gee(X_mu = X2, x_eval = x_eval, Y = Y2, fun_mu = fun_mu_gee,
                     bet_start = c(0, 1), gam_start = 1)
out_gee_3 <- reg_gee(X_mu = X3, x_eval = x_eval, Y = Y3, fun_mu = fun_mu_gee,
                     bet_start = c(0, 1), gam_start = 1)

out_gee_1_bst <- reg_gee_bts_1(X_mu = X1, x_eval = x_eval, Y = Y1, fun_mu = fun_mu_gee,
                               bet_start = out_gee_1$bet_mu, gam_start = out_gee_1$gam_sig2, R = R)
out_gee_2_bst <- reg_gee_bts_1(X_mu = X2, x_eval = x_eval, Y = Y2, fun_mu = fun_mu_gee,
                               bet_start = out_gee_2$bet_mu, gam_start = out_gee_2$gam_sig2, R = R)
out_gee_3_bst <- reg_gee_bts_1(X_mu = X3, x_eval = x_eval, Y = Y3, fun_mu = fun_mu_gee,
                               bet_start = out_gee_3$bet_mu, gam_start = out_gee_3$gam_sig2, R = R)
res_gee_bst <- vus_np2_bst(out_gee_1_bst$y_hat, out_gee_2_bst$y_hat, out_gee_3_bst$y_hat)


res_gee_bst
