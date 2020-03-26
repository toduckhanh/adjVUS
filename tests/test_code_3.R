### ---- Simulating the data ----
### Normal distributions

prob_dise <- c(1/3, 1/3, 1/3)
beta_true <- c(1, 0.5, 2, 0.8, 3.5, 1.2)
sigma_e <- sqrt(c(0.3, 0.8, 1.3))

x <- 1

mu_true <- c(beta_true[1:2] %*% c(1, x), beta_true[3:4] %*% c(1, x), beta_true[5:6] %*% c(1, x))
vus_normal(mu_true, sigma_e)

n <- 150

set.seed(308)
D <- t(rmultinom(n, 1, prob_dise))
D_f <- apply(D, 1, function(x) which(x == 1))
X <- runif(n, -2, 2)
X1 <- X[D_f == 1]
X2 <- X[D_f == 2]
X3 <- X[D_f == 3]

Y1 <- as.numeric(cbind(1, X1) %*% beta_true[1:2] + rnorm(sum(D_f == 1), mean = 0, sd = sigma_e[1]))
Y2 <- as.numeric(cbind(1, X2) %*% beta_true[3:4] + rnorm(sum(D_f == 2), mean = 0, sd = sigma_e[2]))
Y3 <- as.numeric(cbind(1, X3) %*% beta_true[5:6] + rnorm(sum(D_f == 3), mean = 0, sd = sigma_e[3]))

fun_mu <- function(par, z) as.numeric(z %*% par)

out_gee_1 <- reg_gee(X_mu = cbind(1, X1), x_eval = cbind(1, x), Y = Y1, fun_mu = fun_mu)
out_gee_1

out_gee_1_bst <- reg_gee_bts_1(X_mu = cbind(1, X1), x_eval = cbind(1, x), Y = Y1, fun_mu = fun_mu,
                               bet_start = out_gee_1$bet_start, gam_start = out_gee_1$gam_start, R = 250)
out_gee_1_bst

out_gee_2 <- reg_gee(X_mu = cbind(1, X2), x_eval = cbind(1, x), Y = Y2, fun_mu = fun_mu)
out_gee_2

out_gee_2_bst <- reg_gee_bts_1(X_mu = cbind(1, X2), x_eval = cbind(1, x), Y = Y2, fun_mu = fun_mu,
                               bet_start = out_gee_2$bet_start, gam_start = out_gee_2$gam_start, R = 250)
out_gee_2_bst

out_gee_3 <- reg_gee(X_mu = cbind(1, X3), x_eval = cbind(1, x), Y = Y3, fun_mu = fun_mu)
out_gee_3

out_gee_3_bst <- reg_gee_bts_1(X_mu = cbind(1, X3), x_eval = cbind(1, x), Y = Y3, fun_mu = fun_mu,
                               bet_start = out_gee_3$bet_start, gam_start = out_gee_3$gam_start, R = 250)
out_gee_3_bst

vus_np1(mu = c(out_gee_1$mu_x, out_gee_2$mu_x, out_gee_3$mu_x),
        sig = c(out_gee_1$var_x, out_gee_2$var_x, out_gee_3$var_x),
        errors = list(out_gee_1$error, out_gee_2$error, out_gee_3$error))

vus_np1_est <- vus_np1_fy(Y1 = out_gee_1$y_hat, Y2 = out_gee_2$y_hat, Y3 = out_gee_3$y_hat)

vus_np2(Y1 = out_gee_1$y_hat, Y2 = out_gee_2$y_hat, Y3 = out_gee_3$y_hat)

vus_np1_fy_bst(Y1_B = out_gee_1_bst$y_hat, Y2_B = out_gee_2_bst$y_hat, Y3_B = out_gee_3_bst$y_hat,
               bw1 = vus_np1_est$bws[1], bw2 = vus_np1_est$bws[2], bw3 = vus_np1_est$bws[3])

vus_np2_bst(Y1_B = out_gee_1_bst$y_hat, Y2_B = out_gee_2_bst$y_hat, Y3_B = out_gee_3_bst$y_hat)

mu_gee_x_B <- rbind(out_gee_1_bst$mu_x_B, out_gee_2_bst$mu_x_B, out_gee_3_bst$mu_x_B)
var_gee_x_B <- rbind(out_gee_1_bst$var_x_B, out_gee_2_bst$var_x_B, out_gee_3_bst$var_x_B)
error_gee_B <- list()
for(i in 1:250){
  error_gee_B[[i]] <- list(out_gee_1_bst$error[,i], out_gee_2_bst$error[,i], out_gee_3_bst$error[,i])
}

vus_1_bst <- vus_np1_bst(mu_gee_x_B, var_gee_x_B, error_gee_B)
vus_1_bst$t_bar
vus_1_bst$t_sd
hist(vus_1_bst$t)
sort(vus_1_bst$t)[c(floor(250*0.05/2), floor(250*(1 - 0.05/2)))]

vus_2_bst <- vus_np2_bst(mu_gee_x_B, var_gee_x_B, error_gee_B)
vus_2_bst$t_bar
vus_2_bst$t_sd
hist(vus_2_bst$t)
sort(vus_2_bst$t)[c(floor(250*0.05/2), floor(250*(1 - 0.05/2)))]




fit_1 <- glm(Y1 ~ X1, family = gaussian(link = "identity"))
res_1 <- summary(fit_1)
fit_2 <- glm(Y2 ~ X2, family = gaussian(link = "identity"))
res_2 <- summary(fit_2)
fit_3 <- glm(Y3 ~ X3, family = gaussian(link = "identity"))
res_3 <- summary(fit_3)
mu_est <- c(fit_1$coeff %*% c(1, x), fit_2$coeff %*% c(1, x), fit_3$coeff %*% c(1, x))
sigma_est <- sqrt(c(res_1$dispersion, res_2$dispersion, res_3$dispersion))
vus_3 <- vus_normal(mu_est, sigma_est)


## constructing the working samples form parametric estimates
vus_3_np <- vus_np2(mu = mu_est, sig = sigma_est^2,
                    errors = list(fit_1$resid/sigma_est[1], fit_2$resid/sigma_est[2], fit_3$resid/sigma_est[3]))

## nonparametric estimation
out_1 <- reg_locl(X = X1, Y = Y1, x = x)
out_var_1 <- reg_var_locl(X = X1, resid = out_1$resid, x = x)

out_2 <- reg_locl(X = X2, Y = Y2, x = x)
out_var_2 <- reg_var_locl(X = X2, resid = out_2$resid, x = x)

out_3 <- reg_locl(X = X3, Y = Y3, x = x)
out_var_3 <- reg_var_locl(X = X3, resid = out_3$resid, x = x)

vus_1 <- vus_np1(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x), sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                 errors = list(out_var_1$error, out_var_2$error, out_var_3$error))

vus_2 <- vus_np2(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x), sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                 errors = list(out_var_1$error, out_var_2$error, out_var_3$error))
c(vus_1, vus_2, vus_3)

## bootstrap
out_1_bst <- reg_locl_bts(X = X1, Y = Y1, resid = out_1$resid, x = x, h1 = out_1$h, h2 = out_var_1$h, R = 250)
c(mean(out_1_bst$mu_x_B), out_1$mu_x)
c(mean(out_1_bst$var_x_B), out_var_1$var_x)

out_2_bst <- reg_locl_bts(X = X2, Y = Y2, resid = out_2$resid, x = x, h1 = out_2$h, h2 = out_var_2$h, R = 250)
c(mean(out_2_bst$mu_x_B), out_2$mu_x)
c(mean(out_2_bst$var_x_B), out_var_2$var_x)

out_3_bst <- reg_locl_bts(X = X3, Y = Y3, resid = out_3$resid, x = x, h1 = out_3$h, h2 = out_var_3$h, R = 250)
c(mean(out_3_bst$mu_x_B), out_3$mu_x)
c(mean(out_3_bst$var_x_B), out_var_3$var_x)

mu_x_B <- rbind(out_1_bst$mu_x_B, out_2_bst$mu_x_B, out_3_bst$mu_x_B)
var_x_B <- rbind(out_1_bst$var_x_B, out_2_bst$var_x_B, out_3_bst$var_x_B)
error_B <- list()
for(i in 1:250){
  error_B[[i]] <- list(out_1_bst$error[,i], out_2_bst$error[,i], out_3_bst$error[,i])
}

vus_1_bst <- vus_np1_bst(mu_x_B, var_x_B, error_B)
vus_1_bst$t_bar
vus_1_bst$t_sd
hist(vus_1_bst$t)
sort(vus_1_bst$t)[c(floor(250*0.05/2), floor(250*(1 - 0.05/2)))]

vus_2_bst <- vus_np2_bst(mu_x_B, var_x_B, error_B)
vus_2_bst$t_bar
vus_2_bst$t_sd
hist(vus_2_bst$t)
sort(vus_2_bst$t)[c(floor(250*0.05/2), floor(250*(1 - 0.05/2)))]


require(boot)
fit_1_bst <- boot(data = data.frame(X = X1, Y = Y1), statistic = bs_normal, R = 250, formula = Y ~ X, x = x)
fit_1_bst

fit_2_bst <- boot(data = data.frame(X = X2, Y = Y2), statistic = bs_normal, R = 250, formula = Y ~ X, x = x)
fit_2_bst

fit_3_bst <- boot(data = data.frame(X = X3, Y = Y3), statistic = bs_normal, R = 250, formula = Y ~ X, x = x)
fit_3_bst

vus_3_bst <- vus_normal_bst(mu_B = rbind(fit_1_bst$t[,1], fit_2_bst$t[,1], fit_3_bst$t[,1]),
                            sigma_B = rbind(fit_1_bst$t[,2], fit_2_bst$t[,2], fit_3_bst$t[,2]))
vus_3_bst$t
vus_3_bst$t_bar
vus_3_bst$t_sd
hist(vus_3_bst$t)

