### ---- Simulating the data ----
### Gamma dítributions
prob_dise <- c(1/3, 1/3, 1/3)
beta_true <- c(-2, 1, -1, 0.5, 0.9, 0.3)
shape_true <- c(1, 1.2, 1.6)

n <- 150

set.seed(285)
D <- t(rmultinom(n, 1, prob_dise))
D_f <- apply(D, 1, function(x) which(x == 1))
X <- runif(n, -2, 2)

X1 <- X[D_f == 1]
X2 <- X[D_f == 2]
X3 <- X[D_f == 3]

Y1 <- rgamma(sum(D_f == 1), shape = shape_true[1], scale = exp(cbind(1, X1) %*% beta_true[1:2])/(shape_true[1]))
Y2 <- rgamma(sum(D_f == 2), shape = shape_true[2], scale = exp(cbind(1, X2) %*% beta_true[3:4])/(shape_true[2]))
Y3 <- rgamma(sum(D_f == 3), shape = shape_true[3], scale = exp(cbind(1, X3) %*% beta_true[5:6])/(shape_true[3]))

x <- 0.5

scale_true <- c(exp(cbind(1, x) %*% beta_true[1:2])/(shape_true[1]),
                exp(cbind(1, x) %*% beta_true[3:4])/(shape_true[2]),
                exp(cbind(1, x) %*% beta_true[5:6])/(shape_true[3]))
vus_gamma(shape = shape_true, scale = scale_true)

### ---- Local-linear regressions ----
out_1 <- reg_locl(X = X1, Y = Y1, x = x)
out_var_1 <- reg_var_locl(X = X1, resid = out_1$resid, x = x)
# adaptive_nn , bwtype = "fixed", bwscaling = TRUE

fit_1 <- glm(Y1 ~ X1, family = Gamma(link = "log"))
res_1 <- summary(fit_1)
mu_1_x <- predict(fit_1, newdata = data.frame(X1 = x), type = "response")

shape_1_est <- 1/res_1$dispersion
scale_1_est <- fit_1$fitted.values/shape_1_est
scale_1_est_x <- mu_1_x/shape_1_est


plot(X1, Y1, pch = 16, ylab = "Y1")
grid()
lines(sort(X1), out_1$mu_X[order(X1)])
lines(sort(X1), fit_1$fitted.values[order(X1)], col = "forestgreen")
lines(sort(X1), exp(cbind(1, sort(X1)) %*% beta_true[1:2]), col = "blue")

plot(X1, abs(out_1$resid), pch = 16, ylab = "Absolute residuals")
grid()
abline(h = 0, lty = 2, col = "gray")
lines(sort(X1), sqrt(out_var_1$var_X[order(X1)]))
lines(sort(X1), sqrt(shape_1_est)*scale_1_est[order(X1)], col = "forestgreen")
lines(sort(X1), exp(cbind(1, sort(X1)) %*% beta_true[1:2])/sqrt(shape_true[1]), col = "blue")


out_2 <- reg_locl(X = X2, Y = Y2, x = x)
out_var_2 <- reg_var_locl(X = X2, resid = out_2$resid, x = x)

fit_2 <- glm(Y2 ~ X2, family = Gamma(link = "log"))
res_2 <- summary(fit_2)
mu_2_x <- predict(fit_2, newdata = data.frame(X2 = x), type = "response")

shape_2_est <- 1/res_2$dispersion
scale_2_est <- fit_2$fitted.values/shape_2_est
scale_2_est_x <- mu_2_x/shape_2_est


plot(X2, Y2, pch = 16, ylab = "Y2")
grid()
lines(sort(X2), out_2$mu_X[order(X2)])
lines(sort(X2), fit_2$fitted.values[order(X2)], col = "forestgreen")
lines(sort(X2), exp(cbind(1, sort(X2)) %*% beta_true[3:4]), col = "blue")

plot(X2, abs(out_2$resid), pch = 16, ylab = "Absolute residuals")
grid()
abline(h = 0, lty = 2, col = "gray")
lines(sort(X2), sqrt(out_var_2$var_X[order(X2)]))
lines(sort(X2), sqrt(shape_2_est)*scale_2_est[order(X2)], col = "forestgreen")
lines(sort(X2), exp(cbind(1, sort(X2)) %*% beta_true[3:4])/sqrt(shape_true[2]), col = "blue")


out_3 <- reg_locl(X = X3, Y = Y3, x = x)
out_var_3 <- reg_var_locl(X = X3, resid = out_3$resid, x = x)

fit_3 <- glm(Y3 ~ X3, family = Gamma(link = "log"))
res_3 <- summary(fit_3)
mu_3_x <- predict(fit_3, newdata = data.frame(X3 = x), type = "response")

shape_3_est <- 1/res_3$dispersion
scale_3_est <- fit_3$fitted.values/shape_3_est
scale_3_est_x <- mu_3_x/shape_3_est


plot(X3, Y3, pch = 16, ylab = "Y3")
grid()
lines(sort(X3), out_3$mu_X[order(X3)])
lines(sort(X3), fit_3$fitted.values[order(X3)], col = "forestgreen")
lines(sort(X3), exp(cbind(1, sort(X3)) %*% beta_true[5:6]), col = "blue")

plot(X3, abs(out_3$resid), pch = 16, ylab = "Absolute residuals")
grid()
abline(h = 0, lty = 2, col = "gray")
lines(sort(X3), sqrt(out_var_3$var_X[order(X3)]))
lines(sort(X3), sqrt(shape_3_est)*scale_3_est[order(X3)], col = "forestgreen")
lines(sort(X3), exp(cbind(1, sort(X3)) %*% beta_true[5:6])/sqrt(shape_true[3]), col = "blue")


vus_1 <- vus_np1(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x), sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                 errors = list(out_var_1$error, out_var_2$error, out_var_3$error))

vus_2 <- vus_np2(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x), sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                 errors = list(out_var_1$error, out_var_2$error, out_var_3$error))

vus_3 <- vus_gamma(shape = c(shape_1_est, shape_2_est, shape_3_est),
                   scale = c(scale_1_est_x, scale_2_est_x, scale_3_est_x))

## constructing the working samples form parametric estimates
mu_gamma <- c(mu_1_x, mu_2_x, mu_3_x)
sig_gamma <- c(shape_1_est*scale_1_est_x^2, shape_2_est*scale_2_est_x^2, shape_3_est*scale_3_est_x^2)
vus_3_np <- vus_np2(mu = mu_gamma, sig = sig_gamma,
                    errors = list(fit_1$resid/sqrt(sig_gamma[1]), fit_2$resid/sqrt(sig_gamma[2]),
                                  fit_3$resid/sqrt(sig_gamma[3])))

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
  error_B[[i]] <- list(out_1_bst$error[,i], out_3_bst$error[,i], out_3_bst$error[,i])
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
fit_1_bst <- boot(data = data.frame(X = X1, Y = Y1), statistic = bs_gamma, R = 250, formula = Y ~ X, x = x)
fit_1_bst

fit_2_bst <- boot(data = data.frame(X = X2, Y = Y2), statistic = bs_gamma, R = 250, formula = Y ~ X, x = x)
fit_2_bst

fit_3_bst <- boot(data = data.frame(X = X3, Y = Y3), statistic = bs_gamma, R = 250, formula = Y ~ X, x = x)
fit_3_bst

vus_3_bst <- vus_gamma_bst(shape_B = rbind(fit_1_bst$t[,1], fit_2_bst$t[,1], fit_3_bst$t[,1]),
                           scale_B = rbind(fit_1_bst$t[,2], fit_2_bst$t[,2], fit_3_bst$t[,2]))
vus_3_bst$t
vus_3_bst$t_bar
vus_3_bst$t_sd
hist(vus_3_bst$t)

