###---- Try with 1000 times ----
### Gamma dítributions
prob_dise <- c(1/3, 1/3, 1/3)
beta_true <- c(-2, 1, -1, 0.5, 0.9, 0.3)
shape_true <- c(1, 1.2, 1.6)

n <- 150
x <- 0.5
trials <- 1000
res_test <- matrix(0, nrow = 3, ncol = trials)
for(i in 1:1000){
  set.seed(i)
  D <- t(rmultinom(n, 1, prob_dise))
  D_f <- apply(D, 1, function(x) which(x == 1))
  X <- runif(n, -2, 2)
  X1 <- X[D_f == 1]
  X2 <- X[D_f == 2]
  X3 <- X[D_f == 3]
  Y1 <- rgamma(sum(D_f == 1), shape = shape_true[1], scale = exp(cbind(1, X1) %*% beta_true[1:2])/(shape_true[1]))
  Y2 <- rgamma(sum(D_f == 2), shape = shape_true[2], scale = exp(cbind(1, X2) %*% beta_true[3:4])/(shape_true[2]))
  Y3 <- rgamma(sum(D_f == 3), shape = shape_true[3], scale = exp(cbind(1, X3) %*% beta_true[5:6])/(shape_true[3]))
  ##
  out_1 <- reg_locl(X = X1, Y = Y1, x = x)
  out_var_1 <- reg_var_locl(X = X1, resid = out_1$resid, x = x)
  out_2 <- reg_locl(X = X2, Y = Y2, x = x)
  out_var_2 <- reg_var_locl(X = X2, resid = out_2$resid, x = x)
  out_3 <- reg_locl(X = X3, Y = Y3, x = x)
  out_var_3 <- reg_var_locl(X = X3, resid = out_3$resid, x = x)
  res_test[1,i] <- vus_np1(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x),
                           sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                           errors = list(out_var_1$error, out_var_2$error, out_var_3$error))
  res_test[2,i] <- vus_np2(mu = c(out_1$mu_x, out_2$mu_x, out_3$mu_x),
                           sig = c(out_var_1$var_x, out_var_2$var_x, out_var_3$var_x),
                           errors = list(out_var_1$error, out_var_2$error, out_var_3$error))
  ##
  fit_1 <- glm(Y1 ~ X1, family = Gamma(link = "log"))
  res_1 <- summary(fit_1)
  mu_1_x <- predict(fit_1, newdata = data.frame(X1 = x), type = "response")
  fit_2 <- glm(Y2 ~ X2, family = Gamma(link = "log"))
  res_2 <- summary(fit_2)
  mu_2_x <- predict(fit_2, newdata = data.frame(X2 = x), type = "response")
  fit_3 <- glm(Y3 ~ X3, family = Gamma(link = "log"))
  res_3 <- summary(fit_3)
  mu_3_x <- predict(fit_3, newdata = data.frame(X3 = x), type = "response")
  shape_est <- c(1/res_1$dispersion, 1/res_2$dispersion, 1/res_3$dispersion)
  scale_est <- c(mu_1_x/(shape_est[1]), mu_2_x/(shape_est[2]), mu_3_x/(shape_est[3]))
  res_test[3,i] <- vus_gamma(shape = shape_est, scale = scale_est)
}

rowMeans(res_test)

boxplot(res_test[1, -which.max(res_test[1,])], res_test[2,], res_test[3,], ylim = c(0,1))
abline(h = c(1/6, 1), lty = 2, col = "gray")
