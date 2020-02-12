### ---- 1. mu(x) = beta0 + x*beta1, sig^2(x) = sig^2 ----
fun_mu <- function(par, z) as.numeric(z %*% par)

x <- runif(100, 0, 1)
y <- fun_mu(par = c(1, 0.5), z = cbind(1, x)) + 1*rnorm(100, 0, 1)
temp1 <- lm(y ~ x)
system.time({
  out1 <- gee_const(beta_0 = temp1$coeff, sig2_0 = 0.5, z = cbind(1, x), y = y, fun_mu = fun_mu)
})
out1

## try to fit with sig^2(x) = (gamma0 + x*gamma1)^2
fun_sig2 <- function(par, z) as.numeric((z %*% par)^2)

system.time({
  out2 <- gee_fiscor(bet_0 = temp1$coeff, gam_0 = c(sd(temp1$resid), 0), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                     z_mu = cbind(1, x), z_sig2 = cbind(1, x), y = y)
})
out2

res_const <- matrix(0, nrow = 3, ncol = 1000)
res_lm <- matrix(0, nrow = 4, ncol = 1000)
for(i in 1:1000){
  x <- runif(500, 0, 1)
  y <- fun_mu(par = c(1, 0.5), z = cbind(1, x)) + 1*rnorm(500, 0, 1)
  temp1 <- lm(y ~ x)
  out1 <- gee_const(beta_0 = temp1$coeff, sig2_0 = 0.5, z = cbind(1, x), y = y, fun_mu = fun_mu)
  out2 <- gee_fiscor(bet_0 = temp1$coeff, gam_0 = c(sd(temp1$resid), 0), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                     z_mu = cbind(1, x), z_sig2 = cbind(1, x), y = y)
  res_const[,i] <- out1$par
  res_lm[,i] <- c(out2$beta_mu, out2$gam_sig2)
}

rowMeans(res_const)
rowMeans(res_lm)

apply(res_const, 1, sd)
apply(res_lm, 1, sd)

boxplot(cbind(t(res_const), t(res_lm)))
hist(res_lm[4,])
mean(res_lm[4,]) + c(-1, 1)*qnorm((0.05/2), lower.tail = FALSE)*sd(res_lm[4,])

### ---- 2. mu(x) = beta0 + x*beta1, sig^2(x) = (gamma0 + x*gamma1)^2 ----
## comparing the performance of parameter estimator and of nonparametric one.

fun_mu <- function(par, z) as.numeric(z %*% par)
fun_sig2 <- function(par, z) as.numeric((z %*% par)^2)

x <- runif(100, 0, 1)
y <- fun_mu(par = c(1, 0.5), z = cbind(1, x)) + sqrt(fun_sig2(par = c(0.5, 0.5), z = cbind(1, x)))*rnorm(100, 0, 1)
temp1 <- lm(y ~ x)

system.time({
  out3 <- gee_fiscor(bet_0 = temp1$coeff, gam_0 = c(sd(temp1$resid), 0), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                     z_mu = cbind(1, x), z_sig2 = cbind(1, x), y = y)
})
out3


x_eval <- seq(0, 1, length.out = 101)
out_locl <- reg_locl(X = x, Y = y, x = x_eval)
out_var_locl <- reg_var_locl(X = x, resid = out_locl$resid, x = x_eval)

mu_true <- fun_mu(c(1, 0.5), cbind(1, x_eval))
mu_par <- fun_mu(out3$beta_mu, cbind(1, x_eval))
mu_nonpar <- out_locl$mu_x

mu_par - mu_true
mu_nonpar - mu_true

(mu_par - mu_true)^2
(mu_nonpar - mu_true)^2


var_true <- fun_sig2(c(0.5, 0.5), cbind(1, x_eval))
var_par <- fun_sig2(out3$gam_sig2, cbind(1, x_eval))
var_nonpar <- out_var_locl$var_x

mean(abs(sqrt(var_par) - sqrt(var_true)))
mean(abs(sqrt(var_nonpar) - sqrt(var_true)))

mean((sqrt(var_par) - sqrt(var_true))^2)
mean((sqrt(var_nonpar) - sqrt(var_true))^2)

###

x_eval <- seq(0, 1, length.out = 51)
par_est <- matrix(0, nrow = 4, ncol = 1000)
made_est <- msde_est <- matrix(0, nrow = 2, ncol = 1000)
bias_par <- bias_nonpar <- matrix(0, nrow = 51, ncol = 1000)
mse_par <- mse_nonpar <- matrix(0, nrow = 51, ncol = 1000)
n <- 200
for(i in 1:1000){
  x <- runif(n, 0, 1)
  y <- fun_mu(par = c(1, 0.5), z = cbind(1, x)) + sqrt(fun_sig2(par = c(0.5, 0.5), z = cbind(1, x)))*rnorm(n, 0, 1)
  temp1 <- lm(y ~ x)
  out3 <- gee_fiscor(bet_0 = temp1$coeff, gam_0 = c(sd(temp1$resid), 0), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                     z_mu = cbind(1, x), z_sig2 = cbind(1, x), y = y)
  out_locl <- reg_locl(X = x, Y = y, x = x_eval)
  out_var_locl <- reg_var_locl(X = x, resid = out_locl$resid, x = x_eval)
  #
  mu_true <- fun_mu(c(1, 0.5), cbind(1, x_eval))
  mu_par <- fun_mu(out3$beta_mu, cbind(1, x_eval))
  mu_nonpar <- out_locl$mu_x
  bias_par[,i] <- mu_par - mu_true
  bias_nonpar[,i] <- mu_nonpar - mu_true
  mse_par[,i] <- (mu_par - mu_true)^2
  mse_nonpar[,i] <- (mu_nonpar - mu_true)^2
  #
  var_true <- fun_sig2(c(0.5, 0.5), cbind(1, x_eval))
  var_par <- fun_sig2(out3$gam_sig2, cbind(1, x_eval))
  var_nonpar <- out_var_locl$var_x
  made_est[,i] <- c(mean(abs(sqrt(var_par) - sqrt(var_true))), mean(abs(sqrt(var_nonpar) - sqrt(var_true))))
  msde_est[,i] <- c(mean((sqrt(var_par) - sqrt(var_true))^2), mean((sqrt(var_nonpar) - sqrt(var_true))^2))
  par_est[,i] <- c(out3$beta_mu, out3$gam_sig2)
}

rowMeans(par_est)
boxplot(t(par_est))

rowMeans(bias_par)
rowMeans(bias_nonpar)

rowMeans(mse_par)
rowMeans(mse_nonpar)

rowMeans(made_est)
rowMeans(msde_est)

boxplot(t(made_est), names = c("Parametric", "Nonparametric"), outline = FALSE)
boxplot(t(msde_est), names = c("Parametric", "Nonparametric"), outline = FALSE)

plot(x_eval, rowMeans(bias_par), type = "l", ylim = c(-0.02, 0.02), ylab = "Bias", xlab = "X")
lines(x_eval, rowMeans(bias_nonpar), lty = 2)
grid()
legend("topleft", legend = c("Parametric", "Nonparametric"), lty = c(1, 2))

plot(x_eval, rowMeans(mse_par), type = "l", ylim = c(0, 0.22), ylab = "MSE", xlab = "X")
lines(x_eval, rowMeans(mse_nonpar), lty = 2)
grid()
legend("topleft", legend = c("Parametric", "Nonparametric"), lty = c(1, 2))

### ---- 3. mu(x) = x + 2*exp(-16*x^2), sig^2(x) = (0.4*exp(-2x^2) + 0.2)^2 ----

fun_mu <- function(par, z) as.numeric(z*par[1] + par[2]*exp(par[3]*z^2))
fun_sig2 <- function(par, z) as.numeric((par[1]*exp(par[2]*z^2) + par[3])^2)

x <- runif(200, -2, 2)
y <- fun_mu(par = c(1, 2, -16), z = x) + sqrt(fun_sig2(par = c(0.4, -2, 0.2), z = x))*rnorm(100, 0, 1)
temp1 <- lm(y ~ x)

system.time({
  out3 <- gee_fiscor(bet_0 = c(1, 2, -16), gam_0 = c(0.4, -2, 0.2), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                     z_mu = x, z_sig2 = x, y = y)
})
out3

x_eval <- seq(-2, 2, length.out = 101)
out_locl <- reg_locl(X = x, Y = y, x = x_eval)
out_var_locl <- reg_var_locl(X = x, resid = out_locl$resid, x = x_eval)

mu_true <- fun_mu(c(1, 2, -16), z = x_eval)
mu_par <- fun_mu(out3$beta_mu, z = x_eval)
mu_nonpar <- out_locl$mu_x

mu_par - mu_true
mu_nonpar - mu_true

(mu_par - mu_true)^2
(mu_nonpar - mu_true)^2

var_true <- fun_sig2(c(0.4, -2, 0.2), z = x_eval)
var_par <- fun_sig2(out3$gam_sig2, z = x_eval)
var_nonpar <- out_var_locl$var_x

mean(abs(sqrt(var_par) - sqrt(var_true)))
mean(abs(sqrt(var_nonpar) - sqrt(var_true)))

mean((sqrt(var_par) - sqrt(var_true))^2)
mean((sqrt(var_nonpar) - sqrt(var_true))^2)

##
x_eval <- seq(-2, 2, length.out = 101)
par_est <- matrix(0, nrow = 6, ncol = 1000)
made_est <- msde_est <- matrix(0, nrow = 2, ncol = 1000)
bias_par <- bias_nonpar <- matrix(0, nrow = 101, ncol = 1000)
mse_par <- mse_nonpar <- matrix(0, nrow = 101, ncol = 1000)

n <- 200
for(i in 1:1000){
  x <- runif(n, -2, 2)
  y <- fun_mu(par = c(1, 2, -16), z = x) + sqrt(fun_sig2(par = c(0.4, -2, 0.2), z = x))*rt(n, df = 3)/sqrt(3)
  out3 <- try(gee_fiscor(bet_0 = c(1, 2, -16), gam_0 = c(0.4, -2, 0.2), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                         z_mu = x, z_sig2 = x, y = y), silent = TRUE)
  out_locl <- reg_locl(X = x, Y = y, x = x_eval)
  out_var_locl <- reg_var_locl(X = x, resid = out_locl$resid, x = x_eval)
  mu_true <- fun_mu(c(1, 2, -16), z = x_eval)
  var_true <- fun_sig2(c(0.4, -2, 0.2), z = x_eval)
  if(class(out3) == "try-error"){
    mu_par <- var_par <- rep(NA, length(x_eval))
    par_est[,i] <- rep(NA, 6)
  } else{
    mu_par <- fun_mu(out3$beta_mu, z = x_eval)
    var_par <- fun_sig2(out3$gam_sig2, z = x_eval)
    par_est[,i] <- c(out3$beta_mu, out3$gam_sig2)
  }
  mu_nonpar <- out_locl$mu_x
  bias_par[,i] <- mu_par - mu_true
  bias_nonpar[,i] <- mu_nonpar - mu_true
  mse_par[,i] <- (mu_par - mu_true)^2
  mse_nonpar[,i] <- (mu_nonpar - mu_true)^2
  #
  var_nonpar <- out_var_locl$var_x
  made_est[,i] <- c(mean(abs(sqrt(var_par) - sqrt(var_true))), mean(abs(sqrt(var_nonpar) - sqrt(var_true))))
  msde_est[,i] <- c(mean((sqrt(var_par) - sqrt(var_true))^2), mean((sqrt(var_nonpar) - sqrt(var_true))^2))
}


rowMeans(par_est, na.rm = TRUE)
boxplot(t(par_est))

rowMeans(bias_par, na.rm = TRUE)
rowMeans(bias_nonpar)

rowMeans(mse_par, na.rm = TRUE)
rowMeans(mse_nonpar)

rowMeans(made_est, na.rm = TRUE)
rowMeans(msde_est, na.rm = TRUE)

boxplot(t(made_est), names = c("Parametric", "Nonparametric"), outline = FALSE)
boxplot(t(msde_est), names = c("Parametric", "Nonparametric"), outline = FALSE)

plot(x_eval, rowMeans(bias_par, na.rm = TRUE), type = "l", ylim = c(-0.21, 0.1), ylab = "Bias", xlab = "X")
lines(x_eval, rowMeans(bias_nonpar), lty = 2)
grid()
legend("topleft", legend = c("Parametric", "Nonparametric"), lty = c(1, 2))

plot(x_eval, rowMeans(mse_par, na.rm = TRUE), type = "l", ylim = c(0, 0.1), ylab = "MSE", xlab = "X")
lines(x_eval, rowMeans(mse_nonpar), lty = 2)
grid()
legend("topleft", legend = c("Parametric", "Nonparametric"), lty = c(1, 2))

