### GEE fo location-scale regression model with parametric form
### ---- 1. mu(x) = x*beta, sig^2(x) = sig^2 ----

gee_1 <- function(par, x, y, intercept = TRUE){
  n_par <- length(par)
  par[n_par] <- exp(par[n_par])
  if(intercept) {
    mu <- as.numeric(cbind(1, x) %*% par[-n_par])
    u1 <- colSums(cbind(1, x) * (y - mu)/par[n_par])
  } else {
    if(inherits(x, "matrix")){
      mu <- as.numeric(x %*% par[-n_par])
      u1 <- colSums(x * (y - mu)/par[n_par])
    } else{
      mu <- x*par[-n_par]
      u1 <- sum(x * (y - mu)/par[n_par])
    }
  }
  u2 <- sum(((y - mu)^2 - par[n_par]))
  return(sum(u1^2) + u2^2)
}

##
x <- runif(100, -1, 1)
z <- rnorm(100, 0, 1)
y <- 0.5*x - 1*z + 0.5*rt(100, df = 3)/sqrt(3)

temp <- lm(y ~ x + z - 1)

c(temp$coeff, var(temp$residuals))

out <- optim(par = c(temp$coeff, 1), fn = gee_1, method = "L-BFGS-B", x = cbind(x, z), y = y,
             intercept = FALSE)
out
c(out$par[1:2], exp(out$par[3]))

out <- nlminb(c(temp$coeff, var(temp$residuals)), gee_1, x = cbind(x, z), y = y, intercept = FALSE)
out
c(out$par[1:2], exp(out$par[3]))

res <- matrix(0, nrow = 3, ncol = 1000)
for(i in 1:1000){
  x <- runif(100, -1, 1)
  y <- 1 + 0.5*x + 1*rnorm(100, 0, 1)
  out <- optim(par = c(mean(y), 1, var(y)), fn = gee_1, method = "L-BFGS-B", x = x, y = y)
  res[,i] <- c(out$par[1:2], exp(out$par[3]))
}
rowMeans(res)
apply(res, 1, sd)
boxplot(t(res))


### ---- 2. mu(x) = x*beta, sig^2(x) = (x*alpha)^2 ----
library(numDeriv)

fun_mu <- function(par, z) as.numeric(z %*% par)
fun_sig2 <- function(par, z) as.numeric((z %*% par)^2)

set.seed(1)
x <- runif(100, 0, 1)
y <- 1 + 0.5*x + (1 + 0.5*x)*rnorm(100, 0, 1)

temp1 <- lm(y ~ x)
temp1$coeff
sd(temp1$resid)

bet_j <- temp1$coeff # c(0, 0.2)
gam_j <- c(sd(temp1$resid), 0) # c(1, 0)
gam_err <- 1
bet_err <- 1
iter <- 0
while((gam_err > 1e-9 | bet_err > 1e-9) & iter < 200){
  mu_j <- fun_mu(bet_j, cbind(1, x))
  sig2_j <- fun_sig2(gam_j, cbind(1, x))
  jac_mu_j <- jacobian(func = fun_mu, x = bet_j, z = cbind(1, x))
  jac_sig2_j <- jacobian(func = fun_sig2, x = gam_j, z = cbind(1, x))
  deta_sig2 <- solve(crossprod(jac_sig2_j/sig2_j)) %*% colSums(jac_sig2_j*((y - mu_j)^2 - sig2_j)/sig2_j^2)
  gam_j <- as.numeric(gam_j + deta_sig2)
  gam_err <- sqrt(sum(deta_sig2^2))
  deta_mu <- solve(crossprod(jac_mu_j/sqrt(sig2_j))) %*% colSums(jac_mu_j*(y - mu_j)/sig2_j)
  bet_j <- as.numeric(bet_j + deta_mu)
  bet_err <- sqrt(sum(deta_mu^2))
  iter <- iter + 1
}
iter
gam_err
bet_err
bet_j
gam_j


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

res <- matrix(0, nrow = 4, ncol = 500)
for(i in 1:500){
  x <- runif(1000, 0, 1)
  y <- 1 + 0.5*x + (1 + 0.5*x)*rnorm(1000, 0, 1)
  temp1 <- lm(y ~ x)
  out <- gee_fiscor(bet_0 = temp1$coeff, gam_0 = c(sd(temp1$resid), 0), fun_mu = fun_mu, fun_sig2 = fun_sig2,
                    z_mu = cbind(1, x), z_sig2 = cbind(1, x), y = y)
  res[,i] <- c(out$beta_mu, out$gam_sig2)
}
rowMeans(res)
boxplot(t(res))



gee_lg <- function(par, x, y, n_par1, n_par2){
  par1 <- par[1:n_par1]
  par2 <- par[(n_par1 + 1):(n_par1 + n_par2)]
  mu <- as.numeric(cbind(1, x) %*% par1)
  sig <- as.numeric(cbind(1, x) %*% par2)
  if(any(sig) <= 0) return(NA)
  else return(sum(((y - mu)^2/sig^2 + 1)*sig)/2)
}

x <- runif(1000, 0, 1)
y <- 1 + 0.5*x + (0.25 + 0.5*x)*rnorm(1000, 0, 1)

temp1 <- lm(y ~ x)
temp1
optim(par = c(temp1$coeff, 0.4, 0.3), fn = gee_lg, method = "BFGS", x = x, y = y, n_par1 = 2, n_par2 = 2)


gee_2 <- function(par, x, y, n_par1, n_par2, intercept = TRUE){
  par1 <- par[1:n_par1]
  par2 <- par[(n_par1 + 1):(n_par1 + n_par2)]
  if(intercept) {
    mu <- as.numeric(cbind(1, x) %*% par1)
    sig <- as.numeric(cbind(1, x) %*% par2)
    u1 <- colSums(cbind(1, x) * (y - mu)/sig^2)
    u2 <- colSums(cbind(1, x) * ((y - mu)^2 - sig^2)/sig^4) #
  } else {
    if(inherits(x, "matrix")){
      mu <- as.numeric(x %*% par1)
      sig <- as.numeric(x %*% par2)
      u1 <- colSums(x * (y - mu)/sig^2)
      u2 <- colSums(x * ((y - mu)^2 - sig^2)/sig^4)
    } else{
      mu <- x*par1
      sig <- x*par2
      u1 <- sum(x * (y - mu)/sig^2)
      u2 <- sum(x*((y - mu)^2 - sig^2)/sig^4)
    }
  }
  return(sum(u1^2) + sum(u2^2))
}

x <- runif(100, 0, 1)
y <- 1 + 0.5*x + (0.5 + 1*x)*rnorm(100, 0, 1)

temp1 <- lm(y ~ x)
optim(par = c(temp1$coeff, 0.2, 0.4), fn = gee_2, method = "L-BFGS-B",
      x = x, y = y, n_par1 = 2, n_par2 = 2, intercept = TRUE)

optim(par = c(temp1$coeff, 0.2, 0.4), fn = gee_2, method = "BFGS",
      x = x, y = y, n_par1 = 2, n_par2 = 2, intercept = TRUE)

gee_21 <- function(par, X, y, n_par1, n_par2, intercept = TRUE){
  par1 <- par[1:n_par1]
  par2 <- par[(n_par1 + 1):(n_par1 + n_par2)]
  if(intercept) {
    mu <- as.numeric(cbind(1, X) %*% par1)
    sig <- as.numeric(cbind(1, X) %*% par2)
    if(any(sig < 0)) u1 <- u2 <- NA
    else{
      u1 <- colSums(cbind(1, X) * (y - mu)/sig^2)
      u2 <- colSums(cbind(1, X) * ((y - mu)^2 - sig^2)/sig^4)
    }
  } else {
    if(inherits(X, "matrix")){
      mu <- as.numeric(X %*% par1)
      sig <- as.numeric(X %*% par2)
      u1 <- colSums(X * (y - mu)/sig)
      u2 <- colSums(X * ((y - mu)^2 - sig)/sig^2)
    } else{
      mu <- X*par1
      sig <- X*par2
      u1 <- sum(X * (y - mu)/sig)
      u2 <- sum(X*((y - mu)^2 - sig)/sig^2)
    }
  }
  return(c(u1, u2))
}

library(BB)
x <- runif(1000, 0, 1)
y <- 1 + 0.5*x + (1 + 0.5*x)*rnorm(1000, 0, 1)

temp1 <- lm(y ~ x)
temp1
mult_par <- rbind(c(temp1$coeff, 0.2, 0.4), c(temp1$coeff, 0.2, 1), c(temp1$coeff, 0, 0.4), c(temp1$coeff, 0, 0.8),
                  c(temp1$coeff, sd(temp1$resid), 0), c(0, 0, 1, 0), c(1, 1, 1, 1))
multiStart(par = mult_par, fn = gee_21, action = "solve",
           X = x, y = y, n_par1 = 2, n_par2 = 2, intercept = TRUE)



library(nleqslv)
nleqslv(x = c(temp1$coeff, 0, 1), fn = gee_21, X = x, y = y, n_par1 = 2, n_par2 = 2, intercept = TRUE,
        method = "Newton", control = list(trace = 1))


gee_3 <- function(par, x, y, n_par1, n_par2){
  par1 <- par[1:n_par1]
  par2 <- par[(n_par1 + 1):(n_par1 + n_par2)]
  mu <- as.numeric(cbind(1, x) %*% par1)
  sig <- as.numeric(cbind(1, x) %*% par2)
  u1 <- colSums(cbind(1, x) * (y - mu)/sig^2)
  u2 <- colSums(2 * cbind(1, x) * (log((y - mu)^2) - log(sig^2))/sig)
  return(sum(u1^2) + sum(u2^2))
}


optim(par = c(temp1$coeff, 0.5, 0.5), fn = gee_3, method = "BFGS",
      x = x, y = y, n_par1 = 2, n_par2 = 2, control = list(maxit = 10000))





