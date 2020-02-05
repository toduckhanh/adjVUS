#' @import utils

### ---- VUS with three Normal distributions ----
#' @export
vus_normal <- function(mu, sigma){
  a <- sigma[2]/sigma[1]
  b <- (mu[1] - mu[2])/sigma[1]
  c <- sigma[2]/sigma[3]
  d <- (mu[3] - mu[2])/sigma[3]
  integrate(function(x, a, b, c, d) pnorm(a*x - b)*pnorm(-c*x + d)*dnorm(x), a = a, b = b, c = c, d = d,
            lower = -Inf, upper = Inf)$value
}

#' @export
bs_normal <- function(formula, data, ind, x) {
  d <- data[ind,]
  fit <- glm(formula, data = d, family = gaussian(link = "identity"))
  res <- summary(fit)
  mu_x <- predict(fit, newdata = data.frame(X = x), type = "response")
  sig_est <- sqrt(res$dispersion)
  return(c(mu_x, sig_est))
}

#' @export
vus_normal_bst <- function(mu_B, sigma_B){
  res <- list()
  res$t <- numeric(ncol(mu_B))
  for(i in 1:ncol(mu_B)){
    res$t[i] <- vus_normal(mu = c(mu_B[1,i], mu_B[2,i], mu_B[3,i]),
                          sigma = c(sigma_B[1,i], sigma_B[2,i], sigma_B[3,i]))
  }
  res$t_bar <- mean(res$t)
  res$t_sd <- sd(res$t)
  return(res)
}

### ---- VUS with three Gamma distributions ----
#' @export
vus_gamma <- function(shape, scale){
  integrate(function(x) pgamma(x, shape = shape[1], scale = scale[1])*
              pgamma(x, shape = shape[3], scale = scale[3], lower.tail = FALSE)*
              dgamma(x, shape = shape[2], scale = scale[2]),
            lower = -Inf, upper = Inf)$value
}

#' @export
bs_gamma <- function(formula, data, ind, x) {
  d <- data[ind,]
  fit <- glm(formula, data = d, family = Gamma(link = "log"))
  res <- summary(fit)
  mu_x <- predict(fit, newdata = data.frame(X = x), type = "response")
  shape_est <- 1/res$dispersion
  scale_est_x <- mu_x/shape_est
  return(c(shape_est, scale_est_x))
}

#' @export
vus_gamma_bst <- function(shape_B, scale_B){
  res <- list()
  res$t <- numeric(ncol(shape_B))
  for(i in 1:ncol(shape_B)){
    res$t[i] <- vus_gamma(shape = c(shape_B[1,i], shape_B[2,i], shape_B[3,i]),
                          scale = c(scale_B[1,i], scale_B[2,i], scale_B[3,i]))
  }
  res$t_bar <- mean(res$t)
  res$t_sd <- sd(res$t)
  return(res)
}


