# Load necessary package
library(stats)

# Function to sample from the conditional distribution of a bivariate Gaussian copula
rGaussian_cond <- function(v, rho, n_samples = 1) {
  # Step 1: Compute z1 (the inverse CDF of the standard normal)
  z1 <- qnorm(v)
  
  # Step 2: Sample from a standard normal distribution
  z2_samples <- rnorm(n_samples)
  
  # Step 3: Compute the conditional mean and standard deviation
  conditional_mean <- rho * z1
  conditional_sd <- sqrt(1 - rho^2)
  
  # Step 4: Calculate z_u
  z_u_samples <- conditional_mean + conditional_sd * z2_samples
  
  # Step 5: Transform back to uniform scale
  u_samples <- pnorm(z_u_samples)
  
  return(u_samples)
}

# Example usage:
# v <- 0.7  # given U1 value
# rho <- 1  # correlation coefficient
# n_samples <- 10 # number of samples to generate
# 
# samples <- rGaussian_cond (v, rho, n_samples)
# print(samples)

lam2rho <- function(lam) {
  rho <- (exp(2*lam)-1)/(exp(2*lam)+1)
}

dGaussian <- function(rho, data) {
  x <- qnorm(data[1])
  y <- qnorm(data[2])
  (1-rho^2)^{-1/2}*exp(-0.5*(rho^2*x^2-2*rho*x*y+rho^2*y^2)/(1-rho^2))
}


# # Example usage:
# install.packages("copula")
# library(copula)
# lam <- 1
# # Define the Gaussian copula with a correlation parameter rho
# rho <- lam2rho(lam) 
# rho# example value
# gaussian_cop <- normalCopula(param = rho, dim = 2)
# 
# # Define the point (u, v) where you want to compute the density
# u <- 0.4
# v <- 0.6
# 
# # Compute and print the density of the Gaussian copula at (u, v)
# density <- dCopula(c(u, v), gaussian_cop)
# print(density)
# dGaussian(rho, c(u,v))

lik_at_t <- function(rho_t, data_t) {
  prod(apply(data_t, 1, function(row) dGaussian(rho_t, row)))
}

loglik_at_t <- function(rho_t, data_t) {
  sum(apply(data_t, 1, function(row) log(dGaussian(rho_t, row))))
}


lik_at_t <- function(lam_t, data_t) {
  rho_t <- lam2rho(lam_t)
  prod(apply(data_t, 1, function(row) dGaussian(rho_t, row)))
}

loglik_at_t <- function(lam_t, data_t) {
  rho_t <- lam2rho(lam_t)
  sum(apply(data_t, 1, function(row) log(dGaussian(rho_t, row))))
}
# u <- runif(10)
# v <- runif(10)
# data_t <- cbind(u,v)
# rho_t <- 0.5
# lik_at_t(rho_t, data_t)
# exp(loglik_at_t(rho_t, data_t))
########################
prior_lam_t <- function(lam_t_mius, lam_t, lam_t_plus, alp, beta, tau) {
  if (is.na(lam_t_mius)) {
    dnorm(lam_t_plus, alp+beta*lam_t, sqrt(1/tau))*
      dnorm(lam_t,0,sqrt(1/tau))
  } else if (is.na(lam_t_plus)) {
    dnorm(lam_t, alp+beta*lam_t_mius, sqrt(1/tau))
  } else {
    dnorm(lam_t, alp+beta*lam_t_mius, sqrt(1/tau))*
      dnorm(lam_t_plus, alp+beta*lam_t, sqrt(1/tau))
  }
  
}

log_prior_lam_t <- function(lam_t_mius, lam_t, lam_t_plus, alp, beta, tau) {
  if (is.na(lam_t_mius)) {
    dnorm(lam_t_plus, alp+beta*lam_t, sqrt(1/tau),log=T)+
      dnorm(lam_t,0,sqrt(1/tau), log=T)
  } else if (is.na(lam_t_plus)) {
    dnorm(lam_t, alp+beta*lam_t_mius, sqrt(1/tau), log = T)
  } else {
    dnorm(lam_t, alp+beta*lam_t_mius, sqrt(1/tau), log = T)+
      dnorm(lam_t_plus, alp+beta*lam_t, sqrt(1/tau), log=T)
  }
}

########################
prior_alp <- function(alp, tau_alp) {
  dnorm(alp,0, sqrt(1/tau_alp))
}

log_prior_alp <- function(alp, tau_alp) {
  dnorm(alp,0, sqrt(1/tau_alp),log=T)
}

########################
prior_beta <- function(beta, tau_beta) {
  dnorm(beta,0, sqrt(1/tau_beta))
}

log_prior_beta <- function(beta, tau_beta) {
  dnorm(beta,0, sqrt(1/tau_beta), log=T)
}

########################
prior_tau <- function(tau, shape_tau, rate_tau) {
  dgamma(tau, shape_tau,rate_tau)
}

log_prior_tau <- function(tau, shape_tau, rate_tau) {
  dgamma(tau, shape_tau,rate_tau, log=T)
}

########################
prior_lam_all <- function(lam_all, alp, beta, tau) {
  dense <- dnorm(lam_all[1],0,sqrt(1/tau))
  len <- length(lam_all)
  for (i in 2:len) {
    dense <- dense* dnorm(lam_all[i], alp+beta*lam_all[i-1], sqrt(1/tau))
  }
}


log_prior_lam_all <- function(lam_all, alp, beta, tau) {
  dense <- dnorm(lam_all[1],0,sqrt(1/tau),log=T)
  len <- length(lam_all)
  for (i in 2:len) {
    dense <- dense + dnorm(lam_all[i], alp+beta*lam_all[i-1], sqrt(1/tau), log=T)
  }
}
