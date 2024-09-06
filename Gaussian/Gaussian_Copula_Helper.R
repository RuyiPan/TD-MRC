# Load necessary package
library(stats)

# Function to sample from the conditional distribution of a bivariate Gaussian copula
rGaussian_cond <- function(n_samples = 1, v, rho) {
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
  (exp(2*lam)-1)/(exp(2*lam)+1)
}

dGaussian <- function(rho, data, log=T) {
  x <- qnorm(data[1])
  y <- qnorm(data[2])
  if (any(data %in% c(0, 1))) {
    if (log==T) {
      log(0)
    } else {
      0
    }
  } else {
    if (log) {
      -1/2*log(1-rho^2)+ (-0.5*(rho^2*x^2-2*rho*x*y+rho^2*y^2)/(1-rho^2))
    } else {
      (1-rho^2)^{-1/2}*exp(-0.5*(rho^2*x^2-2*rho*x*y+rho^2*y^2)/(1-rho^2))
    }
  }
 
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
# u <- 1
# v <- 0
# 
# # Compute and print the density of the Gaussian copula at (u, v)
# density <- dCopula(c(u, v), gaussian_cop)
# print(density)
# dGaussian(rho, c(u,v))

# lik_at_t <- function(rho_t, data_t) {
#   prod(apply(data_t, 1, function(row) dGaussian(rho_t, row)))
# }
# 
# loglik_at_t <- function(rho_t, data_t) {
#   sum(apply(data_t, 1, function(row) log(dGaussian(rho_t, row))))
# }


lik_at_t <- function(lam_t, data_t) {
  rho_t <- lam2rho(lam_t)
  prod(apply(data_t, 1, function(row) dGaussian(rho_t, row)))
}

loglik_at_t <- function(lam_t, data_t) {
  rho_t <- lam2rho(lam_t)
  sum(apply(data_t, 1, function(row) dGaussian(rho_t, row, log=T)))
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
prior_lam_all <- function(lam_all, alp, beta, tau) {
  dense <- dnorm(lam_all[1],0,sqrt(1/tau))
  len <- length(lam_all)
  for (i in 2:len) {
    dense <- dense* dnorm(lam_all[i], alp+beta*lam_all[i-1], sqrt(1/tau))
  }
  dense
}


log_prior_lam_all <- function(lam_all, alp, beta, tau) {
  dense <- dnorm(lam_all[1],0,sqrt(1/tau),log=T)
  len <- length(lam_all)
  for (i in 2:len) {
    dense <- dense + dnorm(lam_all[i], alp+beta*lam_all[i-1], sqrt(1/tau), log=T)
  }
  dense
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



#LPML
library(foreach)
library(doParallel)
CPO.component <- function(U, rho_t) {
  1/ dGaussian(rho_t, U, log=F)
}

CPO <- function(U, rho_all) {
  L <- length(rho_all)
  1/mean(do.call(rbind,lapply(1:L, 
                              function(l) CPO.component(U,rho_all[l]))))
}


LPML <- function(U_all, rho_all) {
  TT <- length(U_all)
  GlobalFunctions = ls(globalenv())
  ncores <- detectCores()
  cl <- parallel::makeCluster(ncores-2)
  doParallel::registerDoParallel(cl)
  
  res_CPO <- foreach(t=1:TT, .combine="c",
                     .export = GlobalFunctions)%dopar%{
                       nt <- nrow(U_all[[t]])
                       CPO_t <- c(1:nt)
                       for (i in 1:nt) {
                         CPO_t[i] <- CPO(U_all[[t]][i,], rho_all[,t])
                       }
                       list(CPO_t)
                     }
  parallel::stopCluster(cl)
  
  log_CPOs <- log(unlist(res_CPO)) 
  LPML <- sum(log_CPOs[log_CPOs !=-Inf])
  ALPML <- mean(log_CPOs[log_CPOs !=-Inf])
  return(LPML)
}




rho_all <- apply(lam_all[9000:9200,], c(1,2), lam2rho)
U_all <- data$test

CPO(U_all[[32]][4,],rho_all[, 32])

dGaussian(rho_all[1,32], U_all[[32]][4,], log=F)

LPML(U_all, rho_all)




ll.mg.U <- function(U, rho_t) {
  dGaussian(rho_t, U, log=T)
}


Devi <- function(U_all, rho_all) {
  TT <- length(U_all)
  D.single <- 0
  for (t in 1:TT) {
    nt <- nrow(U_all[[t]])
    for (i in 1:nt) {
      temp_t <- ll.mg.U(U_all[[t]][i,], rho_all[t])
      if (temp_t != -Inf) {
        D.single <-  D.single+ll.mg.U(U_all[[t]][i,], rho_all[t])
      }
    }
  }
  
  return (-2*D.single)
}



DIC <- function(U_all, lam_all) {
  rho_all <- apply(lam_all, c(1,2), lam2rho)
  rho_par <- lam2rho(colMeans(lam_all))
  L <- dim(rho_all)[1]
  GlobalFunctions = ls(globalenv())
  ncores <- detectCores()
  cl <- parallel::makeCluster(ncores-2)
  doParallel::registerDoParallel(cl)
  D_all <- foreach(l=1:L, .combine = "rbind",.export = GlobalFunctions)%dopar% {
    Devi(U_all, rho_all[l, ])
  }
  parallel::stopCluster(cl)
  PMD <- mean(D_all)
  
  D_at_mean <- Devi(U_all, rho_par)
  
  DIC <- 2*PMD-D_at_mean
  return (DIC)
}


DIC(U_all, lam_all[9000:10000,])

