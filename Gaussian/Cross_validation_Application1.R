
source("Gaussian_Copula_Helper.R")

#Read data
set.seed(20240902)
job_num <- 0
# for local 
data <- readRDS(paste0("../Data/","split", job_num,".rds"))
# for server
# data <- readRDS(paste0("/home/panruyi/scratch/MixCopula/Data/","split", job_num,".rds"))
U <- data$train
TT<- length(U)
nts <- c(1:TT)
for (day in c(1:TT)) {
  index <- apply(U[[day]],1, function(row) any(row %in% c(0, 1)))
  U[[day]] <- U[[day]][!index,]
  nts[[day]] <- nrow(U[[day]])
}





#hyper-prameters
tau_alp <- 0.01
tau_beta <- 0.01
shape_tau <- 0.01
rate_tau <- 0.01

#tuning parameters, acceptance rate
#15, 3, 80 , 0.1915617, 0.2563487, 0.05718856
#15, 4, 400, 8
delta_alp <- 12 #15 soso
acc_alp <- 0
delta_beta <- 4 #4
acc_beta <- 0
delta_tau <- 600  #80
acc_tau <- 0
delta_lam <- 8

#initialization
alp <- 0
beta <- 1
tau <- 10
lam <- rep(0.5, TT)
alp_all <- NULL
beta_all <- NULL
tau_all <- NULL
lam_all <- NULL

#mcmc settings
burn_in=100
B=200
batch.size=50
M = B*batch.size

for (b in 1:B) {
  print(b)
  #sampling alpha
  for (it in 1:batch.size) {
      prop_alp <- rnorm(1, alp, sqrt(1/delta_alp))
      log_alp_rate <- log_prior_lam_all(lam, prop_alp, beta, tau)+
        log_prior_alp(prop_alp, tau_alp)+
        dnorm(alp, prop_alp, sqrt(1/delta_alp),log = T)-
        log_prior_lam_all(lam, alp, beta, tau)-
        log_prior_alp(alp, tau_alp)-
        dnorm(prop_alp, alp, sqrt(1/delta_alp),log = T)
      if (log(runif(1)) <= log_alp_rate) {
        alp <- prop_alp
        acc_alp <- c(acc_alp, 1)
      } else {
        acc_alp <- c(acc_alp, 0)
      }
      alp_all <- c(alp_all, alp)
      
      #sampling beta
      prop_beta <- rnorm(1, beta, sqrt(1/delta_beta))
      log_beta_rate <- log_prior_lam_all(lam, alp, prop_beta, tau)+
        log_prior_beta(prop_beta, tau_beta)+
        dnorm(beta, prop_beta, sqrt(1/delta_beta),log = T)-
        log_prior_lam_all(lam, alp, beta, tau)-
        log_prior_alp(beta, tau_beta)-
        dnorm(prop_beta,beta, sqrt(1/delta_beta),log = T)
      if (log(runif(1)) <= log_beta_rate) {
        beta <- prop_beta
        acc_beta <- c(acc_beta, 1)
      } else {
        acc_beta <- c(acc_beta, 0)
      }
      beta_all <- c(beta_all, beta)
      
      #sampling tau
      prop_tau <- rgamma(1, delta_tau, delta_tau/tau)
      log_tau_rate <- log_prior_lam_all(lam, alp, beta, prop_tau)+
        log_prior_tau(prop_tau, shape_tau, rate_tau)+
        dgamma(tau, delta_tau, delta_tau/prop_tau,log = T)-
        log_prior_lam_all(lam, alp, beta, tau)-
        log_prior_tau(tau, shape_tau, rate_tau)+
        dgamma(prop_tau, delta_tau, delta_tau/tau,log = T)
      if (log(runif(1)) <= log_tau_rate) {
        tau <- prop_tau
        acc_tau <- c(acc_tau, 1)
      } else {
        acc_tau <- c(acc_tau, 0)
      }
      tau_all <- c(tau_all, tau)
      
      #sampling lambda t
      for (t in 1:TT) {
        if (t==1) {
          lam_t_mius <-NA
          lam_t <- lam[t]
          lam_t_plus <- lam[t+1]
        } else if (t==TT) {
          lam_t_mius <-lam[t-1]
          lam_t <- lam[t]
          lam_t_plus <- NA
         
        } else {
          lam_t_mius <- lam[t-1]
          lam_t <- lam[t]
          lam_t_plus <- lam[t+1]
        }
        prop_lam <- rnorm(1, lam_t, sqrt(1/delta_lam))
        log_lam_rate <- loglik_at_t(prop_lam, U[[t]]) +
          log_prior_lam_t(lam_t_mius, prop_lam, lam_t_plus, alp, beta, tau) +
          dnorm(lam_t, prop_lam, sqrt(1/delta_lam), log=T)-
          loglik_at_t(lam_t, U[[t]]) -
          log_prior_lam_t(lam_t_mius, lam_t, lam_t_plus, alp, beta, tau) +
          dnorm(prop_lam, lam_t,sqrt(1/delta_lam), log=T)
        if (log(runif(1)) <= log_lam_rate) {
          lam[t] <- prop_lam
        } 
      }
      lam_all <- rbind(lam_all, lam)
  }
}


lam_all[,2]
emprical_cor <- 0
for (t in 1:TT) {
  emprical_cor[t] <- cor(U[[t]][,1], U[[t]][,2])
}
emprical_cor

est_cor <- 0
for (t in 1:TT) {
  est_cor[t] <- mean(lam2rho(lam_all[,t]))
}
est_cor
plot(1:TT, emprical_cor, type="l", col="red")
lines(1:TT, est_cor, col="blue")

range <- (burn_in*batch.size):(B*batch.size)
plot(alp_all[range], type="l")
plot(beta_all[range], type="l")
plot(tau_all[range], type="l")

mean(acc_alp[range])
mean(acc_beta[range])
mean(acc_tau[range])
for(t in 1:TT) {
  plot(lam_all[range,t], type="l")
}
range <- (burn_in*batch.size):(B*batch.size)

