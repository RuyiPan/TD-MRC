library(mvtnorm)
library(sparr)
library(MCMCprecision)
library(foreach)
library(doParallel)
source("mcmc_helper.R")
args=(commandArgs(TRUE))
job_name=args[1]
job_num=as.numeric(args[2])
path=args[3]
job_num <- 1
t <- c(30:35)[job_num]
data  <- readRDS("../Data/PM_O3_2017_2019.rds")

set.seed(20240902)
t <- 31
#remove the data contain 1 or 0 because density equal to 0 on boundary
U_train <- data[1:t]
U_test <- data[t+1]

TT<- length(U_train)
nts <- c(1:TT)
for (day in c(1:TT)) {
  index <- apply(U_train[[day]],1, function(row) any(row %in% c(0, 1)))
  U_train[[day]] <- U_train[[day]][!index,]
  nts[[day]] <- nrow(U_train[[day]])
}

id_test <- apply(U_test[[1]],1, function(row) any(row %in% c(0, 1)))
U_test[[1]] <- U_test[[1]][!id_test ,]


#mixtures
ak <- 1
ct <- 10
q  <- 2
c0 <- 1
p <- rep(0.25, 4)  #sum p =1
cts <- rep(ct, TT)


burn_in=5
B=10
batch.size=50
M = B*batch.size


a<-c(0.4,0.25,0.1,0.25) 

#(a)
Z <- vector("list", length = TT)
for (t in 1:TT) {
  if (t==1) {
    Z[[t]] <- matrix(rmultinom(nrow(U_train[[t]]),size=1,prob=a), nrow=nrow(U_train[[t]]), ncol=4)
  } else {
    Z[[t]] <- matrix(NA, nrow=nrow(U_train[[t]]), ncol=4)
  }
  
}

# parameter for beta
ak <- ak
bk <- ak
a.beta <- c(1,1,1,1)*ak
b.beta <- c(1,1,1,1)*bk
beta <- rgamma(4, shape=1, rate=1)
beta.all <- matrix(0,nrow=M,ncol=4)
beta.all[1,] <- beta
omega <- rdirichlet(1, c0*p)
omega.all <- matrix(0,nrow=M,ncol=4)
omega.all[1,] <- omega

Pi <- array(dim=c(M, TT, 4))

Pi[1,,] <- rdirichlet(TT, a=c0*p)

d <- c(1,1,1,1) #parameter for theta
Theta <- array(dim=c(M, TT, 4))
Theta[1,,] <- rgamma(TT*4, shape=1, rate=1)

Eta <- t(rmultinom(TT, size=ct, c(0.25,0.25,0.25,0.25)))

Eta.all <- array(dim=c(M, TT, 4))
Eta.all[1,,]<- Eta
#record acceptance rate for each theta_tk
acc = matrix(0.3, nrow=TT, ncol=4) 
acc.all <- array(dim=c(B,TT,4))
ada.shape <- matrix(1, nrow=TT, ncol=4)  # initial value
kappa.all <- array(dim=c(B,TT,4))
C <- 1.01 #1.01,   1.1 (still fluctuate )
for (b in 1:B) {
  print(b)
  ## adaptive 
  if (b <=100) {
  ada.shape[acc < 0.3] <- ada.shape[acc < 0.3]*C^(sqrt(b))
  ada.shape[acc > 0.4] <- ada.shape[acc > 0.4]*C^(-sqrt(b))
  } 
  # print(acc)
  # print(ada.shape)
  kappa.all[b,,] <- ada.shape
  acc.all[b,,] <- acc
  count = matrix(0, nrow=TT, ncol=4)   #record acceptance for each batch
  for (it in 1:batch.size) {
    j=it+ batch.size*(b-1) #the jth iteration
    if (j == 1) next
    ll.olds <- matrix(0, nrow=TT, ncol=4)
    
    # (a) Posterior conditional for Z_ti i=1,...,nt
    # can improve efficient by parallel
    for (t in 1:TT) {
      #record for (e) use
      for (i in 1:nts[t]) {
        #Obtain the  pi_t^*
        ll <- ll.mg(U_train[[t]][i,], Theta[j-1, t,])
        weights <- Pi[j-1,t,]*exp(ll)
        pi_star <- weights/sum(weights)
        Z[[t]][i,]<- rmultinom(1, size=1, prob=pi_star)
      }
    }
    
    
    # Posterior conditional for theta_tk (possible parallel)
    
    for (t in 1:TT) {
      current_theta <- Theta[j-1, t, ]
      prop_theta <- c(rgamma(1, ada.shape[t,1], ada.shape[t,1]/current_theta[1]),
                      rgamma(1, ada.shape[t,2], ada.shape[t,2]/current_theta[2]),
                      rgamma(1, ada.shape[t,3], ada.shape[t,3]/current_theta[3]),
                      rgamma(1, ada.shape[t,4], ada.shape[t,4]/current_theta[4]))
      
      ll.new <- colSums(do.call(rbind, lapply(1:nts[t], function(i) {
        ll.mg2(U_train[[t]][i,], prop_theta, Z[[t]][i,])
      })))+(d-1)*log(prop_theta)-beta*prop_theta
      
      ll.old <- colSums(do.call(rbind, lapply(1:nts[t], function(i) {
        ll.mg2(U_train[[t]][i,], current_theta,Z[[t]][i,])
      })))+(d-1)*log(current_theta)-beta*current_theta
      #can change to log
      g.old <- c(dgamma(current_theta[1],ada.shape[t,1], ada.shape[t,1]/prop_theta[1], log=T),
                 dgamma(current_theta[2],ada.shape[t,2], ada.shape[t,2]/prop_theta[2], log=T),
                 dgamma(current_theta[3],ada.shape[t,3], ada.shape[t,3]/prop_theta[3], log=T),
                 dgamma(current_theta[4],ada.shape[t,4], ada.shape[t,4]/prop_theta[4], log=T))
      g.new <-c(dgamma(prop_theta[1],ada.shape[t,1], ada.shape[t,1]/current_theta[1],  log=T),
                dgamma(prop_theta[2],ada.shape[t,2], ada.shape[t,2]/current_theta[2],log=T),
                dgamma(prop_theta[3],ada.shape[t,3], ada.shape[t,3]/current_theta[3],log=T),
                dgamma(prop_theta[4],ada.shape[t,4], ada.shape[t,4]/current_theta[4],log=T))
      
      rate <- ll.new+g.old-ll.old-g.new
      
      # rate <- ifelse(is.na(rate), 0, rate)
      v_theta <- log(runif(4))
      
      count[t,] <- count[t,] +  (v_theta <= rate)
      prop_theta[!v_theta <= rate] <- current_theta[!v_theta <= rate]
      Theta[j, t, ] <- prop_theta
      
    }
    
    # (b) Posterior conditional for pi_t
    for (t in 1:TT) {
      tempPar <-  c0*p + colSums(rbind(Eta[max(1, t-q):t,])) + colSums(Z[[t]])
      Pi[j,t,] <- rdirichlet(1, tempPar)
    }
    
    
    
    # (c) Posterior conditional for eta_t
    for (t in 1:TT) {
      
      ##propose from RWM, uniform(eta - L, eta + L) , L = c_t/2
      L <- cts[t] / 2
      current_eta <- Eta
      for (k in 1:3) {
        if (cts[t] == 1) {
          prop_etak <- sample(c(1, 0), 
                              size=1)
        } else {
          prop_etak <- sample(c(as.integer(current_eta[t,k]-L):
                                  as.integer(current_eta[t,k]+L)), 
                              size=1)
        }
        
        if (prop_etak <= (cts[t]-sum(current_eta[t, 1:3][-k])) & prop_etak >=0 ) {
          current_etak <- current_eta[t,k]
          prop_eta <- current_eta 
          prop_eta[t, k] <- prop_etak 
          prop_eta[t, 4] <- cts[t]-sum(prop_eta[t, 1:3])
          prop_denom_k <-  sum(do.call(rbind, lapply(0:q, function(l) log(gamma(c0*p[k] + sum(prop_eta[max(t+l-q,1):min(t+l,TT), k]))))))
          prop_denom_4 <-  sum(do.call(rbind, lapply(0:q, function(l) log(gamma(c0*p[4] + sum(prop_eta[max(t+l-q,1):min(t+l,TT), 4]))))))
          prop_dense <- (log(gamma(prop_etak+1))+prop_denom_k+log(gamma(prop_eta[t,4]+1))+prop_denom_4)
          
          current_denom_k <-  sum(do.call(rbind, lapply(0:q, function(l) log(gamma(c0*p[k] + sum(current_eta[max(t+l-q,1):min(t+l,TT), k]))))))
          current_denom_4 <-  sum(do.call(rbind, lapply(0:q, function(l) log(gamma(c0*p[4] + sum(current_eta[max(t+l-q,1):min(t+l,TT), 4]))))))
          current_dense <- (log(gamma(current_etak+1))+current_denom_k+log(gamma(current_eta[t,4]+1))+current_denom_4)
          
          
          
          ratio <-(prop_etak-current_etak)*(log(omega[k])+sum(log(Pi[j, t:min(TT,t+q) ,k]))-log(omega[4])-sum(log(Pi[j, t:min(TT,t+q),4])))+
            current_dense-prop_dense
          
          # print(ratio)
          v<-log(runif(1))
          
          if(v<=ratio){
            current_eta <- prop_eta
          }
          
        } 
        
      }
      Eta <- current_eta
      
    }
    Eta.all[j,,]<-Eta
    
    # (d) Posterior conditional for omega
    omega <- rdirichlet(1, c0*p +colSums(Eta))
    omega.all[j,] <- omega
    
    # (f) Posterior conditional for beta_k
    shapes <- a.beta+TT*d
    rates <- b.beta+colSums(Theta[j,,])
    beta <- c(rgamma(1, shape=shapes[1], rate=rates[1]),
              rgamma(1, shape=shapes[2], rate=rates[2]),
              rgamma(1, shape=shapes[3], rate=rates[3]),
              rgamma(1, shape=shapes[4], rate=rates[4]))
    beta.all[j,] <- beta
  }
  
  
  acc <- count/50
  
  
  
}
range <- (burn_in*batch.size):(B*batch.size)
predictive_dist <- NULL
theta <- NULL
weight <- NULL
eta_sum <- NULL
for (k in range) {
  theta_t_plus <- c(rgamma(1, d[1], beta.all[k,1]),
                    rgamma(1, d[2], beta.all[k,2]),
                    rgamma(1, d[3], beta.all[k,3]),
                    rgamma(1, d[4], beta.all[k,4]))
  theta <- rbind(theta, theta_t_plus)
  eta_t_plus <- as.vector(rmultinom(1, ct, omega.all[k,]))
  dir_w <- c0*p+colSums(Eta.all[k, (TT+1-q):TT, ])+eta_t_plus
  eta_sum <- rbind(eta_sum, dir_w)
  weights_t_plus <- rdirichlet(1, dir_w)
  weight <- rbind(weight, weights_t_plus)
  temp_dist <- apply(U_test[[1]], 1, function(row) sum(l.mg(row,theta_t_plus)*weights_t_plus))
  predictive_dist <- cbind(predictive_dist, temp_dist)
}
predictive_dist_mean <- rowMeans(predictive_dist)
LPS <- sum(log(predictive_dist_mean)) ##-9.996587
res <- list(nts=nts, U_train=U_train, U_test=U_test, Theta=Theta, Pi=Pi, Eta.all=Eta.all, omega.all=omega.all,
            beta.all=beta.all, kappa.all=kappa.all, acc.all=acc.all,
            predictive_dist=predictive_dist,
            LPS=LPS)

filename<- paste0("job_name=", job_name,"job_num=", job_num, 
                  "LPS_mix",".rds")

saveRDS(res, paste0(path,"/",filename))


