library(mvtnorm)
library(sparr)
library(MCMCprecision)
library(foreach)
library(doParallel)
source("mcmc_helper.R")
source("LPML_DIC_WAIC.R")
args=(commandArgs(TRUE))
job_name=args[1]
job_num=as.numeric(args[2])
path=args[3]

set.seed(20231213)
U_train <- readRDS("../Data/PM_O3_2017_2019.rds")

#removing the boudary points
TT<- length(U_train)
nts <- c(1:TT)
for (day in c(1:TT)) {
  index <- apply(U_train[[day]],1, function(row) any(row %in% c(0, 1)))
  U_train[[day]] <- U_train[[day]][!index,]
  nts[[day]] <- nrow(U_train[[day]])
}



#mixtures
#dependence type

#
Parameters <- expand.grid(c(0,1,2,3,4,5),c(0,1,2,3,4), c(0,1,2))

ak <- 1
ct <- Parameters[job_num, 1]
MA <- Parameters[job_num, 2]
Season <- Parameters[job_num, 3]
c0 <- 1
p <- rep(0.25, 4)  #sum p =1
cts <- rep(ct, TT)

dep_type <- list(MA=MA,Season=c(12,Season))

burn_in=100
B=200
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
  ## adaptive , diminishing  
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
      eta_set <- eta_subset(TT, t, dep_type)
      tempPar <-  c0*p + colSums(rbind(Eta[eta_set,])) + colSums(Z[[t]])
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
          inv_set <- inv_subset(TT, t, dep_type)
          prop_denom_k <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(c0*p[k] + sum(prop_eta[eta_subset(TT, l, dep_type), k]))))))
          prop_denom_4 <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(c0*p[4] + sum(prop_eta[eta_subset(TT, l, dep_type), 4]))))))
          prop_dense <- (log(gamma(prop_etak+1))+prop_denom_k+log(gamma(prop_eta[t,4]+1))+prop_denom_4)
          
          current_denom_k <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(c0*p[k] + sum(current_eta[eta_subset(TT, l, dep_type), k]))))))
          current_denom_4 <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(c0*p[4] + sum(current_eta[eta_subset(TT, l, dep_type), 4]))))))
          current_dense <- (log(gamma(current_etak+1))+current_denom_k+log(gamma(current_eta[t,4]+1))+current_denom_4)
          
          
          
          ratio <-(prop_etak-current_etak)*(log(omega[k])+sum(log(Pi[j, inv_set ,k]))-log(omega[4])-sum(log(Pi[j, inv_set,4])))+
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
range <- seq(burn_in*batch.size, B*batch.size, 5)
#obtain DIC and LPML, WAIC

#WAIC
WAIC <- WAIC(U_train, Theta[range,,], Pi[range,,])
#LPML
LPML <- LPML(U_train, Theta[range,,], Pi[range,,])
#DIC
DIC <- DIC(U_train, Theta[range,,], Pi[range,,])


res <- list(nts=nts, U_train=U_train, Theta=Theta, Pi=Pi, Eta.all=Eta.all, omega.all=omega.all,
            beta.all=beta.all, kappa.all=kappa.all, acc.all=acc.all,
            WAIC=WAIC, DIC=DIC, LPML=LPML)

filename<- paste0("job_name=", job_name,"job_num=", job_num, 
                  "LPS_mix_new",".rds")

saveRDS(res, paste0(path,"/",filename))