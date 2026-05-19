library(mvtnorm)
library(sparr)
library(MCMCprecision)
library(foreach)
library(doParallel)
source("mcmc_helper.R")
source("simulate_data.R")
args=(commandArgs(TRUE))
job_name=args[1]
job_num=as.numeric(args[2])
path=args[3]
MA=as.numeric(args[4])
Season=as.numeric(args[5])


set.seed(20231213)
#get train and test data
U_train <- data$train
U_test <- data$test

#sample a data
U_train <- simulate_mixC_dim2(20, 300)
# U_train <- simulate_mixC_dim3(13,200)
dim <- ncol(U_train[[1]])
comp_num <- 2^dim




TT<- length(U_train)
nts <- c(1:TT)
for (day in c(1:TT)) {
  index <- apply(U_train[[day]],1, function(row) any(row %in% c(0, 1)))
  U_train[[day]] <- U_train[[day]][!index,]
  nts[[day]] <- nrow(U_train[[day]])
}

id_test <- apply(U_test[[1]],1, function(row) any(row %in% c(0, 1)))
U_test[[1]] <- U_test[[1]][!id_test ,]

U_train_transformed_list <- precompute_U_transformed(U_train) #precompute the rotated data for rotated copula
U_test_transformed_list <- precompute_U_transformed(U_test)

#mixtures
#dependence type
MA <- 3; Season <- 0
dep_type <- list(MA=MA,Season=c(12,Season))
eta_subsets <- vector("list", length = TT) # Precompute inv_set once for the current t
eta_inv_subsets <- vector("list", length = TT)
for (t in 1:TT) {
  eta_subsets[[t]] <- eta_subset(TT, t, dep_type)
  eta_inv_subsets[[t]] <- inv_subset(TT, t, dep_type)
}


#
ek <- 1
at <- 30
a0 <- 1
p <- rep(1/comp_num, comp_num)  #sum p =1
ats <- rep(at, TT)

burn_in=100
B=160
batch.size=50
M = B*batch.size



#(a)
Z <- vector("list", length = TT)
for (t in 1:TT) {
    Z[[t]] <- matrix(rmultinom(nrow(U_train[[t]]),size=1,prob=rep(1/comp_num, comp_num)), nrow=nrow(U_train[[t]]), ncol=comp_num)
}

# parameter for beta
ek <- ek
gk <- ek
e.beta <- rep(1, comp_num)*ek
g.beta <- rep(1, comp_num)*gk
beta <- rgamma(comp_num, shape=1, rate=1)
beta.all <- matrix(0, nrow=M, ncol=comp_num)
beta.all[1,] <- beta

omega <- rdirichlet(1, a0*p)
omega.all <- matrix(0,nrow=M,ncol=comp_num)
omega.all[1,] <- omega

Pi <- array(dim=c(M, TT, comp_num))
Pi[1,,] <- rdirichlet(TT, a=a0*p)


d <- rep(1,comp_num) #parameter for theta
Theta <- array(dim=c(M, TT, comp_num))
Theta[1,,] <- rgamma(TT*comp_num, shape=1, rate=1)


Eta <- t(rmultinom(TT, size=at, rep(1/comp_num,comp_num)))
Eta.all <- array(dim=c(M, TT, comp_num))
Eta.all[1,,]<- Eta
#record acceptance rate for each theta_tk
acc = matrix(0.3, nrow=TT, ncol=comp_num) 
acc.all <- array(dim=c(B,TT,comp_num))
ada.shape <- matrix(1, nrow=TT, ncol=comp_num)  # initial value
kappa.all <- array(dim=c(B,TT,comp_num))
C <- 1.01 #1.01,   1.1 (still fluctuate )


for (b in 1:B) {
  print(b)
  ## adaptive , diminishing  
  if (b <=burn_in) {
    ada.shape[acc < 0.3] <- ada.shape[acc < 0.3]*C^(sqrt(b))
    ada.shape[acc > 0.4] <- ada.shape[acc > 0.4]*C^(-sqrt(b))
  } 
  # print(acc)
  # print(ada.shape)
  kappa.all[b,,] <- ada.shape
  acc.all[b,,] <- acc
  count = matrix(0, nrow=TT, ncol=comp_num)   #record acceptance for each batch
  for (it in 1:batch.size) {
    j=it+ batch.size*(b-1) #the jth iteration
    if (j == 1) next
    ll.olds <- matrix(0, nrow=TT, ncol=comp_num)
    
    # (a) Posterior conditional for Z_ti i=1,...,nt
    # can improve efficient by parallel
    for (t in 1:TT) {
      #record for (e) use
      for (i in 1:nts[t]) {
        #Obtain the  pi_t^*
        ll <- ll.mg(U_train_transformed_list[[t]][i, , ], Theta[j-1, t,])
        weights <- Pi[j-1,t,]*exp(ll)
        pi_star <- weights/sum(weights)
        Z[[t]][i,]<- rmultinom(1, size=1, prob=pi_star)
      }
    }
    
    
    # Posterior conditional for theta_tk (possible parallel)
     
    for (t in 1:TT) {
      current_theta <- Theta[j-1, t, ]
      
      # prop_theta <- c(rgamma(1, ada.shape[t,1], ada.shape[t,1]/current_theta[1]),
      #                 rgamma(1, ada.shape[t,2], ada.shape[t,2]/current_theta[2]),
      #                 rgamma(1, ada.shape[t,3], ada.shape[t,3]/current_theta[3]),
      #                 rgamma(1, ada.shape[t,4], ada.shape[t,4]/current_theta[4]))
      
      # Generate prop_theta using vectorized rgamma
      prop_theta <- rgamma(comp_num, ada.shape[t, ], ada.shape[t, ] / current_theta)
      
      
      # ll.new <- colSums(do.call(rbind, lapply(1:nts[t], function(i) {
      #   ll.mg2(U_train_transformed_list[[t]][i,,], prop_theta, Z[[t]][i,])
      # })))+(d-1)*log(prop_theta)-beta*prop_theta
      # 
      # ll.old <- colSums(do.call(rbind, lapply(1:nts[t], function(i) {
      #   ll.mg2(U_train_transformed_list[[t]][i,,], current_theta, Z[[t]][i,])
      # })))+(d-1)*log(current_theta)-beta*current_theta
      
      # ll.new_dim2 <- colSums(do.call(rbind, lapply(1:nts[t], function(i) {
      #   ll.mg2_dim2(U_train[[t]][i,], prop_theta, Z[[t]][i,])
      # })))+(d-1)*log(prop_theta)-beta*prop_theta
      
      # Compute ll.new and ll.old using vectorized operations
      ll.new <- rowSums(sapply(1:nts[t], function(i) {
        ll.mg2(U_train_transformed_list[[t]][i, , ], prop_theta, Z[[t]][i, ])
      })) + (d - 1) * log(prop_theta) - beta * prop_theta
      # ll.new_dim2;ll.new
      ll.old <- rowSums(sapply(1:nts[t], function(i) {
        ll.mg2(U_train_transformed_list[[t]][i, , ], current_theta, Z[[t]][i, ])
      })) + (d - 1) * log(current_theta) - beta * current_theta
      
      #can change to log
      # g.old <- c(dgamma(current_theta[1],ada.shape[t,1], ada.shape[t,1]/prop_theta[1], log=T),
      #            dgamma(current_theta[2],ada.shape[t,2], ada.shape[t,2]/prop_theta[2], log=T),
      #            dgamma(current_theta[3],ada.shape[t,3], ada.shape[t,3]/prop_theta[3], log=T),
      #            dgamma(current_theta[4],ada.shape[t,4], ada.shape[t,4]/prop_theta[4], log=T))
      # g.new <-c(dgamma(prop_theta[1],ada.shape[t,1], ada.shape[t,1]/current_theta[1],  log=T),
      #           dgamma(prop_theta[2],ada.shape[t,2], ada.shape[t,2]/current_theta[2],log=T),
      #           dgamma(prop_theta[3],ada.shape[t,3], ada.shape[t,3]/current_theta[3],log=T),
      #           dgamma(prop_theta[4],ada.shape[t,4], ada.shape[t,4]/current_theta[4],log=T))
      
      # Compute g.old and g.new using vectorized dgamma
      g.old <- dgamma(current_theta, ada.shape[t, ], ada.shape[t, ] / prop_theta, log = TRUE)
      g.new <- dgamma(prop_theta, ada.shape[t, ], ada.shape[t, ] / current_theta, log = TRUE)
      
      rate <- ll.new+g.old-ll.old-g.new
      
      # rate <- ifelse(is.na(rate), 0, rate)
      v_theta <- log(runif(comp_num))
      
      count[t,] <- count[t,] +  (v_theta <= rate)
      prop_theta[!v_theta <= rate] <- current_theta[!v_theta <= rate]
      Theta[j, t, ] <- prop_theta
      
    }
    
    # (b) Posterior conditional for pi_t
    for (t in 1:TT) {
      # eta_set <- eta_subset(TT, t, dep_type)
      eta_set <- eta_subsets[[t]]
      tempPar <-  a0*p + colSums(rbind(Eta[eta_set,])) + colSums(Z[[t]])
      Pi[j,t,] <- rdirichlet(1, tempPar)
    }
    
    
    
    # (c) Posterior conditional for eta_t
    for (t in 1:TT) {
      
      ##propose from RWM, uniform(eta - L, eta + L) , L = c_t/2
      L <- ats[t] / 2
      current_eta <- Eta
      # inv_set <- inv_subset(TT, t, dep_type)
      inv_set <- eta_inv_subsets[[t]]
      for (k in 1:(comp_num-1)) {
        if (ats[t] == 1) {
          prop_etak <- sample(c(1, 0), 
                              size=1)
        } else {
          prop_etak <- sample(c(as.integer(current_eta[t,k]-L):
                                  as.integer(current_eta[t,k]+L)), 
                              size=1)
        }
        
        if (prop_etak <= (ats[t]-sum(current_eta[t, 1:(comp_num-1)])) & prop_etak >=0 ) {
          current_etak <- current_eta[t,k]
          prop_eta <- current_eta 
          prop_eta[t, k] <- prop_etak 
          prop_eta[t, comp_num] <- ats[t]-sum(prop_eta[t, 1:(comp_num-1)])
          
          prop_denom_k <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(a0*p[k] + sum(prop_eta[eta_subsets[[l]], k]))))))
          prop_denom_last <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(a0*p[comp_num] + sum(prop_eta[eta_subsets[[l]], comp_num]))))))
          prop_dense <- (log(gamma(prop_etak+1))+prop_denom_k+log(gamma(prop_eta[t,comp_num]+1))+prop_denom_last)
          
          current_denom_k <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(a0*p[k] + sum(current_eta[eta_subsets[[l]], k]))))))
          current_denom_last <-  sum(do.call(rbind, lapply(inv_set, function(l) log(gamma(a0*p[comp_num] + sum(current_eta[eta_subsets[[l]], comp_num]))))))
          current_dense <- (log(gamma(current_etak+1))+current_denom_k+log(gamma(current_eta[t,comp_num]+1))+current_denom_last)
          
          
          
          ratio <-(prop_etak-current_etak)*(log(omega[k])+sum(log(Pi[j, inv_set,k]))-log(omega[comp_num])-sum(log(Pi[j, inv_set,comp_num])))+
            current_dense-prop_dense
          
          # print(ratio)
          
          if(log(runif(1))<=ratio){
            current_eta <- prop_eta
          }
          
        } 
        
      }
      Eta <- current_eta
      
    }
    Eta.all[j,,]<-Eta
    
    
    # (d) Posterior conditional for omega
    omega <- rdirichlet(1, a0*p +colSums(Eta))
    omega.all[j,] <- omega
    
    # (f) Posterior conditional for beta_k
    shapes <- e.beta+TT*d
    rates <- g.beta+colSums(Theta[j,,])
    # beta <- c(rgamma(1, shape=shapes[1], rate=rates[1]),
    #           rgamma(1, shape=shapes[2], rate=rates[2]),
    #           rgamma(1, shape=shapes[3], rate=rates[3]),
    #           rgamma(1, shape=shapes[4], rate=rates[4]))
    beta <- rgamma(comp_num, shape = shapes, rate = rates)
    beta.all[j,] <- beta
  }
  
  
  acc <- count/50
  
  
  
}

range <- seq(burn_in*batch.size, B*batch.size, 10)
#obtain DIC and LPML, WAIC


source("LPML_DIC_WAIC.R")
#WAIC
WAIC.value <- WAIC(U_train_transformed_list, Theta[range,,], Pi[range,,])
#LPML
LPML.value <- LPML(U_train_transformed_list, Theta[range,,], Pi[range,,])
#DIC
DIC.value <- DIC(U_train_transformed_list, Theta[range,,], Pi[range,,])


res <- list(nts=nts, U_train=U_train, U_test=U_test, Theta=Theta, Pi=Pi, Eta.all=Eta.all, omega.all=omega.all,
            beta.all=beta.all, kappa.all=kappa.all, acc.all=acc.all,
            WAIC=WAIC, DIC=DIC, LPML=LPML,
            MSE=MSE, MSEt=MSEt)



# dim(acc.all)
# plot(acc.all[,1,2],type="l")
# plot(Theta[range,13,5], type="l")
# range <- seq(burn_in*batch.size,B*batch.size)
# apply(Theta[range,13,], 2, function(x) c(mean(x), quantile(x, probs=c(0.025, 0.975))))
# apply(Pi[range,12,], 2, function(x) c(mean(x), quantile(x, probs=c(0.025, 0.975))))


#get LPS
predictive_dist <- NULL
eta_set <- eta_subset(TT, TT+1, dep_type)
for (k in range) {
  # theta_t_plus <- c(max(10^{-10},rgamma(1, d[1], beta.all[k,1])),
  #                   max(10^{-10},rgamma(1, d[2], beta.all[k,2])),
  #                   max(10^{-10},rgamma(1, d[3], beta.all[k,3])),
  #                   max(10^{-10},rgamma(1, d[4], beta.all[k,4])))
  theta_t_plus <- pmax(10^{-10},rgamma(comp_num, d, beta.all[k,]))
  eta_t_plus <- as.vector(rmultinom(1, at, omega.all[k,]))
  
  dir_w <- a0*p+colSums(rbind(Eta.all[k, eta_set, ]))+eta_t_plus
  weights_t_plus <- rdirichlet(1, dir_w)
  temp_dist <- apply(U_test_transformed_list[[1]], 1, function(row) sum(l.mg(row,theta_t_plus)*weights_t_plus))
  predictive_dist <- cbind(predictive_dist, temp_dist)
}

predictive_dist_mean <- rowMeans(predictive_dist)
LPS.value <- sum(log(predictive_dist_mean))

temp_res <- list(nts=nts,
                 U_train=U_train, U_train_transformed_list=U_train_transformed_list,
                 U_test=U_test, U_test_transformed_list=U_test_transformed_list,
                 Theta=Theta, Pi=Pi, Eta.all=Eta.all, omega.all=omega.all,
                 beta.all=beta.all, kappa.all=kappa.all, acc.all=acc.all,
                 WAIC=WAIC.value, DIC=DIC.value, LPML=LPML.value, 
                 predictive_dist=predictive_dist, LPS=LPS)
