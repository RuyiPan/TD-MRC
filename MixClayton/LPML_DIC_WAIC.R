#three gof measures
library(foreach)
library(doParallel)
##LPML
CPO.component <- function(U, Theta, Pi) {
  U1 <- c(U[1], 1- U[1], 1- U[1], U[1])
  U2 <- c(U[2], U[2], 1- U[2], 1- U[2])
  l.4 <- exp(-(Theta+1)*log((U1*U2))-(1/Theta+2)*log(U1^{-Theta}+U2^{-Theta}-1)+log(Theta+1))

  l.4 <- ifelse(U1*U2 == 0, 0, l.4)
  
  comp <- sum(Pi*l.4)
  
  return (1/comp)
}

CPO <- function(U, Theta_all, Pi_all) {
  L <- dim(Theta_all)[1]
  1/mean(do.call(rbind,lapply(1:L, 
                function(l) CPO.component(U, Theta_all[l,], Pi_all[l,]))))
}


LPML <- function(U_all, Theta_all, Pi_all) {
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
                         CPO_t[i] <- CPO(U_all[[t]][i,], Theta_all[,t,], Pi_all[,t,])
                       }
                       list(CPO_t)
                     }
  parallel::stopCluster(cl)
  
  LPML <- sum(log(unlist(res_CPO)))
  ALPML <- mean(log(unlist(res_CPO)))
  return(LPML)
}




ll.mg.U <- function(U, Theta, Pi) {
  U1 <- c(U[1], 1- U[1], 1- U[1], U[1])
  U2 <- c(U[2], U[2], 1- U[2], 1- U[2])
  
  l.4 <- exp(-(Theta+1)*log((U1*U2))-(1/Theta+2)*log(U1^{-Theta}+U2^{-Theta}-1)+log(Theta+1))
  
  l.4 <- ifelse(U1*U2 == 0, 0, l.4)
  
  comp <- sum(Pi*l.4)
  
  return (log(comp))
}


Devi <- function(U_all, Theta, Pi) {
  TT <- length(U_all)
  D.single <- 0
  for (t in 1:TT) {
    nt <- nrow(U_all[[t]])
    for (i in 1:nt) {
      D.single <-  D.single+ll.mg.U(U_all[[t]][i,], Theta[t,], Pi[t,])
    }
  }
  
  return (-2*D.single)
}


#DIC
DIC <- function(U_all, Theta_all, Pi_all) {
  L <- dim(Theta_all)[1]
  GlobalFunctions = ls(globalenv())
  ncores <- detectCores()
  cl <- parallel::makeCluster(ncores-2)
  doParallel::registerDoParallel(cl)
  D_all <- foreach(l=1:L, .combine = "rbind",.export = GlobalFunctions)%dopar% {
    Devi(U_all, Theta_all[l,,], Pi_all[l,,])
  }
  parallel::stopCluster(cl)
  PMD <- mean(D_all)
  
  Theta_bar <- matrix(NA, nrow=TT, ncol=4)
  Pi_bar <- matrix(NA, nrow=TT, ncol=4)
  for(t in 1:TT){
    Theta_bar[t,] <- colMeans(Theta_all[,t,])
    Pi_bar[t,] <-  colMeans(Pi_all[,t,])
  }
  
  D_at_mean <- Devi(U_all, Theta_bar, Pi_bar)
  
  DIC <- 2*PMD-D_at_mean
  return (DIC)
}


#WAIC
WAIC <- function(U_train, Theta, Pi) {
  TT <- length(U_train)
  M <- dim(Theta)[1]
  fit <- 0
  penalty <-0
  for (t in 1:TT) {
    dist_t <- NULL
    for (k in 1:M) {
      temp_dist <-apply(U_train[[t]], 1, function(row) sum(l.mg(row,Theta[k,t,])*Pi[k,t,]))
      dist_t <- cbind(dist_t, temp_dist)
    }
    fit <- fit + sum(log(rowMeans(dist_t )))
    penalty <-  penalty+ sum(apply(log(dist_t), 1, var))
  }
  WAIC <- -2*fit + 2*penalty
  return (list(WAIC=WAIC, fit=fit, penalty=penalty))
}


