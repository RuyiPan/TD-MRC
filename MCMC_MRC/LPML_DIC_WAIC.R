#three gof measures
library(foreach)
library(doParallel)
library(copula)
l.mg <- function(U_trans, Theta) {
  comp_num <- nrow(U_trans)
  
  l_list <- vapply(1:comp_num, function(comp) l.clayton(U_trans[comp,], Theta[comp]),numeric(1))
  return (l_list)
}

ll.mg <- function(U_trans, Theta) {
  comp_num <- nrow(U_trans)
  
  ll_list <- vapply(1:comp_num, function(comp) ll.clayton(U_trans[comp,], Theta[comp]),numeric(1))
  return (ll_list)
}

##LPML
CPO.component <- function(U_trans, Theta, Pi) {
  ll_list <- l.mg(U_trans, Theta)
  
  comp <- sum(Pi*ll_list)
  
  return (1/comp)
}

CPO <- function(U_trans, Theta_all, Pi_all) {
  L <- dim(Theta_all)[1]
  1/mean(do.call(rbind,lapply(1:L, 
                              function(l) CPO.component(U_trans, Theta_all[l,], Pi_all[l,]))))
}


LPML <- function(U_all, Theta_all, Pi_all) {
  TT <- length(U_all)
  GlobalFunctions = ls(globalenv())
  ncores <- detectCores()
  cl <- parallel::makeCluster(ncores-2)
  doParallel::registerDoParallel(cl)
  
  res_CPO <- foreach(t=1:TT, .combine="c", .packages ="copula",
                     .export = GlobalFunctions)%dopar%{
                       nt <- dim(U_all[[t]])[1]
                       CPO_t <- c(1:nt)
                       for (i in 1:nt) {
                         CPO_t[i] <- CPO(U_all[[t]][i,,], Theta_all[,t,], Pi_all[,t,])
                       }
                       list(CPO_t)
                     }
  parallel::stopCluster(cl)
  
  LPML <- sum(log(unlist(res_CPO)))
  ALPML <- mean(log(unlist(res_CPO)))
  return(LPML)
}




ll.mg.U <- function(U_trans, Theta, Pi) {
  ll_list <- l.mg(U_trans, Theta)
  
  comp <- sum(Pi*ll_list)
  
  return (log(comp))
}


Devi <- function(U_all, Theta, Pi) {
  TT <- length(U_all)
  D.single <- 0
  for (t in 1:TT) {
    nt <- dim(U_all[[t]])[1]
    for (i in 1:nt) {
      D.single <-  D.single+ll.mg.U(U_all[[t]][i,,], Theta[t,], Pi[t,])
    }
  }
  
  return (-2*D.single)
}


#DIC
#U_all contains all rotated data
DIC <- function(U_all, Theta_all, Pi_all) {
  L <- dim(Theta_all)[1]
  GlobalFunctions = ls(globalenv())
  ncores <- detectCores()
  cl <- parallel::makeCluster(ncores-2)
  doParallel::registerDoParallel(cl)
  D_all <- foreach(l=1:L, .combine = "rbind",.packages ="copula", .export = GlobalFunctions)%dopar% {
    Devi(U_all, Theta_all[l,,], Pi_all[l,,])
  }
  parallel::stopCluster(cl)
  PMD <- mean(D_all)
  
  comp_num <- ncol(Theta_all[1,,])
  Theta_bar <- matrix(NA, nrow=TT, ncol=comp_num)
  Pi_bar <- matrix(NA, nrow=TT, ncol=comp_num)
  for(t in 1:TT){
    Theta_bar[t,] <- colMeans(Theta_all[,t,])
    Pi_bar[t,] <-  colMeans(Pi_all[,t,])
  }
  
  D_at_mean <- Devi(U_all, Theta_bar, Pi_bar)
  
  dic <- 2*PMD-D_at_mean
  return (dic)
}


#WAIC
#U_all contains all rotated data
WAIC <- function(U_all, Theta, Pi) {
  TT <- length(U_all)
  M <- dim(Theta)[1]
  fit <- 0
  penalty <-0
  for (t in 1:TT) {
    dist_t <- NULL
    for (k in 1:M) {
      temp_dist <-apply(U_all[[t]], 1, function(row) sum(l.mg(row,Theta[k,t,])*Pi[k,t,]))
      dist_t <- cbind(dist_t, temp_dist)
    }
    fit <- fit + sum(log(rowMeans(dist_t )))
    penalty <-  penalty+ sum(apply(log(dist_t), 1, var))
  }
  waic <- -2*fit + 2*penalty
  return (list(waic=waic, fit=fit, penalty=penalty))
}



