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


LPML <- function(U_all, lam_all) {
  rho_all <- apply(lam_all, c(1,2), lam2rho)
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



WAIC <- function(U_train, lam_all) {
  TT <- length(U_train)
  M <- dim(lam_all)[1]
  fit <- 0
  penalty <- 0
  for (t in 1:TT) {
    lam_t <- lam_all[, t]
    dist_t <- NULL
    for (k in 1:M) {
      rho_t <- lam2rho(lam_t[k])
      temp_dist <- apply(U_train[[t]], 1, function(row) dGaussian(rho_t, row))
      dist_t <- cbind(dist_t, temp_dist)
    }
    fit <- fit + sum(log(rowMeans(dist_t)))
    penalty <-  penalty+sum(apply(log(dist_t), 1, var))
  }
  WAIC <- -2*fit + 2*penalty
  return (list(WAIC=WAIC, fit=fit, penalty=penalty))
}


