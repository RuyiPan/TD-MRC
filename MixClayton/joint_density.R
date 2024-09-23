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
#estimated joint distribution form the best fitting
# ll.mg <- function(U, Theta) {
#   
#   U1 <- c(U[1], 1- U[1], 1- U[1], U[1])
#   U2 <- c(U[2], U[2], 1- U[2], 1- U[2])
#   
#   
#   ll.4 <- (log(Theta+1)-(Theta+1)*log(U1*U2)-(1/Theta+2)*log(U1^{-Theta}+U2^{-Theta}-1))
#   ll.4 <- ifelse(U1*U2 == 0, -Inf, ll.4)
#   return (ll.4)
# }

temp_res <- readRDS(paste0("../Fitting_Result/job_name=MixClayton_fitting_realjob_num=",68,"LPS_mix_new.rds"))

burn.in <- 100
B <- 200
batch.size <- 50
range <- seq(burn.in*batch.size,B*batch.size,10)
u1 <- seq(0, 1, length.out=100)
u2 <- seq(0, 1, length.out=100)
points <- expand.grid(u1,u2)
TT <- 36
den_all <- vector('list', TT)
for (t in 1:TT) {
  emp_den <- NULL
  for (i in (500*job_num+1):((job_num+1)*500)) {
    u <- c(points[i,1],points[i,2])
    den <- 0
    for (m in range) {
      den <- den + sum(exp(ll.mg(u,temp_res$Theta[m,t,]))*temp_res$Pi[m,t,])
    }
    emp_den <- c(emp_den, den/length(range))
  }
  den_all[[t]] <- emp_den
}
filename<- paste0("job_name=", job_name,"job_num=", job_num, 
                  "estimated_joint_dist",".rds")

res <- list(den_all=den_all, points=points[(500*job_num+1):((job_num+1)*500),])
saveRDS(res, paste0(path,"/",filename))
