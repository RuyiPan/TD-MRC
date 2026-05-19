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
temp_res <- readRDS("../MCMC_MRC/Result/fitting.rds")
num_points <- 1000
burn.in <- 100
B <- 200
batch.size <- 50
range <- seq(burn.in*batch.size,B*batch.size,10)
u1 <- seq(0, 1, length.out=50)
u2 <- seq(0, 1, length.out=50)
u3 <- seq(0, 1, length.out=50)
points <- as.matrix(expand.grid(u1,u2,u3))
data <- points[(num_points*job_num+1):((job_num+1)*num_points),]
precompute_U_transformed_2 <- function(data) {
  dim <- ncol(data)
  comp_num <- 2^dim
  comp_dict <- generate_index(dim)
  # Initialize a list to store U_transformed for each time step
  nts_t <- nrow(data)
  U_transformed <- array(NA, dim = c(nts_t, comp_num, dim))
    for (day in 1:nts_t) {
      u <- data[day,]
      U_transformed[day,,] <- t(apply(comp_dict, 1, function(b) {
        b * (1 - u) + (1 - b) * u
      }))
    }
  return(U_transformed)
}

U_transformed <- precompute_U_transformed_2(data)
TT <- length(temp_res$nts)
den_all <- vector('list', TT)
range_length <- length(range)
for (t in 1:TT) {
  emp_den <- numeric(num_points)  # Preallocate for efficiency
  for (i in 1:num_points) {
    # Vectorized computation of `den` across `range`
    den <- sapply(range, function(m) {
      sum(l.mg(U_transformed[i,,], temp_res$Theta[m, t, ]) * temp_res$Pi[m, t, ])
    })
    # Average over the range
    emp_den[i] <- sum(den) / range_length
  }
  den_all[[t]] <- emp_den
}
filename<- paste0("job_name=", job_name,"job_num=", job_num, 
                  "estimated_joint_dist",".rds")

res <- list(den_all=den_all, points=data)
saveRDS(res, paste0(path,"/",filename))