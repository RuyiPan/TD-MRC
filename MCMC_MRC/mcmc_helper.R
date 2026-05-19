library(copula)
#define a index dict for component
#@m:dimension for a copula
generate_index <- function(dim) {
  idx <- as.data.frame(expand.grid(rep(list(c(0, 1)), dim)))
  sorted_idx <- idx[do.call(order, idx), ]
  rownames(sorted_idx) <- NULL
  return(sorted_idx)
}
#dim<- 3
# 
# # Example 
# index_dict <- generate_index(dim)
# print(index_dict)
# 000 -- 1
# 001 -- 2
# 010 -- 3
# 011 -- 4
# 100 -- 5
# 101 -- 6
# 110 -- 7
# 111 -- 8
# Precompute U_transformed for all U and components
precompute_U_transformed <- function(U_train) {
  TT <- length(U_train)
  dim <- ncol(U_train[[1]])
  comp_num <- 2^dim
  comp_dict <- generate_index(dim)
  # Initialize a list to store U_transformed for each time step
  U_transformed_list <- vector("list", TT)
  
  # for (t in 1:TT) {
  #   nts_t <- nrow(U_train[[t]])
  #   # Precompute U_transformed for all samples and components
  #   U_transformed <- array(NA, dim = c(nts_t, comp_num, dim))
  #   for (comp in 1:comp_num) {
  #     U_transformed[, comp,] <- as.matrix(do.call(rbind,apply(U_train[[t]], 1, 
  #                                     function(row) row^(1 - comp_dict[comp, ]) * (1-row)^comp_dict[comp, ]))) 
  #   }
  #   U_transformed_list[[t]] <- U_transformed
  # }
  
  for (t in 1:TT) {
    nts_t <- nrow(U_train[[t]])
    U_transformed <- array(NA, dim = c(nts_t, comp_num, dim))
    for (day in 1:nts_t) {
      u <- U_train[[t]][day,]
      U_transformed[day,,] <- t(apply(comp_dict, 1, function(b) {
        b * (1 - u) + (1 - b) * u
      }))
    }
    U_transformed_list[[t]] <- U_transformed
  }
  
  return(U_transformed_list)
}

# Example usage:
# U_transformed <- precompute_U_transformed(U_train, comp_dict)

#find the subset of eta for weights
#@type: the dependence type, MA(q) or Seasonal S(p), or both
#@TT: the number time points
#@t: the current time
eta_subset <- function(TT, t, type=list(MA=q, 
                                        Season=c(s,p))) {
  partial_MA <- t:(t-type[['MA']])
  partial_Season <- t-type[['Season']][1]*(0:type[['Season']][2])
  candidate <- union(partial_MA,partial_Season)
  candidate[candidate%in%(1:TT)]
}

#example
# eta_subset(36, 25, type=list(MA=2, Season=c(12,0)))
# type<-list(MA=2, Season=c(12,2))


#Find the inverse subset, weight for eta
inv_subset <- function(TT, t, type=list(MA=q, 
                                        Season=c(s,p))) {
  partial_MA <- t:(t+type[['MA']])
  partial_Season <- t+type[['Season']][1]*(0:type[['Season']][2])
  candidate <- union(partial_MA,partial_Season)
  candidate[candidate%in%(1:TT)]
}
#example
# inv_subset(36, 13, type=list(MA=2, Season=c(12,2)))
# type<-list(MA=2, Season=c(12,2))


#likelihood, density for a specific component
# l.c <- function(U, theta, comp){
#   if (sum(U==0) | sum(U==1)) {
#     return (0)
#   }
#   if (theta < 10^{-10}) {
#     return (1)
#   }
#   U <- U^{1-comp_dict[comp,]}*(1-U)^comp_dict[comp,]
#   dim <- length(U)
#   S <- sum(U^{-theta})-dim + 1
#   l <- prod(1+1:(dim-1)*theta)*S^{-1/theta-dim}*prod(U^{-theta-1})
#   return (l)
# }

# #loglikelihood for a specific component
# ll.c <- function(U, theta, comp) {
#   if (sum(U==0) | sum(U==1)) {
#     return (-Inf)
#   }
#   if (theta < 10^{-10}) {
#     return (0)
#   }
#   U <- U^{1-comp_dict[comp,]}*(1-U)^comp_dict[comp,]
#   dim <- length(U)
#   S <- sum(U^{-theta}) - dim + 1
#   ll <- sum(log(1+1:(dim-1)*theta))+{-1/theta-dim}*log(S)+sum((-theta-1)*log(U))
#   return (ll)
# }

#likelihood for all 2^m components
# l.mg <- function(U, Theta) {
#   dim <- length(U)
#   comp_num <- 2^dim
#   
#   l_list <- do.call(rbind,lapply(1:comp_num, function(comp) l.c(U, Theta[j], comp)))
#   
#   return (l_list)
# }

# library(copula)
# theta <- 10^{-10}
# dim = 4
# clayton_cop <- claytonCopula(param = theta, dim = dim)  # Define the Clayton copula
# U <- matrix(c(0.9, 0.7, 0.9,0.4), nrow = 1)               # Input as a matrix
# density <- dCopula(U,clayton_cop,log = T)               # Compute the density
# print(density)
# 
# comp_dict <- generate_index(dim)
# comp=1
# ll.c(U, theta, comp=comp)
# U^{1-comp_dict[comp,]}*(1-U)^comp_dict[comp,]

#loglikelihood for a specific component
# ll.c <- function(U, theta, comp) {
#   if (sum(U==0) | sum(U==1)) {
#     return (-Inf)
#   }
#   if (theta < 10^{-10}) {
#     return (0)
#   }
#   U <- U^{1-comp_dict[comp,]}*(1-U)^comp_dict[comp,]
#   dim <- length(U)
#   S <- sum(U^{-theta}) - dim + 1
#   ll <- sum(log(1+1:(dim-1)*theta))+{-1/theta-dim}*log(S)+sum((-theta-1)*log(U))
#   return (ll)
# }

l.clayton <- function(U, theta){
  if (sum(U==0) | sum(U==1)) {
    return (0)
  }
  if (theta < 10^{-10}) {
    return (1)
  }
  dim <- length(U)
  S <- sum(U^{-theta})-dim + 1
  l <- prod(1+1:(dim-1)*theta)*S^{-1/theta-dim}*prod(U^{-theta-1})
  if (is.na(l)) {
    clayton_cop <- claytonCopula(param = theta, dim = dim)  # Define the Clayton copula
    # Input as a matrix
    l <- dCopula(U,clayton_cop,log = F) 
  }
  return (l)
}


ll.clayton  <- function(U, theta) {
  if (sum(U==0) | sum(U==1)) {
    return (-Inf)
  }
  if (theta < 10^{-10}) {
    return (0)
  }
  dim <- length(U)
  S <- sum(U^{-theta}) - dim + 1
  ll <- sum(log(1+1:(dim-1)*theta))+{-1/theta-dim}*log(S)+sum((-theta-1)*log(U))
  if (is.na(ll)) {
    clayton_cop <- claytonCopula(param = theta, dim = dim)  # Define the Clayton copula
    # Input as a matrix
    ll <- dCopula(U,clayton_cop,log = T) 
  }
  return (ll)
}


l.mg <- function(U_trans, Theta) {
  comp_num <- nrow(U_trans)
  
  l_list <- vapply(1:comp_num, function(comp) l.clayton(U_trans[comp,], Theta[comp]),numeric(1))
  return (l_list)
}


#likelihood for all 2^m components
ll.mg <- function(U_trans, Theta) {
  comp_num <- nrow(U_trans)
  
  ll_list <- vapply(1:comp_num, function(comp) ll.clayton(U_trans[comp,], Theta[comp]),numeric(1))
  return (ll_list)
}

# ll.mg_dim2 <- function(U, Theta) {
#   
#   U1 <- c(U[1], U[1], 1- U[1], 1-U[1])
#   U2 <- c(U[2], 1-U[2], U[2], 1- U[2])
#   
#   
#   ll.4 <- (log(Theta+1)-(Theta+1)*log(U1*U2)-(1/Theta+2)*log(U1^{-Theta}+U2^{-Theta}-1))
#   ll.4 <- ifelse(U1*U2 == 0, -Inf, ll.4)
#   index <- which(Theta < 10^{-10}) #when theta is small enough, it's near indepedent
#   ll.4[index] <- 0
#   return (ll.4)
# }

ll.mg2 <- function(U_trans, Theta, Z) {
  comp <- which(Z==1)
  U <- U_trans[comp,]
  dim <- length(U)
  theta <- Theta[comp]
  
  ll.one <- ll.clayton(U,theta)
  
  ll_list <- rep(0, comp_num)
  ll_list[comp] <- ll.one
  return (ll_list)
}



# ll.mg2_dim2 <- function(U, Theta, Z) {
#   
#   ty <- which(Z==1)
#   if (ty == 1) {
#     U1 <- U[1]
#     U2 <- U[2]
#   } else if (ty==2) {
#     U1 <- U[1]
#     U2 <- 1-U[2]
#   } else if (ty==3) {
#     U1 <- 1-U[1]
#     U2 <- U[2]
#   } else {
#     U1 <- 1-U[1]
#     U2 <- 1-U[2]
#   }
#   theta <- Theta[ty]
#   
#   ll <- rep(0, 4)
#   ll[ty] <- (log(theta+1)-(theta+1)*log(U1*U2)-(1/theta+2)*log(U1^{-theta}+U2^{-theta}-1))
#   index <- which(theta < 10^{-10}) #when theta is small enough, it's near indepedent
#   ll[index] <- 0
#   return (ll)
# }
