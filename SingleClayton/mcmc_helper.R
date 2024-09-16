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
eta_subset(36, 25, type=list(MA=2, Season=c(12,0)))
type<-list(MA=2, Season=c(12,2))


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



#mixtures
rClayton_cond1 <- function(u, theta, weights) {
  # comp <- sample(1:4,1,prob=weights)
  # th <- theta[comp]
  p<-u
  q<-runif(1)
  
  v1<-(p^(-theta[1])*(q^(-theta[1]/(theta[1]+1))-1)+1)^(-1/theta[1])
  
  v4<-1-(p^(-theta[2])*((1-q)^(-theta[2]/(theta[2]+1))-1)+1)^(-1/theta[2])
  
  v3<-1-((1-p)^(-theta[3])*((1-q)^(-theta[3]/(theta[3]+1))-1)+1)^(-1/theta[3])
  
  v2<-((1-p)^(-theta[4])*(q^(-theta[4]/(theta[4]+1))-1)+1)^(-1/theta[4])
  
  sum(c(v1, v2, v3, v4)*weights)
}


rClayton_cond2 <- function(u, theta, weights) {
  ty <- sample(1:4,1,prob=weights)
  th <- theta[ty]
  p<-u
  q<-runif(1)
  if (ty==1){
    v<-(p^(-th)*(q^(-th/(th+1))-1)+1)^(-1/th)
  }
  if (ty==4){
    v<-1-(p^(-th)*((1-q)^(-th/(th+1))-1)+1)^(-1/th)
  }
  if (ty==3){
    v<-1-((1-p)^(-th)*((1-q)^(-th/(th+1))-1)+1)^(-1/th)
  }
  if (ty==2){
    v<-((1-p)^(-th)*(q^(-th/(th+1))-1)+1)^(-1/th)
  }
  v
}

l.mg <- function(U, Theta) {
  
  U1 <- c(U[1], 1- U[1], 1- U[1], U[1])
  U2 <- c(U[2], U[2], 1- U[2], 1- U[2])
  l <- ((Theta+1)*(U1*U2)^{-(Theta+1)}*(U1^{-Theta}+U2^{-Theta}-1)^{-(1/Theta+2)})
  index <- which(Theta < 10^{-10}) #when theta is small enough, it's near indepedent
  l[index] <- 1
  return (l)
}

ll.mg <- function(U, Theta) {
  
  U1 <- c(U[1], 1- U[1], 1- U[1], U[1])
  U2 <- c(U[2], U[2], 1- U[2], 1- U[2])
  
  
  ll.4 <- (log(Theta+1)-(Theta+1)*log(U1*U2)-(1/Theta+2)*log(U1^{-Theta}+U2^{-Theta}-1))
  ll.4 <- ifelse(U1*U2 == 0, -Inf, ll.4)
  index <- which(Theta < 10^{-10}) #when theta is small enough, it's near indepedent
  ll.4[index] <- 0
  return (ll.4)
}


ll.mg2 <- function(U, Theta, Z) {
  
  ty <- which(Z==1)
  if (ty == 1) {
    U1 <- U[1]
    U2 <- U[2]
  } else if (ty==2) {
    U1 <- 1-U[1]
    U2 <- U[2]
  } else if (ty==3) {
    U1 <- 1-U[1]
    U2 <- 1-U[2]
  } else {
    U1 <- U[1]
    U2 <- 1-U[2]
  }
  theta <- Theta[ty]
  
  ll <- rep(0, 4)
  ll[ty] <- (log(theta+1)-(theta+1)*log(U1*U2)-(1/theta+2)*log(U1^{-theta}+U2^{-theta}-1))
  index <- which(theta < 10^{-10}) #when theta is small enough, it's near indepedent
  ll[index] <- 0
  return (ll)
}
