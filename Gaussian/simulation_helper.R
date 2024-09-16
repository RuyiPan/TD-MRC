vf<-function(ty,th,p,q){
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

simulate_mixC <- function(TT, nt) {
  U <- vector("list", length = TT)
  nts <- rep(300, TT)
  a.all <- matrix(0, nrow=TT, ncol=4)
  a<-c(0.4,0.25,0.1,0.25)
  a.all[1,]  <- a
  for (t in 1:TT) {
    L<-4
    ty<-c(1,2,3,4)
    th<-c(5,3,4,3)
    if (t != 1) {
      a[1] <- a[1]*0.95
      a[2] <- a[2]*1.05
      a[3] <- a[3]
      a[4] <- 1-sum(a[1:3])
      a.all[t,] <- a
    }
    print(a)
    
    nt <- nts[t]
    U[[t]] <- matrix(NA, nrow=nt, ncol=2)
    u<-1:nt
    v<-1:nt
    p<-matrix(NA,nrow=nt,ncol=L)
    q<-matrix(NA,nrow=nt,ncol=L)
    
    ua<-p
    va<-q
    for (l in 1:L){
      p[,l]<-runif(nt)
      q[,l]<-runif(nt)
      ua[,l]<-p[,l]
      va[,l]<-vf(ty[l],th[l],p[,l],q[,l])
    }
    z<-sample(1:L,nt,replace=TRUE,prob=a)
    for (i in 1:nt){
      for (l in 1:L){
        u[i]<-ua[i,z[i]]
        v[i]<-va[i,z[i]]
      }
    }
    
    U[[t]]<- cbind(u,v)
    
  }
  return (U)
}

