# *************************************************************************
# CMI2NI: Conditional mutual inclusive information(CMI2)-based Network
# Inference method from gene expression data
# *************************************************************************
# This is matlab code for netwrk inference method CMI2NI.
# Input:
# 'data' is expression of variable,in which row is varible and column is the sample;
# 'lamda' is the parameter decide the dependence;
# 'order0' is the parameter to end the program when order=order0;
# If nargin==2,the algorithm will be terminated untill there is no change
# in network toplogy.
# Output:
# 'G' is the 0-1 network or graph after pc algorithm;
# 'Gval' is the network with strenthness of dependence;
# 'order' is the order of the pc algorithm, here is equal to order0;
# Example:
#
# Author: Xiujun Zhang.
# Version: Sept.2014.
#
# Downloaded from: http://www.comp-sysbio.org/cmi2ni/

CMI2NI <- function(data,lamda,order0=NULL){
  data = t(data) # In original implementation, data is given as pxn
  n_gene <- dim(data)[1]
  G <- matrix(1,n_gene,n_gene)
  G[lower.tri(G, diag=T)] <- 0 
  G <- G+t(G)
  Gval <- G
  order <- -1
  t <- 0
  while (t == 0){
     order <- order+1
     if (!is.null(order0)){
       if (order>order0){
           order <- order-1
           return(list(G=G,Gval=Gval,order=order))
       }
     }
     res.temp <- edgereduce(G,Gval,order,data,t,lamda) 
     G = res.temp$G
     Gval = res.temp$Gval
     t = res.temp$t

     if (t==0){
          print('No edge is reduced! Algorithm  finished!')
          return(list(G=G,Gval=Gval,order=order))
     } else {
          t <- 0
     }
  }
  order <- order-1 # The value of order is the last order of the algorithm
  return(list(G=G,Gval=Gval,order=order))
}

## edgereduce
edgereduce <- function(G,Gval,order,data,t,lamda){
  G0 <- G
  if (order==0){
    for (i in 1:dim(G)[1]){
        for (j in 1:dim(G)[1]){
            if (G[i,j]!=0){
                cmiv <- cmi(data[i,],data[j,])
                Gval[i,j] <- cmiv
                Gval[j,i] <- cmiv
                if (cmiv<lamda){
                    G[i,j] <- 0
                    G[j,i] <- 0
                }
            }
        }
    }
    t <- t+1
    return(list(G=G,Gval=Gval,t=t))
  }
  else {
    for (i in 1:dim(G)[1]){
      for (j in 1:dim(G)[1]){
          if (G[i,j]!=0){
              adj <- c()
              for (k in 1:dim(G)[1]){
                  if (G[i,k]!=0 & G[j,k]!=0){
                      adj <- c(adj,k)
                  }
              }
              if (length(adj)>=order){
                   if(length(adj)==1){ # Need a special case as combn does not work allow single element arguments
                     combntnslist <- as.matrix(adj)
                   }
                   else{
                     combntnslist <- t(combn(adj,order))
                   }
                   combntnsrow <- dim(combntnslist)[1] 
                   
                   cmiv <- 0
                   v1 <- data[i,]
                   v2 <- data[j,]
                   for (k in 1:combntnsrow){
                     vcs <- data[combntnslist[k,],]
                     a <- MI2(v1,v2,vcs) 
                     cmiv <- max(cmiv,a)
                   }
                   Gval[i,j] <- cmiv
                   Gval[j,i] <- cmiv
                   if (cmiv<lamda){
                      G[i,j] <- 0
                      G[j,i] <- 0
                   }
                   t <- t+1
              }
          }
      }
    }
    return(list(G=G,Gval=Gval,t=t))
  }
}

## compute conditional mutual information of x and y
cmi <- function(v1,v2,vcs=NULL){
 if (is.null(vcs)){
      c1 <- (var(v1)) # removed det as this is just a scalar
      c2 <- (var(v2))
      c3 <- det(cov(t(rbind(v1,v2)))) # cov in r only gives a matrix if we give one, so must use rbind
      cmiv <- 0.5*log(c1*c2/c3)
  } 
  else if (!is.null(vcs)){
      c1 <- det(cov(t(rbind(v1,vcs))))
      c2 <- det(cov(t(rbind(v2,vcs))))
      c3 <- (var(t(vcs)))
      c4 <- det(cov(t(rbind(v1,v2,vcs))))
      cmiv <- 0.5*log((c1*c2)/(c3*c4))
  }
  if (cmiv == Inf){
      cmiv <- 1.0e+010
  }
  return(cmiv)
}

# Conditional mutual inclusive information (CMI2)
MI2 <- function(x,y,z){
  r_dmi <- (cas(x,y,z) + cas(y,x,z))/2
  return(r_dmi)
}

# x and y are 1*m dimensional vector; z is n1*m dimensional.
cas <- function(x,y,z){
  # x,y,z are row vectors;
  if(is.null(dim(z))){
    n1 = 1
  }
  else{
    n1 <- dim(z)[1] # This is equal to the order
  }
  n <- n1 +2

  Cov <- var(x)
  Covm <- cov(t(rbind(x,y,z)))
  Covm1 <- cov(t(rbind(x,z)))
  
  InvCov <- solve(Cov)
  InvCovm <- solve(Covm)
  InvCovm1 <- solve(Covm1)

  C11 <- InvCovm1[1,1] # 1x1
  C12 <- 0 # 1x1
  C13 <- InvCovm1[1,2:(1+n1)] # 1x n1
  C23 <- InvCovm[2,3:(2+n1)]-(InvCovm[1,2] * as.numeric(1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1]))) %*% (InvCovm[1,3:(2+n1)] - InvCovm1[1,2:(1+n1)]) # 1x n1
  C22 <- InvCovm[2,2]- InvCovm[1,2]^2 * (1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1])) # 1x1

  C33 <- InvCovm[3:(2+n1),3:(2+n1)]- (as.numeric((1/(InvCovm[1,1]-InvCovm1[1,1]+InvCov[1,1]))) * (t(t(InvCovm[1,3:(2+n1)]-
                                          InvCovm1[1,2:(1+n1)]))%*%(InvCovm[1,3:(2+n1)]-InvCovm1[1,2:(1+n1)]))) # n1 x n1

  InvC <- rbind(c(C11,C12,C13),c(C12,C22,C23),cbind(C13,t(C23),C33))

  C0 <- Cov * (InvCovm[1,1] - InvCovm1[1,1] + InvCov[1,1])
  CS <- 0.5 * (sum(diag(InvC%*%Covm))+log(C0)-n)
  return(CS)
}
