###This R-file contains the procedures that were proposed in the paper 
#"ADDIS-Graphs for online error control with application to platform trials"

#Input (at most): alpha, gamma, w,h, tau, lambda,lags,d_j,e, p, n. 

#n:     Number of hypotheses.
#alpha: Overall significance level. Number between 0 and 1.
#gamma: Weights for Alpha-Spending. Non-negative n-dimensional vector with sum less than 1.
#w:     Weights for graphical procedures. Upper triangle n x n matrix with row sum less than 1.
#h:     Additional weights for FDR-ADDIS-Graph. Upper triangle n x n matrix with row sum less than 1.
#tau:   Used for discarding procedures. n-dimensional vector with values between 0 and 1. Besides, it can be chosen as a fixed number between 0 and 1.
#lambda:Used for the adaptive procedures. n-dimensional vector with values between 0 and tau_i. Besides, it can be chosen as a fixed number fulfilling these conditions.
#lags:  Representing the given local dependence structure. n-dimensional vector of natural numbers (including 0) or a fix natural number (including 0).
#d_j:   First future $p$-value that is independent of P_j. n-dimensional vector of natural numbers with d_j>j.
#e:     Stopping times in an asynchronous testing setup. n-dimensional vector of natural numbers (e_1,...,e_n) such that e_i >= i for all i=1,...,n.
#p:     n-dimensional vector of p-values. 
#corr:  Correlation matrix of the normally distributed test statistics. Symmetric and positiv-semidefinite nxn matrix.


#ADDIS-Spending_{local} 
ADDIS_Spending=function(alpha, gamma, tau, lambda,lags, p, n){
  if(length(lags)== 1){lags=rep(lags,n)}
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
    warning("mismatching length")
  }
  S_C=rep(0,n)
  t=rep(1,n)
  alpha_ind=rep(0,n)
  for(i in 1:n){
    if(lags[i]>=i-1){
      t[i]=i
    }else{
      t[i]=1+lags[i]+sum(S_C[1:(i-lags[i]-1)])
    }
    alpha_ind[i]=alpha*gamma[t[i]]*(tau[i]-lambda[i])
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      S_C[i]=1
    }
  }
  return(alpha_ind)
}

#Closed ADDIS-Spending_{local}
closed_ADDIS_Spending=function(alpha, gamma, tau, lambda,lags, p, n){
  if(length(lags)== 1){lags=rep(lags,n)}
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
    warning("mismatching length")
  }
  S_C=rep(0,n)
  R=rep(0,n)
  t=rep(1,n)
  alpha_ind=rep(0,n)
  for(i in 1:n){
    if(lags[i]>=1){
      if(lags[i]>=i-1){
        t[i]=1+sum(1-(R[(i-lags[i]):(i-1)]))
      }else{
        t[i]=1+sum(S_C[1:(i-lags[i]-1)])+sum(1-(R[(i-lags[i]):(i-1)]))
      }
    }else{
      t[i]=1+sum(S_C[1:(i-1)])
    }
    alpha_ind[i]=alpha*gamma[t[i]]*(tau[i]-lambda[i])
    if(p[i]<=tau[i] & p[i]>lambda[i]){
      S_C[i]=1
    }
    if(p[i]<=alpha_ind[i]){
      R[i]=1
    }
  }
  return(alpha_ind)
}

 #ADDIS-Graph_{conf}
 ADDIS_Graph=function(alpha, gamma, w, tau, lambda,lags, p, n){
   if(length(lags)== 1){lags=rep(lags,n)}
   if(length(tau)== 1){tau=rep(tau,n)}
   if(length(lambda)== 1){lambda=rep(lambda,n)}
   if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
     warning("mismatching length")
   }
   C_S=rep(0,n)
   alpha_ind=rep(0,n)
   alpha_ind[1]=alpha*gamma[1]*(tau[1]-lambda[1])
   if(p[1]>tau[1] | p[1]<=lambda[1]){
     C_S[1]=1
   }
   for(i in 2:n){
     if(lags[i]>=i-1){
       alpha_ind[i]=alpha*gamma[i]*(tau[i]-lambda[i])
     }else{
       alpha_ind[i]=(alpha*gamma[i]+sum(C_S[1:(i-lags[i]-1)]*alpha_ind[1:(i-lags[i]-1)]*w[1:(i-lags[i]-1),i]/(tau[1:(i-lags[i]-1)]-lambda[1:(i-lags[i]-1)])))*(tau[i]-lambda[i]) 
     }
     if(p[i]>tau[i] | p[i]<=lambda[i]){
       C_S[i]=1
     }
   }
   return(alpha_ind)
 }
 
 #ADDIS-Graph_{conf-u}
 ADDIS_Graph_imp=function(alpha, gamma, tau, lambda,lags, d_j, p, n){
   if(length(lags)== 1){lags=rep(lags,n)}
   if(length(tau)== 1){tau=rep(tau,n)}
   if(length(lambda)== 1){lambda=rep(lambda,n)}
   if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
     warning("mismatching length")
   }
   S_C=(p>lambda & p<=tau)
   C_S=1-S_C
   w=matrix(0,n,n)
   w[1,2:n]=(gamma[1:(n-1)]-gamma[2:n])/gamma[1]
   for(j in 2:(n-1)){
     w[j,(j+1):n]=(gamma[(sum(S_C[1:(j-1)])+1):(sum(S_C[1:(j-1)])+n-j)]-gamma[(sum(S_C[1:(j-1)])+2):(sum(S_C[1:(j-1)])+n-j+1)])/gamma[sum(S_C[1:(j-1)])+1]
   }
   w_star=w
   w_minus=matrix(0,n,n)
   for(j in 1:(n-2)){
     for(i in (j+1):(n-1)){
     if(i<d_j[j]){
      w_minus[j,i]=w_star[j,i]
      w_star[j,i]=0
      w_star[j,(i+1):n]= w_star[j,(i+1):n]+w_minus[j,i]*w[i,(i+1):n]
     }else{
      w_minus[j,i]=ifelse(lags[i]==0,0,sum(w[(i-lags[i]):(i-1),i]*w_minus[j,(i-lags[i]):(i-1)]))
      w_star[j,i]=w_star[j,i]-w_minus[j,i]
      w_star[j,(i+1):n]= w_star[j,(i+1):n]+w_minus[j,i]*w[i,(i+1):n]
     }  
     }
   }
   for(j in 1:(n-1)){
     if(n<d_j[j]){
       w_star[j,n]=0
     }else{
       w_minus[j,n]=ifelse(lags[n]==0,0,sum(w[(n-lags[n]):(n-1),n]*w_minus[j,(n-lags[n]):(n-1)]))
       w_star[j,n]=w_star[j,n]-w_minus[j,n]
     } 
   }
     
   alpha_ind=rep(0,n)
   alpha_ind[1]=alpha*gamma[1]*(tau[1]-lambda[1])
   for(i in 2:n){
     if(lags[i]>=i-1){
       alpha_ind[i]=alpha*gamma[i]*(tau[i]-lambda[i])
     }else{
       alpha_ind[i]=(alpha*gamma[i]+sum(C_S[1:(i-lags[i]-1)]*alpha_ind[1:(i-lags[i]-1)]*w_star[1:(i-lags[i]-1),i]/(tau[1:(i-lags[i]-1)]-lambda[1:(i-lags[i]-1)])))*(tau[i]-lambda[i]) 
     }
   }
   return(alpha_ind)
 }
 
 
 #ADDIS-Graph_{corr}
 ADDIS_Graph_corr=function(alpha, corr, gamma, w, lambda,lags, p, n){
   if(length(lags)== 1){lags=rep(lags,n)}
   if(length(lambda)== 1){lambda=rep(lambda,n)}
   if(length(gamma)!=n | length(lambda)!=n| length(p)!=n| length(lags)!=n){
     warning("mismatching length")
   }
   Z=-qnorm(p)
   C=(p<=lambda)
   C_idx=which(C==0)
   alpha_ind=rep(0,n)
   alpha_used=rep(0,n)
   
   alpha_ind[1]=alpha*gamma[1]*(1-lambda[1])
   alpha_used[1]=alpha_ind[1]
   
   for(i in 2:n){
     if(lags[i]>=i-1){
       alpha_ind[i]=alpha*gamma[i]*(1-lambda[i])
     }else{
       alpha_ind[i]=(alpha*gamma[i]+sum(C[1:(i-lags[i]-1)]*alpha_ind[1:(i-lags[i]-1)]*w[1:(i-lags[i]-1),i]/(1-lambda[1:(i-lags[i]-1)]))+sum((1-C[1:(i-lags[i]-1)])*(alpha_ind[1:(i-lags[i]-1)]-alpha_used[1:(i-lags[i]-1)])*w[1:(i-lags[i]-1),i]/(1-lambda[1:(i-lags[i]-1)])))*(1-lambda[i]) 
     }
     C_lag=c(C_idx[C_idx>=(i-lags[i]) & C_idx<i],i)
     if(length(C_lag)==1){
       alpha_used[i]=alpha_ind[i]
     }else{
       R_hat=corr[C_lag,C_lag]
       lb=c(rep(-Inf,(length(C_lag)-1)),qnorm(1-alpha_ind[i]))
       ub=c(qnorm(1-alpha_ind[C_lag[1:(length(C_lag)-1)]]),Inf)
       calc=pmvnorm(upper=ub,lower=lb,mean=rep(0,(length(C_lag))),sigma=R_hat)
       alpha_used[i]=calc[1]
     }
   }
   return(data.frame(alpha_ind,alpha_used))
 }

#FDR-ADDIS-Graph_{async}
ADDIS_Graph_fdr=function(alpha, gamma, w,h, tau, lambda, e, p, n){
  if(length(tau)== 1){tau=rep(tau,n)}
  if(length(lambda)== 1){lambda=rep(lambda,n)}
  if(length(gamma)!=n | length(tau)!= n | length(lambda)!=n| length(p)!=n| length(e)!= n){
    warning("mismatching length")
  }
  C_S=rep(0,n)
  R=rep(0,n)
  count=0
  alpha_ind=rep(0,n)
  alpha_graph=rep(0,n)
  alpha_graph=alpha*gamma[1]*(tau[1]-lambda[1])
  alpha_ind[1]=min(alpha_graph[1], lambda[1])
  if(p[1]>tau[1] | p[1]<=lambda[1]){
    C_S[1]=1
  }
  if(p[1]<=alpha_ind[1]){
    count=1
  }
  for(i in 2:n){
    e_ind=e[1:(i-1)]<=rep(i,i-1)
    alpha_graph[i]=(alpha*gamma[i]+sum(C_S[1:(i-1)]*e_ind[1:(i-1)]*alpha_graph[1:(i-1)]*w[1:(i-1),i]/(tau[1:(i-1)]-lambda[1:(i-1)]))+sum(e_ind[1:(i-1)]*h[(1:i-1),i]*R[1:(i-1)])*alpha)*(tau[i]-lambda[i])
    alpha_ind[i]=min(lambda[i],alpha_graph[i]) 
    if(p[i]>tau[i] | p[i]<=lambda[i]){
      C_S[i]=1
    }
    if(p[i]<=alpha_ind[i]){
      R[i]=1
      count=count+1
    }
    if(count==1){
      R[i]=0
    }
  }
  return(alpha_ind)
}



