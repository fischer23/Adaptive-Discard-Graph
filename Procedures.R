###This R-file contains the procedures that were proposed in the paper "An Adaptive-Discard-Graph for online error control"

#Input (at most): alpha, gamma, w, tau, lambda,lags, p, n. 

#n:     Number of hypotheses.
#alpha: Overall significance level. Number between 0 and 1.
#gamma: Weights for Alpha-Spending. Non-negative n-dimensional vector with sum less than 1.
#w:     Weights for graphical procedures. Upper triangle n x n matrix with row sum less than 1.
#h:     Additional weights for FDR-ADDIS-Graph. Upper triangle n x n matrix with row sum less than 1.
#tau:   Used for discarding procedures. n-dimensional vector with values between 0 and 1. Besides, it can be chosen as a fixed number between 0 and 1.
#lambda:Used for the adaptive procedures. n-dimensional vector with values between 0 and tau_i. Besides, it can be chosen as a fixed number fulfilling these conditions.
#lags:  Representing the given local dependence structure. n-dimensional vector of natural numbers (including 0) or a fix natural number (including 0).
#e:     Stopping times in an asynchronous testing setup. n-dimensional vector of natural numbers (e_1,...,e_n) such that e_i >= i for all i=1,...,n.
#p:     n-dimensional vector of p-values. 


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

 #ADDIS-Graph_{local}
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
  alpha_ind[1]=min(alpha*gamma[1]*(tau[1]-lambda[1]), lambda[1])
  if(p[1]>tau[1] | p[1]<=lambda[1]){
    C_S[1]=1
  }
  if(p[1]<=alpha_ind[1]){
    count=1
  }
  for(i in 2:n){
    e_ind=e[1:(i-1)]<=rep(i,i-1)  
    alpha_ind[i]=min(lambda[i],(alpha*gamma[i]+sum(C_S[1:(i-1)]*e_ind[1:(i-1)]*alpha_ind[1:(i-1)]*w[1:(i-1),i]/(tau[1:(i-1)]-lambda[1:(i-1)]))+sum(e_ind[1:(i-1)]*h[(1:i-1),i]*R[1:(i-1)])*alpha)*(tau[i]-lambda[i])) 
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



