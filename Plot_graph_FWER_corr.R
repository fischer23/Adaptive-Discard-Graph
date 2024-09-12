#Generates plots for Figure 9 of the paper 
#"ADDIS-Graphs for online error control with application to platform trials"

rm(list=ls())
library(ggplot2)
library(MASS)
library(patchwork)
library(Matrix)
library(mvtnorm)

###load the procedures
source("Graph_Procedures.R")

###Gaussian testing problem for exact p-values
m=1000   #Number of Trials
n=100    #Number of Hypotheses per Trial
mu_A=3    #Strength of the alternative
mu_N=0   #Conservativeness of null p-values (<0 for conservative null p-values)
#pi_A is defined in the loop below

###Initialise Hyperparameters
seq_1_n=seq(1,n,1)
gamma=6/(pi^2*(seq_1_n^2))
alpha=0.2
tau=1
lambda=0.2
lambda_corr=0.2
#Correlation within batch
rho=0.5

###Set seed to make the results reproducible
set.seed(12345)

#Set the different batch sizes
batch_sizes=c(1,5,10,20)

###Predefine vectors for FWER and power of the different procedures
FWER_corr=matrix(0,9,length(batch_sizes))
power_corr=matrix(0,9,length(batch_sizes))
FWER_Graph=matrix(0,9,length(batch_sizes))
power_Graph=matrix(0,9,length(batch_sizes))

b=1 #Counter 

for(batch_size in batch_sizes){
  
###Parameters for a local dependence structure given by batches
batch_number=ceiling(n/batch_size)
lags=rep(seq(0,(batch_size-1),1),batch_number)
lags=lags[1:n]
sigma=matrix(rho,batch_size,batch_size)+diag((1-rho),batch_size)
mu=rep(0,batch_size)
corr=as.matrix(bdiag(replicate(batch_number,sigma,simplify=FALSE)))

#Calculate d_j
k_lag=seq(1,n,1)-lags
d_j=rep((n+1),n)
for(k in 1:n){
  for(i in k:n){
    if(k_lag[i]>k){
      d_j[k]=i
      break
    }
  }
}

#Calculate local dependence adjusted weights
w_lda=matrix(0,n,n)
for(k in 1:n){
  if(d_j[k]>n){
    break
  }
  else if(d_j[k]==k+1){
    w_lda[k,d_j[k]:n]=gamma[1:(n-k)]
  }
  else{
    gamma_exh=gamma[(d_j[k]-k):(n-k)]/(1-sum(gamma[1:(d_j[k]-k-1)]))
    w_lda[k,d_j[k]:n]=gamma_exh
  }
}

###Generate p-values and compute FWER and power for the desired procedures
for(l in 1:9){
  
  ##Fast way for p-values in batches.
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  Z=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rep(0,batch_number*batch_size)
    for(k in 1:batch_number){
    X[((k-1)*batch_size+1):(k*batch_size)]=mvrnorm(1,mu,sigma)
    }
    X=X[1:n]
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }

##Adaptive_graph_{corr}
  V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
  power=rep(0,m)  #Power within each trial
  
  for(j in 1:m){
    alpha_ind=ADDIS_Graph_corr(alpha, corr, gamma, w_lda, lambda_corr,lags, p[,j], n)$alpha_ind
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
  }
  FWER_corr[l,b]=mean(V)
  power_corr[l,b]=mean(power, na.rm=TRUE)
  

##ADDIS-Graph_{local-u}
V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
power=rep(0,m)  #Power within each trial

for(j in 1:m){
  alpha_ind=ADDIS_Graph(alpha, gamma,w_lda, tau,lambda,lags, p[,j], n)
  hypo_est=alpha_ind>=p[,j]
  V[j]=max((hypo[,j]==0 & hypo_est==1))
  D=(hypo[,j]==1 & hypo[,j]==hypo_est)
  power[j]=sum(D)/sum(hypo[,j])
}
FWER_Graph[l,b]=mean(V)
power_Graph[l,b]=mean(power, na.rm=TRUE)

}
b=b+1
}
###Create Plot for Adaptive-Graph_{corr}

results_df=data.frame(seq(0.1,0.9,0.1), power_corr,FWER_corr)

p_corr=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2", shape = "2")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2", shape = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X3, linetype = "3",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2", shape = "2")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X4, linetype = "4",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2", shape = "2")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2", shape = "2"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("2"=15))+
  scale_colour_manual(guide="none", values=c( "2"="#ff5c33"))+
  scale_linetype_manual(name  ="Batch-size", values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("1","5","10", "20"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("Adaptive-Graph"[corr])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FWER_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4", shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X3, linetype = "3",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X4, linetype = "4",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4", shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  ="Batch-size", values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("1","5","10", "20"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("Adaptive-Graph"[conf]))

combined <- p_corr + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_ADDIS_Graph_corr_lags.pdf", width = 7.64, height = 3.59)


############Simulations for different rho

#Set different rhos
rhos=c(0.3,0.5,0.7,0.9)

###Set seed to make the results reproducible
set.seed(12345)

#Set the batch sizes
batch_size=10

###Predefine vectors for FWER and power of the different procedures
FWER_corr=matrix(0,9,length(batch_sizes))
power_corr=matrix(0,9,length(batch_sizes))
FWER_Graph=matrix(0,9,length(batch_sizes))
power_Graph=matrix(0,9,length(batch_sizes))
b=1 #Counter 
for(rho in rhos){
  ###Parameters for a local dependence structure given by batches
  batch_number=ceiling(n/batch_size)
  lags=rep(seq(0,(batch_size-1),1),batch_number)
  lags=lags[1:n]
  sigma=matrix(rho,batch_size,batch_size)+diag((1-rho),batch_size)
  mu=rep(0,batch_size)
  corr=as.matrix(bdiag(replicate(batch_number,sigma,simplify=FALSE)))
  
  #Calculate d_j
  k_lag=seq(1,n,1)-lags
  d_j=rep((n+1),n)
  for(k in 1:n){
    for(i in k:n){
      if(k_lag[i]>k){
        d_j[k]=i
        break
      }
    }
  }
  
  #Calculate local dependence adjusted weights
  w_lda=matrix(0,n,n)
  
  for(k in 1:n){
    if(d_j[k]>n){
      break
    }
    else if(d_j[k]==k+1){
      w_lda[k,d_j[k]:n]=gamma[1:(n-k)]
    }
    else{
      gamma_exh=gamma[(d_j[k]-k):(n-k)]/(1-sum(gamma[1:(d_j[k]-k-1)]))
      w_lda[k,d_j[k]:n]=gamma_exh
    }
  }
  
  ###Generate p-values and compute FWER and power for the desired procedures
  for(l in 1:9){
    ##Fast way for p-values in batches.
    pi_A=l/10
    p=matrix(,nrow=n,ncol=m)
    Z=matrix(,nrow=n,ncol=m)
    hypo=matrix(,nrow=n,ncol=m)
    for(j in 1:m){
      hypo[,j]=rbinom(n,1,pi_A)
      X=rep(0,batch_number*batch_size)
      for(k in 1:batch_number){
        X[((k-1)*batch_size+1):(k*batch_size)]=mvrnorm(1,mu,sigma)
      }
      X=X[1:n]
      Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
      p[,j]=pnorm(-Z)
    }
    
    ##Adaptive_graph_{corr}
    V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
    power=rep(0,m)  #Power within each trial
    
    for(j in 1:m){
      alpha_ind=ADDIS_Graph_corr(alpha, corr, gamma, w_lda, lambda_corr,lags, p[,j], n)$alpha_ind
      hypo_est=alpha_ind>=p[,j]
      V[j]=max((hypo[,j]==0 & hypo_est==1))
      D=(hypo[,j]==1 & hypo[,j]==hypo_est)
      power[j]=sum(D)/sum(hypo[,j])
    }
    FWER_corr[l,b]=mean(V)
    power_corr[l,b]=mean(power, na.rm=TRUE)
    
    
    ##ADDIS-Graph_{local-u}
    V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
    power=rep(0,m)  #Power within each trial
    
    for(j in 1:m){
      alpha_ind=ADDIS_Graph(alpha, gamma,w_lda, tau,lambda,lags, p[,j], n)
      hypo_est=alpha_ind>=p[,j]
      V[j]=max((hypo[,j]==0 & hypo_est==1))
      D=(hypo[,j]==1 & hypo[,j]==hypo_est)
      power[j]=sum(D)/sum(hypo[,j])
    }
    FWER_Graph[l,b]=mean(V)
    power_Graph[l,b]=mean(power, na.rm=TRUE)
    
  }
  b=b+1
}
###Create Plot for Adaptive-Graph_{corr}

results_df=data.frame(seq(0.1,0.9,0.1), power_corr,FWER_corr)

p_corr=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2", shape = "2")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2", shape = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X3, linetype = "3",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2", shape = "2")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X4, linetype = "4",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2", shape = "2")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2", shape = "2"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("2"=15))+
  scale_colour_manual(guide="none", values=c( "2"="#ff5c33"))+
  scale_linetype_manual(name  =expression(paste(rho, "  ")), values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("0.3","0.5","0.7", "0.9"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("Adaptive-Graph"[corr])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FWER_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4", shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X3, linetype = "3",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X4, linetype = "4",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4", shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  =expression(paste(rho, "  ")), values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("0.3","0.5","0.7", "0.9"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("Adaptive-Graph"[conf]))

combined <- p_corr + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_ADDIS_Graph_corr_rhos.pdf", width = 7.64, height = 3.59)


############Simulations for different mu_N
mu_Ns=c(0,-0.5,-1,-2)
tau=0.8
lambda=0.16
#Set different rhos
rho=0.5

###Set seed to make the results reproducible
set.seed(12345)

#Set the batch sizes
batch_size=10

###Predefine vectors for FWER and power of the different procedures
FWER_corr=matrix(0,9,length(batch_sizes))
power_corr=matrix(0,9,length(batch_sizes))
FWER_Graph=matrix(0,9,length(batch_sizes))
power_Graph=matrix(0,9,length(batch_sizes))
b=1 #Counter 
for(mu_N in mu_Ns){
  ###Parameters for a local dependence structure given by batches
  batch_number=ceiling(n/batch_size)
  lags=rep(seq(0,(batch_size-1),1),batch_number)
  lags=lags[1:n]
  sigma=matrix(rho,batch_size,batch_size)+diag((1-rho),batch_size)
  mu=rep(0,batch_size)
  corr=as.matrix(bdiag(replicate(batch_number,sigma,simplify=FALSE)))
  
  #Calculate d_j
  k_lag=seq(1,n,1)-lags
  d_j=rep((n+1),n)
  for(k in 1:n){
    for(i in k:n){
      if(k_lag[i]>k){
        d_j[k]=i
        break
      }
    }
  }
  
  #Calculate local dependence adjusted weights
  w_lda=matrix(0,n,n)
  
  for(k in 1:n){
    if(d_j[k]>n){
      break
    }
    else if(d_j[k]==k+1){
      w_lda[k,d_j[k]:n]=gamma[1:(n-k)]
    }
    else{
      gamma_exh=gamma[(d_j[k]-k):(n-k)]/(1-sum(gamma[1:(d_j[k]-k-1)]))
      w_lda[k,d_j[k]:n]=gamma_exh
    }
  }
  
  ###Generate p-values and compute FWER and power for the desired procedures
  for(l in 1:9){
    ##Fast way for p-values in batches.
    pi_A=l/10
    p=matrix(,nrow=n,ncol=m)
    Z=matrix(,nrow=n,ncol=m)
    hypo=matrix(,nrow=n,ncol=m)
    for(j in 1:m){
      hypo[,j]=rbinom(n,1,pi_A)
      X=rep(0,batch_number*batch_size)
      for(k in 1:batch_number){
        X[((k-1)*batch_size+1):(k*batch_size)]=mvrnorm(1,mu,sigma)
      }
      X=X[1:n]
      Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
      p[,j]=pnorm(-Z)
    }
    
    ##Adaptive_graph_{corr}
    V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
    power=rep(0,m)  #Power within each trial
    
    for(j in 1:m){
      alpha_ind=ADDIS_Graph_corr(alpha, corr, gamma, w_lda, lambda_corr,lags, p[,j], n)$alpha_ind
      hypo_est=alpha_ind>=p[,j]
      V[j]=max((hypo[,j]==0 & hypo_est==1))
      D=(hypo[,j]==1 & hypo[,j]==hypo_est)
      power[j]=sum(D)/sum(hypo[,j])
    }
    FWER_corr[l,b]=mean(V)
    power_corr[l,b]=mean(power, na.rm=TRUE)
    
    
    ##ADDIS-Graph_{local-u}
    V=rep(0,m)      #Indicates, whether there was atleast one type 1 error in a trial
    power=rep(0,m)  #Power within each trial
    
    for(j in 1:m){
      alpha_ind=ADDIS_Graph(alpha, gamma, w_lda, tau,lambda,lags, p[,j], n)
      hypo_est=alpha_ind>=p[,j]
      V[j]=max((hypo[,j]==0 & hypo_est==1))
      D=(hypo[,j]==1 & hypo[,j]==hypo_est)
      power[j]=sum(D)/sum(hypo[,j])
    }
    FWER_Graph[l,b]=mean(V)
    power_Graph[l,b]=mean(power, na.rm=TRUE)
    
  }
  b=b+1
}
###Create Plot for Adaptive-Graph_{corr}

results_df=data.frame(seq(0.1,0.9,0.1), power_corr,FWER_corr)

p_corr=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2", shape = "2")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2", shape = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X3, linetype = "3",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2", shape = "2")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2", shape = "2"))+
  geom_line(aes(y = X4, linetype = "4",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2", shape = "2")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2", shape = "2"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("2"=15))+
  scale_colour_manual(guide="none", values=c( "2"="#ff5c33"))+
  scale_linetype_manual(name  =expression(paste(mu[N], "  ")), values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("0","-0.5","-1", "-2"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("Adaptive-Graph"[corr])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FWER_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "1",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4", shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "1",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4", shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4", shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X3, linetype = "3",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4", shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "3",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4", shape = "4"))+
  geom_line(aes(y = X4, linetype = "4",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4", shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "4",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4", shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  =expression(paste(mu[N], "  ")), values = c("1"="solid","2"="longdash","3"="dashed","4"="dotted"),
                        labels=c("0","-0.5","-1", "-2"))+
  xlab(expression(pi[A]))+
  ylab("FWER / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("ADDIS-Graph"[conf]))

combined <- p_corr + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FWER_ADDIS_Graph_corr_mu_Ns.pdf", width = 7.64, height = 3.59)

