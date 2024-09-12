#Generates plots for Figure 10 of the paper 
#"ADDIS-Graphs for online error control with application to platform trials"

rm(list=ls())
library(ggplot2)
library(onlineFDR)  #Needed for ADDIS*_{async}
library(MASS)
library(patchwork)

###load the procedures
source("Procedures.R")

###Gaussian testing problem for exact p-values
m=200    #Number of Trials
n=1000    #Number of Hypotheses per Trial
mu_A=3    #Strength of the alternative
mu_N=-0.5    #Conservativeness of null p-values (<0 for conservative null p-values)
#pi_A is defined in the loop below

###Initialise Hyperparameters
seq_1_n=seq(1,n,1)
gamma=(1/((seq_1_n+1)*log(seq_1_n+1)^2))/2.06227      #2.06227 is the approximated value such that the series equals 1
alpha=0.05
tau=0.5
lambda=0.25

#Set stopping times
e_values=c(0,2,5,10)


###Predefine vectors for FWER and power of the different procedures
FDR_ADDIS=matrix(0,9,length(e_values))
power_ADDIS=matrix(0,9,length(e_values))
FDR_Graph=matrix(0,9,length(e_values))
power_Graph=matrix(0,9,length(e_values))

###Set seed to make the results reproducible
set.seed(12345)

b=1 #Counter 
for(e_value in e_values){
  
e=seq(1,n)+rep(e_value,n)

#Adjusted weigths
w_adj=matrix(0,n,n)

for(k in 1:n){
  if(e[k]>=n){
    break
  }
  else if(e[k]==k){
    w_adj[k,(k+1):n]=gamma[1:(n-k)]
  }
  else{
    gamma_exh=gamma[(e[k]-k+1):(n-k)]/(1-sum(gamma[1:(e[k]-k)]))
    w_adj[k,(e[k]+1):n]=gamma_exh
  }
}


###Generate p-values and compute FDR and power for the desired procedures
for(l in 1:9){
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rnorm(n)
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }

##ADDIS*
power=rep(0,m)  #Power within each trial
fdr=rep(0,m)    #FDR within each trial

for(j in 1:m){
  df=data.frame(id=seq_1_n,pval=p[,j],decision.times=e)
  alpha_ind=ADDIS(df,alpha, gamma, tau*lambda*alpha,lambda, tau, TRUE)
  hypo_est=alpha_ind$alphai>=p[,j]
  D=(hypo[,j]==1 & hypo[,j]==hypo_est)
  power[j]=sum(D)/sum(hypo[,j])
  fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
}
FDR_ADDIS[l,b]=mean(fdr)
power_ADDIS[l,b]=mean(power)


##ADDIS-Graph FDR
power=rep(0,m)  #Power within each trial
fdr=rep(0,m)    #FDR within each trial

for(j in 1:m){
  alpha_ind=ADDIS_Graph_fdr(alpha, gamma, w_adj,w_adj, tau, lambda,e, p[,j], n)
  hypo_est=alpha_ind>=p[,j]
  D=(hypo[,j]==1 & hypo[,j]==hypo_est)
  power[j]=sum(D)/sum(hypo[,j])
  fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
}
FDR_Graph[l,b]=mean(fdr)
power_Graph[l,b]=mean(power)
}
b=b+1
}
###Create Plot for exact p-values
Sys.setlocale('LC_CTYPE', 'greek')

results_df=data.frame(seq(0.1,0.9,0.1), power_ADDIS,FDR_ADDIS)

p_ADDIS=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2"))+
  geom_line(aes(y = X3, linetype = "5",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2"))+
  geom_line(aes(y = X4, linetype = "10",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2"))+
  geom_hline(yintercept=alpha)+
  scale_colour_manual(guide="none", values=c( "2"="#f84f4f"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("ADDIS*"[async])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FDR_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4",shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4",shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4",shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X3, linetype = "5",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4",shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X4, linetype = "10",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4",shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4",shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("FDR-ADDIS-Graph"[async]))

combined <- p_ADDIS + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FDR_ADDIS_Graph_logq.pdf")

#############################Results for gamma_i=6/(pi^2*i^2)
gamma=6/(pi^2*(seq_1_n^2))

###Predefine vectors for FWER and power of the different procedures
FDR_ADDIS=matrix(0,9,length(e_values))
power_ADDIS=matrix(0,9,length(e_values))
FDR_Graph=matrix(0,9,length(e_values))
power_Graph=matrix(0,9,length(e_values))

b=1 #Counter 
for(e_value in e_values){
  
e=seq(1,n)+rep(e_value,n)

#Adjusted weigths
w_adj=matrix(0,n,n)
  
for(k in 1:n){
  if(e[k]>=n){
    break
  }
  else if(e[k]==k){
    w_adj[k,(k+1):n]=gamma[1:(n-k)]
  }
  else{
    gamma_exh=gamma[(e[k]-k+1):(n-k)]/(1-sum(gamma[1:(e[k]-k)]))
    w_adj[k,(e[k]+1):n]=gamma_exh
  }
}
  
  
###Generate p-values and compute FDR and power for the desired procedures
for(l in 1:9){
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rnorm(n)
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }
    
  ##ADDIS*
  power=rep(0,m)  #Power within each trial
  fdr=rep(0,m)    #FDR within each trial
  
  for(j in 1:m){
    df=data.frame(id=seq_1_n,pval=p[,j],decision.times=e)
    alpha_ind=ADDIS(df,alpha, gamma, tau*lambda*alpha,lambda, tau, TRUE)
    hypo_est=alpha_ind$alphai>=p[,j]
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
    fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
  }
  FDR_ADDIS[l,b]=mean(fdr)
  power_ADDIS[l,b]=mean(power)
    
    
  ##ADDIS-Graph FDR
  power=rep(0,m)  #Power within each trial
  fdr=rep(0,m)    #FDR within each trial
  
  for(j in 1:m){
    alpha_ind=ADDIS_Graph_fdr(alpha, gamma, w_adj,w_adj, tau, lambda,e, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
    fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
  }
  FDR_Graph[l,b]=mean(fdr)
  power_Graph[l,b]=mean(power)
}
b=b+1
}
###Create Plot for exact p-values
Sys.setlocale('LC_CTYPE', 'greek')

results_df=data.frame(seq(0.1,0.9,0.1), power_ADDIS,FDR_ADDIS)

p_ADDIS=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2"))+
  geom_line(aes(y = X3, linetype = "5",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2"))+
  geom_line(aes(y = X4, linetype = "10",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2"))+
  geom_hline(yintercept=alpha)+
  scale_colour_manual(guide="none", values=c( "2"="#f84f4f"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("ADDIS*"[async])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FDR_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4",shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4",shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4",shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X3, linetype = "5",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4",shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X4, linetype = "10",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4",shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4",shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("FDR-ADDIS-Graph"[async]))

combined <- p_ADDIS + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FDR_ADDIS_Graph_q^2.pdf")

#############################Results for gamma_i\propto 1/(pi^1.6)
gamma=1/(2.28577*seq_1_n^1.6)                      #2.28577 approximated value such that sum of gamma equals one


###Predefine vectors for FWER and power of the different procedures
FDR_ADDIS=matrix(0,9,length(e_values))
power_ADDIS=matrix(0,9,length(e_values))
FDR_Graph=matrix(0,9,length(e_values))
power_Graph=matrix(0,9,length(e_values))

b=1 #Counter 
for(e_value in e_values){
  
e=seq(1,n)+rep(e_value,n)

#Adjusted weigths
w_adj=matrix(0,n,n)
  
for(k in 1:n){
  if(e[k]>=n){
    break
  }
  else if(e[k]==k){
    w_adj[k,(k+1):n]=gamma[1:(n-k)]
  }
  else{
    gamma_exh=gamma[(e[k]-k+1):(n-k)]/(1-sum(gamma[1:(e[k]-k)]))
    w_adj[k,(e[k]+1):n]=gamma_exh
  }
}
  
  
###Generate p-values and compute FDR and power for the desired procedures
for(l in 1:9){
  pi_A=l/10
  p=matrix(,nrow=n,ncol=m)
  hypo=matrix(,nrow=n,ncol=m)
  for(j in 1:m){
    hypo[,j]=rbinom(n,1,pi_A)
    X=rnorm(n)
    Z=mu_N*(hypo[,j]-1)*(-1)+mu_A*hypo[,j]+X
    p[,j]=pnorm(-Z)
  }
    
  ##ADDIS*
  power=rep(0,m)  #Power within each trial
  fdr=rep(0,m)    #FDR within each trial
  
  for(j in 1:m){
    df=data.frame(id=seq_1_n,pval=p[,j],decision.times=e)
    alpha_ind=ADDIS(df,alpha, gamma, tau*lambda*alpha,lambda, tau, TRUE)
    hypo_est=alpha_ind$alphai>=p[,j]
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
    fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
  }
  FDR_ADDIS[l,b]=mean(fdr)
  power_ADDIS[l,b]=mean(power)
    
    
  ##ADDIS-Graph FDR
  power=rep(0,m)  #Power within each trial
  fdr=rep(0,m)    #FDR within each trial
  
  for(j in 1:m){
    alpha_ind=ADDIS_Graph_fdr(alpha, gamma, w_adj,w_adj, tau, lambda,e, p[,j], n)
    hypo_est=alpha_ind>=p[,j]
    V[j]=max((hypo[,j]==0 & hypo_est==1))
    D=(hypo[,j]==1 & hypo[,j]==hypo_est)
    power[j]=sum(D)/sum(hypo[,j])
    fdr[j]=sum((hypo[,j]==0 & hypo_est==1))/max(1,sum(hypo_est))
  }
  FDR_Graph[l,b]=mean(fdr)
  power_Graph[l,b]=mean(power)
}
b=b+1
}
###Create Plot for exact p-values
Sys.setlocale('LC_CTYPE', 'greek')

results_df=data.frame(seq(0.1,0.9,0.1), power_ADDIS,FDR_ADDIS)

p_ADDIS=ggplot(results_df, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "2")) + 
  geom_point(aes(y = X1,colour = "2")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "2"))+
  geom_point(aes(y = X1.1,colour = "2"))+
  geom_line(aes(y = X2, linetype = "2",colour = "2")) + 
  geom_point(aes(y = X2,colour = "2")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "2"))+
  geom_point(aes(y = X2.1,colour = "2"))+
  geom_line(aes(y = X3, linetype = "5",colour = "2")) + 
  geom_point(aes(y = X3,colour = "2")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "2"))+
  geom_point(aes(y = X3.1,colour = "2"))+
  geom_line(aes(y = X4, linetype = "10",colour = "2")) + 
  geom_point(aes(y = X4,colour = "2")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "2"))+
  geom_point(aes(y = X4.1,colour = "2"))+
  geom_hline(yintercept=alpha)+
  scale_colour_manual(guide="none", values=c( "2"="#f84f4f"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("ADDIS*"[async])) 


###Create Plot for ADDIS-Graph
results_df_graph=data.frame(seq(0.1,0.9,0.1), power_Graph,FDR_Graph)

p_graph=ggplot(results_df_graph, aes(seq.0.1..0.9..0.1.)) + 
  geom_line(aes(y = X1, linetype = "0",colour = "4")) + 
  geom_point(aes(y = X1,colour = "4",shape = "4")) +
  geom_line(aes(y = X1.1, linetype = "0",colour = "4"))+
  geom_point(aes(y = X1.1),colour = "4",shape = "4")+
  geom_line(aes(y = X2, linetype = "2",colour = "4")) + 
  geom_point(aes(y = X2,colour = "4",shape = "4")) +
  geom_line(aes(y = X2.1, linetype = "2",colour = "4"))+
  geom_point(aes(y = X2.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X3, linetype = "5",colour = "4")) + 
  geom_point(aes(y = X3,colour = "4",shape = "4")) +
  geom_line(aes(y = X3.1, linetype = "5",colour = "4"))+
  geom_point(aes(y = X3.1,colour = "4",shape = "4"))+
  geom_line(aes(y = X4, linetype = "10",colour = "4")) + 
  geom_point(aes(y = X4,colour = "4",shape = "4")) +
  geom_line(aes(y = X4.1, linetype = "10",colour = "4"))+
  geom_point(aes(y = X4.1,colour = "4",shape = "4"))+
  geom_hline(yintercept=alpha)+
  scale_shape_manual(guide="none", values=c("4"=17))+
  scale_colour_manual(guide="none", values=c( "4"="cornflowerblue"))+
  scale_linetype_manual(name  ="Test duration", values = c("0"="solid","2"="longdash","5"="dashed","10"="dotted"))+
  xlab(expression(pi[A]))+
  ylab("FDR / Power")+
  scale_x_continuous(breaks = seq(0.2,0.8,0.2), limits=c(0.1,0.9),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,1.2,0.2), limits=c(0,1),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1))+
  ggtitle(expression("FDR-ADDIS-Graph"[async]))

combined <- p_ADDIS + p_graph  & theme(legend.position = "bottom")
combined + plot_layout(guides = "collect")
ggsave("Plot_FDR_ADDIS_Graph_q^1,6.pdf")
