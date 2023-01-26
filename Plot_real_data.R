#Generates Figure 10 of the paper "An Adaptive-Discard-Graph for online FWER control"

rm(list=ls())
library(ggplot2)
library(MASS)
library(patchwork)

###load the procedures
source("Procedures.R")

data=read.csv("IMPC_ProcessedData_Continuous.csv")

#Only consider non-missing p-values
data=data[!is.na(data$Sex.P.Val),]

#Number of hypotheses
n=5000

#Sort p-values with respect to their experimental id
mydata=data[order(data$Experimental.Id),]
mydata=mydata[1:n,]
mydata$my_id=seq(1,n)
mydata$lag=rep(0,n)

#Calculate lags for the local dependence structure given by the experimental id
for(i in 1:n){
  mydata_sameexp=mydata[mydata$Experimental.Id==mydata$Experimental.Id[i],]
  first_id=min(mydata_sameexp$my_id)
  mydata$lag[i]=i-first_id
}


###Initialise Hyperparameters
seq_1_n=seq(1,n,1)
gamma=(1/((seq_1_n+1)*log(seq_1_n+1)^2))/2.06227      #2.06227 is the approximated value such that the series equals 1
tau=0.8
lambda=0.16

#Setting of lags and p-values
lags=mydata$lag
p=mydata$Sex.P.Val
  
#Choosing different overall significance levels
alpha_vec=c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4)

b=1 #Counter 

#Weights for ADDIS-Graph_{local} 
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

###Predefine vectors for number of rejections and individual significance levels for the different procedures
n_rej_Spending=matrix(0,length(alpha_vec))    
n_rej_Graph=matrix(0,length(alpha_vec))
alpha_ind_Spending=matrix(0,n,length(alpha_vec)) 
alpha_ind_Graph=matrix(0,n,length(alpha_vec)) 

for(alpha in alpha_vec){
  
##ADDIS-Spending
alpha_ind_Spending[,b]=ADDIS_Spending(alpha, gamma, tau,lambda,lags, p, n)
rej_Spending=alpha_ind_Spending[,b]>=p
n_rej_Spending[b]=sum(rej_Spending)

##ADDIS-Graph_{local}
  alpha_ind_Graph[,b]=ADDIS_Graph(alpha, gamma,w_lda, tau,lambda,lags, p, n)
  rej_Graph=alpha_ind_Graph[,b]>=p
  n_rej_Graph[b]=sum(rej_Graph)

b=b+1

}

###Create Plot for number of rejections vs. FWER level
cols <- c("#f84f4f",  "cornflowerblue")

lab=c(expression("ADDIS-Spending"[local]),expression("ADDIS-Graph"[local]))

results_df=data.frame(seq(0.05,0.4,0.05), n_rej_Spending,n_rej_Graph)

plot_real=ggplot(results_df, aes(seq.0.05..0.4..0.05.)) + 
  geom_line(aes(y = n_rej_Spending, colour = "ADDIS-Spending")) + 
  geom_point(aes(y = n_rej_Spending, colour = "ADDIS-Spending", shape = "ADDIS-Spending")) +
  geom_line(aes(y = n_rej_Graph, colour = "ADDIS-Graph")) + 
  geom_point(aes(y = n_rej_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph")) +
  scale_shape_manual(name  ="Procedure",values=c("ADDIS-Spending"=16,"ADDIS-Graph"=17),labels = c(lab[1],lab[2]))+
  scale_color_manual(name  ="Procedure", labels = c(lab[1],lab[2]), values = c("ADDIS-Spending"=cols[1],"ADDIS-Graph"=cols[2]))+
  xlab("FWER level")+
  ylab("Number of rejections")+
  scale_x_continuous(breaks = seq(0.1,0.35,0.05), limits=c(0.05,0.4),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(150,300,30), limits=c(150,300),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),legend.text.align = 0)

###Create Plot for individual significance levels


results_df_2=data.frame(index=seq(100,n,1), ind_Spending=alpha_ind_Spending[100:n,4],ind_Graph=alpha_ind_Graph[100:n,4])

plot_real_2=ggplot(results_df_2, aes(index)) +
  geom_line(aes(y =  rep(0,n-100+1), colour = "ADDIS-Spending"))+
  geom_point(aes(y = ind_Spending, colour = "ADDIS-Spending", shape = "ADDIS-Spending", size="ADDIS-Spending")) +
  geom_line(aes(y =  rep(0,n-100+1), colour = "ADDIS-Graph"))+
  geom_point(aes(y = ind_Graph, colour = "ADDIS-Graph", shape = "ADDIS-Graph", size="ADDIS-Graph")) +
  scale_size_manual(guide="none", name  ="Procedure",values=c("ADDIS-Spending"=0.6,"ADDIS-Graph"=0.6),labels = c(lab[1],lab[2]))+
  scale_shape_manual(name  ="Procedure",values=c("ADDIS-Spending"=16,"ADDIS-Graph"=17),labels = c(lab[1],lab[2]))+
  scale_color_manual(name  ="Procedure", labels = c(lab[1],lab[2]), values = c("ADDIS-Spending"=cols[1],"ADDIS-Graph"=cols[2]))+
  xlab("Index of hypothesis")+
  ylab("Individual significance level")+
  scale_x_continuous(breaks = seq(0,5000,1000), limits=c(0,5000),expand = c(0, 0))+
  scale_y_continuous(breaks = seq(0,0.0003,0.00005), limits=c(0,0.0003),expand = c(0, 0))+
  theme(panel.grid.major = element_line(color = "grey", size = 0.5,linetype = 1),
        panel.border = element_rect(colour = "black", fill=NA, size=1),legend.text.align = 0)

combined <- plot_real + plot_real_2  & theme(legend.position = "bottom")
combined=combined + plot_layout(guides = "collect")

ggsave("Plot_real_data_combined.pdf",plot=combined)