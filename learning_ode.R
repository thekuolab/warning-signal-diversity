library(deSolve)
library(rootSolve)
library(tidyverse)
library(plotly)
require(parallel)
library(colorspace)

#### Separate K
# functions
learning_ode_separateK<-function(t,state,params){
  
  r1<-params["r1"]
  r2<-params["r2"]
  K1<-params["K1"]
  K2<-params["K2"]
  lambda1<-params["lambda1"]
  lambda2<-params["lambda2"]
  alphaL1<-0.5+abs(lambda1-0.5)
  alphaL2<-0.5+abs(lambda2-0.5)
  alphaF1<-params["alphaF1"]
  alphaF2<-params["alphaF2"]
  H<-params["H"]
  GL<-params["GL"]
  GF<-params["GF"]
  # fL<-params["fL"]
  # Eff_social_L<-params["Eff_social_L"]
  # Eff_social_F<-params["Eff_social_F"]
  
  x1<-state["x1"]
  x2<-state["x2"]
  P1<-state["P1"]
  P2<-state["P2"]
  
  e1<-x1/(x1+x2)
  e2<-x2/(x1+x2)
  
  # alphaL1_hat<-H*fL*alphaL1+H*(1-fL)*alphaL1*Eff_social_L
  # alphaL2_hat<-H*fL*alphaL2+H*(1-fL)*alphaL2*Eff_social_L
  # 
  # alphaF1_hat<-H*fL*alphaF1+H*(1-fL)*alphaF1*Eff_social_F
  # alphaF2_hat<-H*fL*alphaF2+H*(1-fL)*alphaF2*Eff_social_F
  
  dP1<-alphaL1*(e1*P1+GL*e2*P2)*(lambda1-P1)+alphaF1*(e1*(1-P1)+e2)*(0.5-P1)
  dP2<-alphaL2*(e2*P2+GL*e1*P1)*(lambda2-P2)+alphaF2*(e2*(1-P2)+e1)*(0.5-P2)
  
  dx1<-r1*x1*(1-x1/K1)-H*x1*P1
  dx2<-r2*x2*(1-x2/K2)-H*x2*P2
  
  list(c(dx1,dx2,dP1,dP2))
}

mimic_ode_separateK<-function(x2,r,K2,lambda,alphaF,G,H){
  
  state<-c(x1 = 200, x2 = x2,P1 = 0.5,P2 = 0.5)
  params2<-c(r1=r,r2=r,K1 = 200,K2 = K2,lambda1=lambda,lambda2=lambda,alphaF1=alphaF,alphaF2=alphaF,H=H,GL=G,GF=1)
  time <- seq(0, 10000, by=0.1)
  test<-deSolve::ode(y = state,func = learning_ode_separateK,parms = params2,times = time)
  end<-data.frame(tail(test,n = 10000))
  # popsize_end<-c(mean(end$x1),mean(end$x2))
  # cv<-c(sd(end$x1)/popsize_end[1],sd(end$x2)/popsize_end[2])*100
  
  write_csv(end,paste0("D:/Desktop/mimicODE/separateK/rare",x2,"_r",r,"_K2",K2,"_lambda",lambda,"_alphaF",alphaF,"_G",G,"_H",H,".csv"))
}


# solving ODE
x2<-c(20,100,180)
K2<-c(20,100,180,200)
vector_P<-seq(0,0.5,length.out = 20)
vector_F<-seq(0,1,length.out = 20)
baby<-c(0.2,0.3)
H<-1.2
G<-0.5

all_combo<-matrix(unlist(expand.grid(180,180,vector_P,vector_F,baby)),ncol=5)
parameter_space<-all_combo[all_combo[,2]>=all_combo[,1],]

clust<-makeCluster(no.Cores, setup_timeout = 0.5)

clusterEvalQ(clust,library(tidyverse))
clusterExport(clust,"learning_ode_separateK")
clusterExport(clust,"mimic_ode_separateK")
clusterExport(clust,"parameter_space")
clusterExport(clust,"H")
clusterExport(clust,"G")

parApply(cl = clust,parameter_space,1,function(y){mimic_ode_separateK(x2 = y[1],r = y[5],K2 = y[2],lambda = y[3],alphaF = y[4],G = G,H = H)})

stopCluster(clust)

# run the parameter combos that were left out
files<-list.files("D:/Desktop/mimicODE/separateK")

combos<-NULL

for (i in 1:nrow(parameter_space)){
  
  combos[i]<-(paste0("rare",parameter_space[i,1],"_r",parameter_space[i,5],"_K2",parameter_space[i,2],"_lambda",parameter_space[i,3],"_alphaF",parameter_space[i,4],"_G",G,"_H",H,".csv"))
  
  
}

missing<-setdiff(combos,files)

parameter_space_missing<-matrix(999,nrow = length(missing),ncol = 5)

for (i in 1:length(missing)){
  
  parameter_space_missing[i,1]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
  parameter_space_missing[i,2]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][3],pattern = "K2")[[1]][2])
  parameter_space_missing[i,3]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][4],pattern = "lambda")[[1]][2]) 
  parameter_space_missing[i,4]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][5],pattern = "alphaF")[[1]][2])
  parameter_space_missing[i,5]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][2],pattern = "r")[[1]][2])
}

clust<-makeCluster(no.Cores, setup_timeout = 0.5)

clusterEvalQ(clust,library(tidyverse))
clusterExport(clust,"learning_ode_separateK")
clusterExport(clust,"mimic_ode_separateK")
clusterExport(clust,"parameter_space_missing")
clusterExport(clust,"H")
clusterExport(clust,"G")

parApply(cl = clust,parameter_space_missing,1,function(y){mimic_ode_separateK(x2 = y[1],r = y[5],K2 = y[2],lambda = y[3],alphaF = y[4],G = G,H = H)})

stopCluster(clust)

# summarizing results
G<-0.5
H<-1.2

files<-list.files(paste0("D:/Desktop/mimicODE/separateK/G",G,"/H",H))
df<-matrix(999,nrow = length(files),ncol = 13)


for (i in 1:length(files)){
  
  dff<-read_csv(paste0("D:/Desktop/mimicODE/separateK/G",G,"/H",H,"/",files[i]))
  df[i,1]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][2],pattern = "r")[[1]][2])
  df[i,2]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
  df[i,3]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][3],pattern = "K2")[[1]][2])
  df[i,4]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][4],pattern = "lambda")[[1]][2])
  df[i,5]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][5],pattern = "alphaF")[[1]][2])
  df[i,6]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][6],pattern = "G")[[1]][2])
  df[i,7]<-H
  df[i,8]<-mean(dff$x1)
  df[i,9]<-mean(dff$x2)
  df[i,10]<-sd(dff$x1)/mean(dff$x1)
  df[i,11]<-sd(dff$x2)/mean(dff$x2)
  df[i,12]<-mean(dff$P1)
  df[i,13]<-mean(dff$P2)
}

dfcsv<-data.frame(df)
colnames(dfcsv)<-c("r","x2_0","K2","lambda","alphaF","G","H","x1_end","x2_end","x1_cv","x2_cv","P1_end","P2_end")

write_csv(dfcsv,paste0("G",G,"_separateK_H",H,".csv"))

# preparing dataset for visualization
sepK<-read_csv(paste0("G",G,"_separateK_H",H,".csv"))

# sepK$outcome<-99
# 
# for (i in 1:nrow(sepK)){
#   
#   if (is.na(sepK[i,]$x1_end)==TRUE){
#     sepK[i,]$x1_end<-0}
#   if (is.na(sepK[i,]$x2_end)==TRUE){
#     sepK[i,]$x2_end<-0}
#   
#     
#   if (sepK[i,]$x1_end<1&sepK[i,]$x2_end<1){
#     sepK[i,]$outcome<-0
#   } else if (sepK[i,]$x1_end>=1&sepK[i,]$x2_end<1){
#     sepK[i,]$outcome<-1
#   } else {sepK[i,]$outcome<-2}
#   
# }

sepK$x1_end2<-NA
sepK$x2_end2<-NA
sepK$selection<-99

for (i in 1:nrow(sepK)){

  if (sepK[i,]$x1_end>=1){
    sepK[i,]$x1_end2<-sepK[i,]$x1_end} else {
      sepK[i,]$x1_end2<-0
      }
  if (sepK[i,]$x2_end>=1){
    sepK[i,]$x2_end2<-sepK[i,]$x2_end} else {
      sepK[i,]$x2_end2<-0
      }


  sepK[i,]$selection<-(sepK[i,]$x2_end2/sepK[i,]$x2_0)/(sepK[i,]$x1_end2/200)-1

}

sepK[is.na(sepK$selection)==TRUE,]$selection<-0

hist(sepK$selection)

sepK[sepK$selection==max(sepK$selection),]

write_csv(sepK,paste0("G",G,"_separateK_H",H,".csv"))



#### shared K
# functions
learning_ode_sharedK<-function(t,state,params){
  
  r1<-params["r1"]
  r2<-params["r2"]
  K_total<-params["K_total"]
  lambda1<-params["lambda1"]
  lambda2<-params["lambda2"]
  alphaL1<-0.5+abs(lambda1-0.5)
  alphaL2<-0.5+abs(lambda2-0.5)
  alphaF1<-params["alphaF1"]
  alphaF2<-params["alphaF2"]
  H<-params["H"]
  GL<-params["GL"]
  GF<-params["GF"]
  # fL<-params["fL"]
  # Eff_social_L<-params["Eff_social_L"]
  # Eff_social_F<-params["Eff_social_F"]
  
  x1<-state["x1"]
  x2<-state["x2"]
  P1<-state["P1"]
  P2<-state["P2"]
  
  e1<-x1/(x1+x2)
  e2<-x2/(x1+x2)
  
  # alphaL1_hat<-H*fL*alphaL1+H*(1-fL)*alphaL1*Eff_social_L
  # alphaL2_hat<-H*fL*alphaL2+H*(1-fL)*alphaL2*Eff_social_L
  # 
  # alphaF1_hat<-H*fL*alphaF1+H*(1-fL)*alphaF1*Eff_social_F
  # alphaF2_hat<-H*fL*alphaF2+H*(1-fL)*alphaF2*Eff_social_F
  
  dP1<-alphaL1*(e1*P1+GL*e2*P2)*(lambda1-P1)+alphaF1*(e1*(1-P1)+e2)*(0.5-P1)
  dP2<-alphaL2*(e2*P2+GL*e1*P1)*(lambda2-P2)+alphaF2*(e2*(1-P2)+e1)*(0.5-P2)
  
  dx1<-r1*x1*(1-(x1+x2)/K_total)-H*x1*P1
  dx2<-r2*x2*(1-(x1+x2)/K_total)-H*x2*P2
  
  list(c(dx1,dx2,dP1,dP2))
}

mimic_ode_sharedK<-function(x2,r,K,lambda,alphaF,G,H){
  
  state<-c(x1 = 200, x2 = x2,P1 = 0.5,P2 = 0.5)
  params2<-c(r1=r,r2=r,K = K,lambda1=lambda,lambda2=lambda,alphaF1=alphaF,alphaF2=alphaF,H=H,GL=G)
  time <- seq(0, 10000, by=0.1)
  test<-deSolve::ode(y = state,func = learning_ode_sharedK,parms = params2,times = time)
  end<-data.frame(tail(test,n = 10000))
  # popsize_end<-c(mean(end$x1),mean(end$x2))
  # cv<-c(sd(end$x1)/popsize_end[1],sd(end$x2)/popsize_end[2])*100
  
  write_csv(end,paste0("D:/Desktop/mimicODE/sharedK/rare",x2,"_r",r,"_K",K,"_lambda",lambda,"_alphaF",alphaF,"_G",G,"_H",H,".csv"))
}

# solving ODE
x2<-c(20,100,180)
K<-c(220,300,380,400)
vector_P<-seq(0,0.5,length.out = 20)
vector_F<-seq(0,1,length.out = 20)
baby<-c(0.2,0.3)
H<-1.2
G<-0.5

all_combo<-matrix(unlist(expand.grid(x2,K,vector_P,vector_F,baby)),ncol=5)
parameter_space<-all_combo[all_combo[,2]-200==all_combo[,1]|all_combo[,2]-200==200,]

clust<-makeCluster(no.Cores, setup_timeout = 0.5)

clusterEvalQ(clust,library(tidyverse))
clusterExport(clust,"learning_ode_sharedK")
clusterExport(clust,"mimic_ode_sharedK")
clusterExport(clust,"parameter_space")
clusterExport(clust,"H")
clusterExport(clust,"G")

parApply(cl = clust,parameter_space,1,function(y){mimic_ode_sharedK(x2 = y[1],r = y[5],K = y[2],lambda = y[3],alphaF = y[4],G = G,H = H)})

stopCluster(clust)
####

# Run the parameter combos that were left out
files<-list.files("D:/Desktop/mimicODE/sharedK")

combos<-NULL

for (i in 1:4800){
  
  combos[i]<-(paste0("rare",parameter_space[i,1],"_r",parameter_space[i,5],"_K",parameter_space[i,2],"_lambda",parameter_space[i,3],"_alphaF",parameter_space[i,4],"_G",G,"_H",H,".csv"))
  
  
}

missing<-setdiff(combos,files)

parameter_space_missing<-matrix(999,nrow = length(missing),ncol = 5)

for (i in 1:length(missing)){
  
  parameter_space_missing[i,1]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
  parameter_space_missing[i,2]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][3],pattern = "K")[[1]][2])
  parameter_space_missing[i,3]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][4],pattern = "lambda")[[1]][2]) 
  parameter_space_missing[i,4]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][5],pattern = "alphaF")[[1]][2])
  parameter_space_missing[i,5]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][2],pattern = "r")[[1]][2])
}

clust<-makeCluster(no.Cores, setup_timeout = 0.5)

clusterEvalQ(clust,library(tidyverse))
clusterExport(clust,"learning_ode_sharedK")
clusterExport(clust,"mimic_ode_sharedK")
clusterExport(clust,"parameter_space_missing")
clusterExport(clust,"H")
clusterExport(clust,"G")

parApply(cl = clust,parameter_space_missing,1,function(y){mimic_ode_sharedK(x2 = y[1],r = y[5],K = y[2],lambda = y[3],alphaF = y[4],G = G,H = H)})

stopCluster(clust)

# summarizing results
G<-0.5
H<-1.2

files<-list.files(paste0("D:/Desktop/mimicODE/sharedK/G",G,"/H",H))
df<-matrix(999,nrow = length(files),ncol = 13)

for (i in 1:length(files)){
  
  dff<-read_csv(paste0("D:/Desktop/mimicODE/sharedK/G",G,"/H",H,"/",files[i]))
  df[i,1]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][2],pattern = "r")[[1]][2])
  df[i,2]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
  df[i,3]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][3],pattern = "K")[[1]][2])
  df[i,4]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][4],pattern = "lambda")[[1]][2])
  df[i,5]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][5],pattern = "alphaF")[[1]][2])
  df[i,6]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][6],pattern = "G")[[1]][2])
  df[i,7]<-H
  df[i,8]<-mean(dff$x1)
  df[i,9]<-mean(dff$x2)
  df[i,10]<-sd(dff$x1)/mean(dff$x1)
  df[i,11]<-sd(dff$x2)/mean(dff$x2)
  df[i,12]<-mean(dff$P1)
  df[i,13]<-mean(dff$P2)
}

dfcsv<-data.frame(df)
colnames(dfcsv)<-c("r","x2_0","K","lambda","alphaF","G","H","x1_end","x2_end","x1_cv","x2_cv","P1_end","P2_end")

write_csv(dfcsv,paste0("G",G,"_sharedK_H",H,".csv"))

# preparing dataset for visualization
sharedK<-read_csv(paste0("G",G,"_sharedK_H",H,".csv"))

# sharedK$outcome<-99
# 
# for (i in 1:nrow(sharedK)){
#   
#   if (is.na(sharedK[i,]$x1_end)==TRUE){
#     sharedK[i,]$x1_end<-0}
#   if (is.na(sharedK[i,]$x2_end)==TRUE){
#     sharedK[i,]$x2_end<-0}
#   
#   if (sharedK[i,]$x1_end<1&sharedK[i,]$x2_end<1){
#     sharedK[i,]$outcome<-0
#   } else if (sharedK[i,]$x1_end>=1&sharedK[i,]$x2_end<1){
#     sharedK[i,]$outcome<-1
#   } else {sharedK[i,]$outcome<-2}
#   
# }
# write_csv(sharedK,paste0("G",G,"_sharedK_H",H,".csv"))


sharedK$x1_end2<-NA
sharedK$x2_end2<-NA
sharedK$selection<-99

for (i in 1:nrow(sharedK)){
  
  if (sharedK[i,]$x1_end>=1){
    sharedK[i,]$x1_end2<-sharedK[i,]$x1_end} else {
      sharedK[i,]$x1_end2<-0
    }
  if (sharedK[i,]$x2_end>=1){
    sharedK[i,]$x2_end2<-sharedK[i,]$x2_end} else {
      sharedK[i,]$x2_end2<-0
    }
  
  
  sharedK[i,]$selection<-(sharedK[i,]$x2_end2/sharedK[i,]$x2_0)/(sharedK[i,]$x1_end2/200)-1
  
}

sharedK[is.na(sharedK$selection)==TRUE,]$selection<-0

hist(sharedK$selection)

sharedK[sharedK$selection==max(sharedK$selection),]

unique(sharedK$outcome)

write_csv(sharedK,paste0("G",G,"_sharedK_H",H,".csv"))


# visualization
# main figures: rare signal is constrained to be rare

G<-0.5

# separate K
sepK_H0.6<-read_csv(paste0("G",G,"_separateK_H0.6.csv"))
sepK_H1.2<-read_csv(paste0("G",G,"_separateK_H1.2.csv"))
sepK<-rbind(sepK_H0.6,sepK_G0.5_H1.2)
sepK$outcome<-as.factor(sepK$outcome)
sepK$scenario<-paste0("x2 = ",sepK$x2_0,", K2 = ",sepK$K2)
sepK2<-subset(sepK,subset = scenario!= "x2 = 20, K2 = 100")

sepK2$priority<-NA

for (i in 1:nrow(sepK2)){
  
if (sepK2[i,]$scenario=="x2 = 20, K2 = 20"){sepK2[i,]$priority<-"A"
  } else if (sepK2[i,]$scenario=="x2 = 20, K2 = 200"){sepK2[i,]$priority<-"B"
  } else if (sepK2[i,]$scenario=="x2 = 100, K2 = 100"){sepK2[i,]$priority<-"C"
  } else if (sepK2[i,]$scenario=="x2 = 100, K2 = 200"){sepK2[i,]$priority<-"D"
  } else if (sepK2[i,]$scenario=="x2 = 180, K2 = 180"){sepK2[i,]$priority<-"E"
  } else if(sepK2[i,]$scenario=="x2 = 180, K2 = 200"){sepK2[i,]$priority<-"F"}

}

scenario_list<-c(A = "x2 = 20, K2 = 20",
                 B = "x2 = 20, K2 = 200",
                 C = "x2 = 100, K2 = 100",
                 D = "x2 = 100, K2 = 200",
                 E = "x2 = 180, K2 = 180",
                 F = "x2 = 180, K2 = 200")

sepK_outcome_r0.2_H0.6<-ggplot(subset(sepK2,subset = r==0.2 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.2_H1.2<-ggplot(subset(sepK2,subset = r==0.2 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.3_H0.6<-ggplot(subset(sepK2,subset = r==0.3 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.3_H1.2<-ggplot(subset(sepK2,subset = r==0.3 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sepK_outcome_r0.2_H0.6,sepK_outcome_r0.2_H1.2,sepK_outcome_r0.3_H0.6,sepK_outcome_r0.3_H1.2,ncol=1)


sepK_sel_r0.2_H0.6<-ggplot(subset(sepK2,subset = r==0.2 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),breaks=c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
sepK_sel_r0.2_H1.2<-ggplot(subset(sepK2,subset = r==0.2 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),breaks=c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_sel_r0.3_H0.6<-ggplot(subset(sepK2,subset = r==0.3 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),breaks=c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_sel_r0.3_H1.2<-ggplot(subset(sepK2,subset = r==0.3 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),breaks=c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sepK_sel_r0.2_H0.6,sepK_sel_r0.2_H1.2,sepK_sel_r0.3_H0.6,sepK_sel_r0.3_H1.2,ncol=1)

# supplementary figures: rare signal could be as abundant as the common signal
sepK_outcome_r0.2_H0.6<-ggplot(subset(sepK2,subset = r==0.2 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.2_H1.2<-ggplot(subset(sepK2,subset = r==0.2 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.3_H0.6<-ggplot(subset(sepK2,subset = r==0.3 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_outcome_r0.3_H1.2<-ggplot(subset(sepK2,subset = r==0.3 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sepK_outcome_r0.2_H0.6,sepK_outcome_r0.2_H1.2,sepK_outcome_r0.3_H0.6,sepK_outcome_r0.3_H1.2,ncol=1)


sepK_sel_r0.2_H0.6<-ggplot(subset(sepK2,subset = r==0.2 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),breaks=c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_sel_r0.2_H1.2<-ggplot(subset(sepK2,subset = r==0.2 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),breaks=c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_sel_r0.3_H0.6<-ggplot(subset(sepK2,subset = r==0.3 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),breaks=c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sepK_sel_r0.3_H1.2<-ggplot(subset(sepK2,subset = r==0.3 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),breaks=c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sepK_sel_r0.2_H0.6,sepK_sel_r0.2_H1.2,sepK_sel_r0.3_H0.6,sepK_sel_r0.3_H1.2,ncol=1)


# shared K
G<-0.5

sharedK_H0.6<-read_csv(paste0("G",G,"_sharedK_H0.6.csv"))
sharedK_H1.2<-read_csv(paste0("G",G,"_sharedK_H1.2.csv"))
sharedK<-rbind(sharedK_H0.6,sharedK_H1.2)
sharedK$outcome<-as.factor(sharedK$outcome)
sharedK$scenario<-paste0("x2 = ",sharedK$x2_0,", K = ",sharedK$K)
sharedK2<-subset(sharedK,subset = scenario!= "x2 = 20, K = 300")

sharedK2$priority<-NA

for (i in 1:nrow(sharedK2)){
  
  if (sharedK2[i,]$scenario=="x2 = 20, K = 220"){sharedK2[i,]$priority<-"A"
  } else if (sharedK2[i,]$scenario=="x2 = 20, K = 400"){sharedK2[i,]$priority<-"B"
  } else if (sharedK2[i,]$scenario=="x2 = 100, K = 300"){sharedK2[i,]$priority<-"C"
  } else if (sharedK2[i,]$scenario=="x2 = 100, K = 400"){sharedK2[i,]$priority<-"D"
  } else if (sharedK2[i,]$scenario=="x2 = 180, K = 380"){sharedK2[i,]$priority<-"E"
  } else if(sharedK2[i,]$scenario=="x2 = 180, K = 400"){sharedK2[i,]$priority<-"F"}
  
}

scenario_list<-c(A = "x2 = 20, K = 220",
                 B = "x2 = 20, K = 400",
                 C = "x2 = 100, K = 300",
                 D = "x2 = 100, K = 400",
                 E = "x2 = 180, K = 380",
                 F = "x2 = 180, K = 400")


# main figures: rare signal was constrained to be rare
sharedK_outcome_r0.2_H0.6<-ggplot(subset(sharedK2,subset = r==0.2 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.2, G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.2_H1.2<-ggplot(subset(sharedK2,subset = r==0.2 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.2, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.3_H0.6<-ggplot(subset(sharedK2,subset = r==0.3 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.3, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.3_H1.2<-ggplot(subset(sharedK2,subset = r==0.3 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.3, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sharedK_outcome_r0.2_H0.6,sharedK_outcome_r0.2_H1.2,sharedK_outcome_r0.3_H0.6,sharedK_outcome_r0.3_H1.2,ncol=1)


sharedK_sel_r0.2_H0.6<-ggplot(subset(sharedK2,subset = r==0.2 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.2_H1.2<-ggplot(subset(sharedK2,subset = r==0.2 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.3_H0.6<-ggplot(subset(sharedK2,subset = r==0.3 & H==0.6 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.3_H1.2<-ggplot(subset(sharedK2,subset = r==0.3 & H==1.2 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-1,1),c(seq(-1,1,by = 0.2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sharedK_sel_r0.2_H0.6,sharedK_sel_r0.2_H1.2,sharedK_sel_r0.3_H0.6,sharedK_sel_r0.3_H1.2,ncol=1)

# supplementary figures: rare signal could be as abundant as the common signal
sharedK_outcome_r0.2_H0.6<-ggplot(subset(sharedK2,subset = r==0.2 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.2, G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.2_H1.2<-ggplot(subset(sharedK2,subset = r==0.2 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.2, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.3_H0.6<-ggplot(subset(sharedK2,subset = r==0.3 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.3, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_r0.3_H1.2<-ggplot(subset(sharedK2,subset = r==0.3 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_tile()+
  # ggtitle(paste0("Shared K, r0 = 0.3, G = ",G, ", H = 1.2"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_manual(values = c("black","red","blue"),breaks = c(0,1,2))+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sharedK_outcome_r0.2_H0.6,sharedK_outcome_r0.2_H1.2,sharedK_outcome_r0.3_H0.6,sharedK_outcome_r0.3_H1.2,ncol=1)


sharedK_sel_r0.2_H0.6<-ggplot(subset(sharedK2,subset = r==0.2 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.2_H1.2<-ggplot(subset(sharedK2,subset = r==0.2 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.3_H0.6<-ggplot(subset(sharedK2,subset = r==0.3 & H==0.6 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_r0.3_H1.2<-ggplot(subset(sharedK2,subset = r==0.3 & H==1.2 & priority %in% c("B","D","F")),aes(x = lambda,y = alphaF,fill=round(selection,3)))+
  geom_tile(stat = "identity")+
  # ggtitle(paste0("Shared K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-9,9),c(seq(-8,8,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

gridExtra::grid.arrange(sharedK_sel_r0.2_H0.6,sharedK_sel_r0.2_H1.2,sharedK_sel_r0.3_H0.6,sharedK_sel_r0.3_H1.2,ncol=1)


# Examining the relationship between the persistence of the rare signal and selection acting on it
# Separate K
sepK_G0_H0.6<-read_csv("G0_separateK_H0.6.csv")
sepK_G0_H1.2<-read_csv("G0_separateK_H1.2.csv")
sepK_G0.5_H0.6<-read_csv("G0.5_separateK_H0.6.csv")
sepK_G0.5_H1.2<-read_csv("G0.5_separateK_H1.2.csv")

sepK_all<-rbind(sepK_G0_H0.6,sepK_G0_H1.2,sepK_G0.5_H0.6,sepK_G0.5_H1.2)
head(sepK_all)

hist_S_sepK<-ggplot(data = sepK_all[sepK_all$outcome==2,])+
  geom_histogram(aes(x = selection))+
  mdthemes::md_theme_bw()+
  labs(x = "Selection on rare signal (S)", y = "Count")+
  ggtitle("(A) Separate carrying capacities")

# Shared K
sharedK_G0_H0.6<-read_csv("G0_sharedK_H0.6.csv")
sharedK_G0_H1.2<-read_csv("G0_sharedK_H1.2.csv")
sharedK_G0.5_H0.6<-read_csv("G0.5_sharedK_H0.6.csv")
sharedK_G0.5_H1.2<-read_csv("G0.5_sharedK_H1.2.csv")

sharedK_all<-rbind(sharedK_G0_H0.6,sharedK_G0_H1.2,sharedK_G0.5_H0.6,sharedK_G0.5_H1.2)
head(sharedK_all)

hist_S_sharedK<-ggplot(data = sharedK_all[sharedK_all$outcome==2,])+
  geom_histogram(aes(x = selection))+
  mdthemes::md_theme_bw()+
  labs(x = "Selection on rare signal (S)", y = "Count")+
  ggtitle("(B) Shared carrying capacity")

gridExtra::grid.arrange(hist_S_sepK,hist_S_sharedK,ncol = 2)
