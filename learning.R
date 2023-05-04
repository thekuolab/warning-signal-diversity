m(list=ls())

# ### Run simulation by varying P and F
# source("learning_functions.r")
# 
# vector_P<-seq(0.02,0.2,length.out = 20)
# vector_F<-seq(0.01,0.1,length.out = 20)
# 
# parameter_space<-expand.grid(vector_P,vector_F)
# 
# no.Cores<-detectCores()
# 
# for (i in 1:15){
#   
#   # # When run in mac
#   # dir<-paste0("/Volumes/KuoLab-1/learning/ver2/P_by_F/",i)
#   
#   # When run in Windows
#   dir<-paste0("E:/learning/ver2/P_by_F_400/",i)
#   
#   #dir.create(dir)
#   
#   setwd(dir)
#   
#   clust<-makeCluster(no.Cores, setup_timeout = 0.5)
#   
#   clusterEvalQ(clust,library(tidyverse))
#   clusterExport(clust,"learning")
#   clusterExport(clust,"mimic")
#   clusterExport(clust,"parameter_space")
#   
#   
#   
#   parApply(cl = clust,parameter_space,1,function(y){mimic(G = 0,P_OLD = y[1],Abundance = c(200,20),
#                                                           lambda_L = c(y[1],y[1]),alpha_F = y[2],baby = 0.2,K = c(200,20),
#                                                           GenTime = 100,timestep = 50000)})
#   
#   stopCluster(clust)
#   
#   # gm_send_message(email)
#   
# }
# 
# ### Summarizing P_by_F simulation results 
# rm(list=ls())
# setwd("~/Google Drive/learning/P_by_F_20by20_G0.5")
# 
# library(tidyverse)
# library(ggplot2)
# 
# replication<-20
# 
# files<-list.files(paste0("~/Google Drive/learning/P_by_F_20by20_G0.5/",replication))
# setwd(paste0("~/Google Drive/learning/P_by_F_20by20_G0.5/",replication))
# 
# P_by_F<-matrix(0,nrow=length(files),ncol=5)
# 
# 
# for (i in 1:length(files)){
#   
#   # string<-files[i]
#   # first<-str_extract(string = string,pattern = "G[0-1]+")
#   # G<-as.numeric(substring(first,2:))
#   # second<-str_extract(string = string,pattern = "NEW[0-9]+")
#   # abnc_new_t0<-as.numeric(substring(second,4))
#   
#   data<-read_csv(files[i])
#   P<-data$P_OLD[1]
#   abnc_OLD_final<-data[50000,]$abnc_OLD
#   abnc_new_final<-data[50000,]$abnc_NEW
#   
#   P_by_F[i,]<-c(P,0,abnc_OLD_final,abnc_new_final,-99)
#   
# }
# 
# # P_by_F[,2]<-rep(seq(0.01,0.1,length.out = 10),10)
# P_by_F[,2]<-rep(seq(0.01,0.1,length.out = 20),20)
# 
# P_by_F<-data.frame(P_by_F)
# colnames(P_by_F)<-c("P","F","abnc_OLD_final","abnc_NEW_final","outcome")
# 
# for (i in 1:nrow(P_by_F)){
#   
#   if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
#     P_by_F[i,]$outcome<-0
#   } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
#     P_by_F[i,]$outcome<-1
#   } else {P_by_F[i,]$outcome<-2}
#   
# }
# 
# P_by_F$selection_NEW<-(P_by_F$abnc_NEW_final/20)/(P_by_F$abnc_OLD_final/200)-1
# 
# P_by_F[is.na(P_by_F$selection_NEW),]$selection_NEW<--1
# 
# write_csv(P_by_F,paste0("~/Google Drive/learning/P_by_F_20by20_G0.5/P_by_F_rep",replication,".csv"))
# 
# # Heatmap of P x F based on outcome
# 
# ggplot(P_by_F,aes(x = P,y = F,fill=factor(outcome)))+
#   geom_tile()+
#   theme_minimal()+
#   scale_fill_manual(values=c("black","red","blue"))
# # scale_fill_gradientn(colors=c("black","#b5d1e2","red"),values=c(0,1,2))
# 
# # Heatmap of P x F based on selection against the NEW prey
# 
# ggplot(P_by_F,aes(x = P,y = F,fill=selection_NEW))+
#   geom_tile()+
#   scale_fill_distiller(palette = "YlGnBu")+
#   theme_minimal()  
# 
# ### Summarizing simulation results
# 
# rm(list=ls())
# 
# setwd("~/Google Drive/learning/")
# 
# source("learning_functions.r")
# 
# files<-list.files("~/Google Drive/learning/P_by_F_20by20_G0.5/summary")
# 
# outcome<-matrix(-99,nrow = 400,ncol = 20)
# selection_NEW<-matrix(-99,nrow = 400,ncol = 20)
# 
# for (i in 1:length(files)){
#   # for (i in 1:1){
#   
#   df<-read_csv(paste0("~/Google Drive/learning/P_by_F_20by20_G0.5/summary/",files[i]))
#   
#   outcome[,i]<-df$outcome
#   selection_NEW[,i]<-df$selection_NEW
#   
# }
# 
# P_by_F<-data.frame(matrix(0,nrow = 400,ncol = 6))
# colnames(P_by_F)<-c("P","F","outcome","outcome_sd","selection_NEW","selection_NEW_sd")
# 
# P_by_F[,1:2]<-df[,1:2]
# 
# for (i in 1:400){
#   
#   P_by_F[i,3]<-mean(outcome[i,1:20])
#   P_by_F[i,4]<-sd(outcome[i,1:20])
#   P_by_F[i,5]<-mean(selection_NEW[i,1:20])
#   P_by_F[i,6]<-sd(selection_NEW[i,1:20])
#   
# }
# 
# write_csv(P_by_F,"P_by_F_G0.5.csv")
# 
# ggplot(P_by_F,aes(x = P,y = F,fill=selection_NEW))+
#   geom_tile()+
#   scale_fill_distiller(palette = "YlGnBu")+
#   theme_minimal()
# 
# 
# ggplot(P_by_F,aes(x = P,y = F,fill=outcome))+
#   geom_tile()+
#   theme_minimal()+
#   # scale_color_gradientn(colors=c("black","#b5d1e2","red"),values=c(0,1,2))
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
# 
# 
# ## Sensitivity analysis on R0 (baby in the function) on server
# ## Use learning_for_server.R for code
# 
# # Summarizing sensitivity analyses results
# rm(list=ls())
# 
# library(tidyverse)
# library(ggplot2)
# 
# replication<-20
# 
# # dir<-paste0("/Volumes/KuoLab-1/learning/ver2/sensitivity/",replication)
# # dir<-paste0("~/Google Drive/learning/sensitivity/G0.5/",replication)
# 
# dir<-paste0("~/Desktop/",replication)
# 
# setwd(dir)
# 
# files<-list.files(dir)
# 
# P_by_F<-matrix(0,nrow=length(files),ncol=5)
# 
# baby<-NULL
# 
# for (i in 1:length(files)){
#   
#   # Getting baby value from file names
#   filename<-files[i]
#   # This step creates a "baby value.csv" string
#   a<-strsplit(filename,"_baby_")[[1]][2]
#   # This step isolates the baby value from the ".csv" and turns it into a number
#   baby[i]<-as.numeric(str_remove(a,".csv"))
#   # b would be the baby value used in the simulation
#   
#   data<-read_csv(files[i])
#   P<-data$P_OLD[1]
#   abnc_OLD_final<-data[50000,]$abnc_OLD
#   abnc_new_final<-data[50000,]$abnc_NEW
#   
#   P_by_F[i,]<-c(P,0,abnc_OLD_final,abnc_new_final,-99)
#   
# }
# 
# P_by_F[,2]<-rep(rep(seq(0.01,0.1,length.out = 10),each = 10),10)
# 
# P_by_F<-data.frame(P_by_F)
# colnames(P_by_F)<-c("P","F","abnc_OLD_final","abnc_NEW_final","outcome")
# 
# P_by_F$baby<-baby
# 
# for (i in 1:nrow(P_by_F)){
#   
#   if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
#     P_by_F[i,]$outcome<-0
#   } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
#     P_by_F[i,]$outcome<-1
#   } else {P_by_F[i,]$outcome<-2}
#   
# }
# 
# P_by_F$selection_NEW<-(P_by_F$abnc_NEW_final/20)/(P_by_F$abnc_OLD_final/200)-1
# 
# P_by_F[is.na(P_by_F$selection_NEW),]$selection_NEW<--1
# 
# write_csv(P_by_F,paste0("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5/",replication,".csv"))
# 
# # Plotting the outcomes for different baby values
# library(tidyverse)
# 
# rm(list=ls())
# 
# replication<-20
# 
# setwd("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5")
# 
# P_by_F<-read_csv(paste0(replication,".csv"))
# 
# baby<-unique(P_by_F$baby)
# 
# # Heatmap of P x F based on outcome
# 
# # Create two lists to store ggplot objects in the loop
# plot1_list<-list(NA)
# plot2_list<-list(NA)
# 
# for (i in 1:length(baby)){
#   
#   df<-subset(P_by_F,subset = P_by_F$baby==baby[i])
#   
#   plot1_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=factor(outcome)))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+
#     theme(legend.position = "none")+
#     scale_fill_manual(values=c("black","red","blue"),breaks = c(0,1,2))
#   
#   # Heatmap of P x F based on selection against the NEW prey
#   
#   plot2_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=selection_NEW))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+  
#     theme(legend.position = "none")+
#     scale_fill_distiller(palette = "YlGnBu")
#   
# }
# 
# library(gridExtra)
# 
# grid.arrange(plot1_list[[1]],plot1_list[[2]],plot1_list[[3]],plot1_list[[4]],
#              plot1_list[[5]],plot1_list[[6]],plot1_list[[7]],plot1_list[[8]],
#              plot1_list[[9]],plot1_list[[10]],ncol=5)
# 
# grid.arrange(plot2_list[[1]],plot2_list[[2]],plot2_list[[3]],plot2_list[[4]],
#              plot2_list[[5]],plot2_list[[6]],plot2_list[[7]],plot2_list[[8]],
#              plot2_list[[9]],plot2_list[[10]],ncol=5)
# 
# 
# # Summarizing results from different replications
# rm(list=ls())
# 
# # Change the working directory as appropriate
# setwd("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5")
# 
# files<-list.files("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5")
# 
# outcome<-matrix(-99,nrow = 1000,ncol = 20)
# selection_NEW<-matrix(-99,nrow = 1000,ncol = 20)
# 
# for (i in 1:length(files)){
#   # for (i in 1:1){
#   
#   df<-read_csv(paste0("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5/",files[i]))
#   
#   outcome[,i]<-df$outcome
#   selection_NEW[,i]<-df$selection_NEW
#   
# }
# 
# P_by_F<-data.frame(matrix(0,nrow = 1000,ncol = 7))
# colnames(P_by_F)<-c("P","F","baby","outcome","outcome_sd","selection_NEW","selection_NEW_sd")
# 
# P_by_F[,1:3]<-df[,c(1,2,6)]
# 
# for (i in 1:1000){
#   
#   # Change X in outcome[i,1:X] as new replications are done
#   P_by_F[i,4]<-mean(outcome[i,1:20])
#   P_by_F[i,5]<-sd(outcome[i,1:20])
#   P_by_F[i,6]<-mean(selection_NEW[i,1:20])
#   P_by_F[i,7]<-sd(selection_NEW[i,1:20])
#   
# }
# 
# baby<-unique(P_by_F$baby)
# 
# # Heatmap of P x F based on outcome and selection against the new prey
# 
# # Create two lists to store ggplot objects in the loop
# plot1_list<-list(NA)
# plot2_list<-list(NA)
# 
# for (i in 1:length(baby)){
#   
#   df<-subset(P_by_F,subset = P_by_F$baby==baby[i])
#   
#   plot1_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=outcome))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+
#     theme(legend.position = "none")+
#     scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
#   
#   # Heatmap of P x F based on selection against the NEW prey
#   
#   plot2_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=selection_NEW))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+  
#     theme(legend.position = "none")+
#     scale_fill_distiller(palette = "YlGnBu")
#   
# }
# 
# library(gridExtra)
# 
# grid.arrange(plot1_list[[1]],plot1_list[[2]],plot1_list[[3]],plot1_list[[4]],
#              plot1_list[[5]],plot1_list[[6]],plot1_list[[7]],plot1_list[[8]],
#              plot1_list[[9]],plot1_list[[10]],ncol=5)
# 
# grid.arrange(plot2_list[[1]],plot2_list[[2]],plot2_list[[3]],plot2_list[[4]],
#              plot2_list[[5]],plot2_list[[6]],plot2_list[[7]],plot2_list[[8]],
#              plot2_list[[9]],plot2_list[[10]],ncol=5)
# 
# write_csv(P_by_F,"G0.5.csv")
# 
# # Plot the ecological outcomes and selection against the new signal
# # Only using baby from 0.122 - 0.2
# # For G = 0 and G = 0.5
# 
# plot1_list<-list(NA)
# plot2_list<-list(NA)
# 
# setwd("/Volumes/GoogleDrive/My Drive/Predator learning and warning signal diversity/sensitivity/G0.5")
# P_by_F<-read_csv("G0.5.csv")
# 
# baby<-unique(P_by_F$baby)
# 
# for (i in 1:length(baby)){
#   
#   df<-subset(P_by_F,subset = P_by_F$baby==baby[i])
  
#   plot1_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=outcome))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+
#     theme(legend.position = "none")+
#     scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
#   
#   # Heatmap of P x F based on selection against the NEW prey
#   
#   plot2_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=selection_NEW))+
#     geom_tile()+
#     ggtitle(paste0("baby = ",baby[i]))+
#     theme_minimal()+  
#     theme(legend.position = "none")+
#     scale_fill_distiller(palette = "YlGnBu")
#   
# }
# 
# P_by_F_baby0.2<-read_csv("/Volumes/GoogleDrive/My Drive/Predator learning and warning signal diversity/P_by_F_20by20_G0.5/P_by_F_G0.5.csv")
# 
# plot1<-ggplot(P_by_F_baby0.2,aes(x = P,y = F,fill=outcome))+
#   geom_tile()+
#   ggtitle(paste0("baby = 0.2"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
# 
# plot2<-ggplot(P_by_F_baby0.2,aes(x = P,y = F,fill=selection_NEW))+
#   geom_tile()+
#   ggtitle(paste0("baby = 0.2"))+
#   theme_minimal()+  
#   theme(legend.position = "none")+
#   scale_fill_distiller(palette = "YlGnBu")
# 
# library(gridExtra)
# 
# grid.arrange(plot1_list[[2]],plot1_list[[3]],plot1_list[[4]],plot1_list[[5]],plot1,
#              plot2_list[[2]],plot2_list[[3]],plot2_list[[4]],plot2_list[[5]],plot2,ncol=5)
# 
# rm(list=ls())
# 
# # Plot the best- and worst-case scenario with added data points
# library(tidyverse)
# df<-read_csv("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0/G0.csv")
# 
# G0baby0.122<-subset(df,subset = df$baby==0.122)
# 
# G0baby0.2<-read_csv("~/Google Drive/Predator learning and warning signal diversity/P_by_F_20by20_G0/P_by_F_G0.csv")
# 
# Ps<-c(0.02,0.2,0.05,0.02,0.08)
# Fs<-c(0.01,0.06,0.01,0.01,0.035)
# data<-data.frame(cbind(Ps,Fs))
# 
# plot1<-ggplot(G0baby0.122,aes(x = P,y = F))+
#   geom_tile(aes(fill=outcome))+
#   # ggtitle(paste0("baby = 0.122"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
#   geom_point(data = data,mapping = aes(x=Ps,y=Fs,size=3),shape=21,color="black",fill="white")
# 
# plot2<-ggplot(G0baby0.2,aes(x = P,y = F))+
#   geom_tile(aes(fill=outcome))+
#   # ggtitle(paste0("baby = 0.2"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
#   geom_point(data = data,mapping = aes(x=Ps,y=Fs,size=3),shape=21,color="black",fill="white")
# 
# df2<-read_csv("~/Google Drive/Predator learning and warning signal diversity/sensitivity/G0.5/G0.5.csv")
# 
# G0.5baby0.122<-subset(df2,subset = df$baby==0.122)
# 
# G0.5baby0.2<-read_csv("~/Google Drive/Predator learning and warning signal diversity/P_by_F_20by20_G0.5/P_by_F_G0.5.csv")
# 
# plot3<-ggplot(G0.5baby0.122,aes(x = P,y = F))+
#   geom_tile(aes(fill=outcome))+
#   # ggtitle(paste0("baby = 0.122"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
#   geom_point(data = data,mapping = aes(x=Ps,y=Fs,size=3),shape=21,color="black",fill="white")
# 
# plot4<-ggplot(G0.5baby0.2,aes(x = P,y = F))+
#   geom_tile(aes(fill=outcome))+
#   # ggtitle(paste0("baby = 0.2"))+
#   theme_minimal()+
#   theme(legend.position = "none")+
#   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
#   geom_point(data = data,mapping = aes(x=Ps,y=Fs,size=3),shape=21,color="black",fill="white")
# 
# gridExtra::grid.arrange(plot1,plot2,plot3,plot4,ncol=2)


# Social learning (perennial prey)
library(tidyverse)
library(parallel)

source("learning_functions.R")

vector_P<-seq(0,0.5,length.out = 10)
vector_F<-seq(0,0.5,length.out = 10)

PbyF<-expand.grid(vector_P,vector_F)

v_PbyF<-PbyF %>% slice(rep(1:n(),each = 1))

baby<-rep(0.28,times=100)

parameter_space<-cbind(v_PbyF,baby)

no.Cores<-detectCores()

# For perennial prey, use mimic_social_fast
# For annual prey, use mimic_social_fast2

for (i in 9:20){
  
  dir<-paste0("D:/Desktop/",i)
  
  dir.create(dir)
  
  setwd(dir)
  
  clust<-makeCluster(no.Cores, setup_timeout = 0.5)
  
  clusterEvalQ(clust,library(tidyverse))
  clusterExport(clust,"learning")
  clusterExport(clust,"mimic_social_fast")
  clusterExport(clust,"parameter_space")
  
  
  
  parApply(cl = clust,parameter_space,1,function(y){mimic_social_fast(G = 0,P_OLD = y[1],Abundance = c(200,10),
                                                                 lambda_L = c(y[1],y[1]),alpha_F = y[2],baby = y[3],
                                                                 K = c(200,10),social_spread = 0.5,social_efficacy = 0.5,
                                                                 inter_learning = TRUE,
                                                                 GenTime = 100,timestep = 50000)})
  
  stopCluster(clust)
  
}

# Summarizing social learning results
library(tidyverse)
library(gridExtra)

rm(list=ls())

K<-c(200,10)

G<-0

# baby<-c(0.19,0.21,0.24,0.26,0.28)
baby<-0.28

for (j in 1:20){
  
  replication<-j
  # On Windows
  # dir<-paste0("G:/我的雲端硬碟/Predator learning and warning signal diversity/social_annualprey/K",K[1],"_",K[2],"/",replication)
  
  # On Mac
  # dir<-paste0("~/Desktop/K200_",K[2],"/",replication)
  # dir<-paste0("~/Desktop/K200_",K[2],"_G",G,"/",replication)
  dir<-paste0("~/Desktop/K200_",K[2],"_social0.5/",replication)
  
  setwd(dir)
  
  files<-unique(list.files(dir))
  
  P_by_F<-matrix(0,nrow=length(files),ncol=6)
  
  P<-seq(0,0.5,length=10)
  F<-seq(0,0.5,length=10)
  
  # P_by_F[,1]<-rep(P,each=50)
  # P_by_F[,2]<-rep(F,each=5)
  # P_by_F[,3]<-rep(baby)
  
  P_by_F[,1]<-rep(P,each=10)
  P_by_F[,2]<-rep(F,10)
  P_by_F[,3]<-baby
  
  for (i in 1:length(files)){
    
    # string<-files[i]
    # first<-str_extract(string = string,pattern = "G[0-1]+")
    # G<-as.numeric(substring(first,2:))
    # second<-str_extract(string = string,pattern = "NEW[0-9]+")
    # abnc_new_t0<-as.numeric(substring(second,4))
    
    data<-read_csv(files[i])
    abnc_OLD_final<-data[1,]$abnc_OLD
    abnc_new_final<-data[1,]$abnc_NEW
    
    P_by_F[i,c(4:6)]<-c(abnc_OLD_final,abnc_new_final,-99)
    
  }
  
  P_by_F<-data.frame(P_by_F)
  colnames(P_by_F)<-c("P","F","baby","abnc_OLD_final","abnc_NEW_final","outcome")
  
  for (i in 1:nrow(P_by_F)){
    
    if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-0
    } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-1
    } else {P_by_F[i,]$outcome<-2}
    
  }
  
  P_by_F$selection_NEW<-(P_by_F$abnc_NEW_final/K[2])/(P_by_F$abnc_OLD_final/K[1])-1
  
  P_by_F[is.na(P_by_F$selection_NEW),]$selection_NEW<--1
  
  ggplot(P_by_F,aes(x = P,y = F,fill=outcome))+
      geom_tile()+
      ggtitle(paste0(paste0("r0 = ",baby[1])))+
      theme_minimal()+
      theme(legend.position = "none")+
      scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  
  
  # P_by_F_baby1<-subset(P_by_F,subset = P_by_F$baby==baby[1])
  # P_by_F_baby2<-subset(P_by_F,subset = P_by_F$baby==baby[2])
  # P_by_F_baby3<-subset(P_by_F,subset = P_by_F$baby==baby[3])
  # P_by_F_baby4<-subset(P_by_F,subset = P_by_F$baby==baby[4])
  # P_by_F_baby5<-subset(P_by_F,subset = P_by_F$baby==baby[5])
  # 
  # plot1<-ggplot(P_by_F_baby1,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[1])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot2<-ggplot(P_by_F_baby2,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[2])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot3<-ggplot(P_by_F_baby3,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[3])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot4<-ggplot(P_by_F_baby4,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[4])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot5<-ggplot(P_by_F_baby5,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[5])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # library(gridExtra)
  # 
  # grid.arrange(plot1,plot2,plot3,plot4,plot5,ncol=5)
  
  
  
  # plot2<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=selection_NEW))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[1])))+
  #   theme_minimal()+  
  #   theme(legend.position = "none")+
  #   scale_fill_distiller(palette = "YlGnBu")
  # 
  # 
  # plot4<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=selection_NEW))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[2])))+
  #   theme_minimal()+  
  #   theme(legend.position = "none")+
  #   scale_fill_distiller(palette = "YlGnBu")
  
  
  
  # dir2<-paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_no social learning/baby0.28/")
  dir2<-paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_social learning/baby0.28/")
  
  dir.create(dir2)
  
  setwd(dir2)
  
  write_csv(P_by_F,paste0(replication,".csv"))
  
}

# Summarizing results from different replications
rm(list=ls())

K<-c(200,40)

G<-0

baby<-c(0.19,0.21,0.24,0.26,0.28)
# baby<-0.28

##### Change the working directory as appropriate
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_no social learning/baby_0.19to0.28/"))
# setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_no social learning/baby0.28/"))
# setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_social learning/baby0.28/"))
files<-list.files()

outcome<-matrix(-99,nrow = 100,ncol = 20)
selection_NEW<-matrix(-99,nrow = 100,ncol = 20)

# outcome<-matrix(-99,nrow = 500,ncol = 20)
# selection_NEW<-matrix(-99,nrow = 500,ncol = 20)

for (i in 1:length(files)){
  
  df<-read_csv(files[i])
  
  outcome[,i]<-df$outcome
  selection_NEW[,i]<-df$selection_NEW
  
}

P_by_F<-data.frame(matrix(0,nrow = 100,ncol = 7))
colnames(P_by_F)<-c("P","F","baby","outcome","outcome_sd","selection_NEW","selection_NEW_sd")

P_by_F[,1:3]<-df[,1:3]

# for (i in 1:500){
for (i in 1:100){
  
  # Change X in outcome[i,1:X] as new replications are done
  P_by_F[i,4]<-mean(outcome[i,1:length(files)])
  P_by_F[i,5]<-sd(outcome[i,1:length(files)])
  P_by_F[i,6]<-mean(selection_NEW[i,1:length(files)])
  P_by_F[i,7]<-sd(selection_NEW[i,1:length(files)])
  
}

baby<-unique(P_by_F$baby)

# Heatmap of P x F based on outcome and selection against the new prey

# Run this line when directly reading the summary data frame
P_by_F<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G0_no social learning/baby_0.19to0.28/K200_",K[2],"_baby_0.19to0.28_G0.csv"))

# Create two lists to store ggplot objects in the loop
plot1_list<-list(NA)
plot2_list<-list(NA)

for (i in 1:length(baby)){
  
  df<-subset(P_by_F,subset = P_by_F$baby==baby[i])
  
  plot1_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=outcome))+
    geom_tile()+
    ggtitle(paste0("r0 = ",baby[i]))+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    xlab("")+
    ylab("")
  
  # Heatmap of P x F based on selection against the NEW prey
  
  plot2_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=selection_NEW))+
    geom_tile()+
    ggtitle(paste0("r0 = ",baby[i]))+
    theme_minimal()+  
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "YlGnBu")
  
}

library(gridExtra)

plot1_list[[1]]
plot2_list[[1]]

grid.arrange(plot1_list[[1]],plot1_list[[2]],plot1_list[[3]],plot1_list[[4]],plot1_list[[5]],ncol=5)
# 
# grid.arrange(plot2_list[[1]],plot2_list[[2]],plot2_list[[4]],plot2_list[[4]],plot2_list[[5]],ncol=5)

P_by_F_40<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_40/G0_no social learning/baby_0.19to0.28/K200_40_baby_0.19to0.28_G0.csv"))
P_by_F_20<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_20/G0_no social learning/baby_0.19to0.28/K200_20_baby_0.19to0.28_G0.csv"))
P_by_F_10<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_10/G0_no social learning/baby_0.19to0.28/K200_10_baby_0.19to0.28_G0.csv"))

plot1_list_40<-list(NA)
plot1_list_20<-list(NA)
plot1_list_10<-list(NA)

for (i in 1:length(baby)){
  
  df40<-subset(P_by_F_40,subset = P_by_F_40$baby==baby[i])
  df20<-subset(P_by_F_20,subset = P_by_F_20$baby==baby[i])
  df10<-subset(P_by_F_10,subset = P_by_F_10$baby==baby[i])
  
  plot1_list_40[[i]]<-ggplot(df40,aes(x = P,y = F,fill=outcome))+
    geom_tile()+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    xlab("")+
    ylab("")
  
  plot1_list_20[[i]]<-ggplot(df20,aes(x = P,y = F,fill=outcome))+
    geom_tile()+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    xlab("")+
    ylab("")
  
  plot1_list_10[[i]]<-ggplot(df10,aes(x = P,y = F,fill=outcome))+
    geom_tile()+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    xlab("")+
    ylab("")
}

tiff("figure1.tiff",width = 10,height = 5,units = "in",res = 300,)
gridExtra::grid.arrange(plot1_list_40[[1]],plot1_list_40[[2]],plot1_list_40[[3]],plot1_list_40[[4]],plot1_list_40[[5]],
                        plot1_list_20[[1]],plot1_list_20[[2]],plot1_list_20[[3]],plot1_list_20[[4]],plot1_list_20[[5]],
                        plot1_list_10[[1]],plot1_list_10[[2]],plot1_list_10[[3]],plot1_list_10[[4]],plot1_list_10[[5]],
                        ncol=5)
dev.off()

# G = 0
write_csv(P_by_F,paste0("K200_",K[2],"_baby_0.19to0.28_G",G,".csv"))

# G = 0.5
write_csv(P_by_F,paste0("K200_",K[2],"_G",G,".csv"))

# Social = 0.5
write_csv(P_by_F,paste0("K200_",K[2],"_social0.5.csv"))


# G = 0, G = 0.5, and social = 0.5 together
df<-read_csv("~/Google Drive/My Drive/Predator learning and warning signal diversity/P&F.csv")

K200_40_G0<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_40/G0_no social learning/baby_0.19to0.28/K200_40_baby_0.19to0.28_G0.csv"))
K200_40_G0b<-subset(K200_40_G0,subset = K200_40_G0$baby==0.28)

K200_20_G0<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_20/G0_no social learning/baby_0.19to0.28/K200_20_baby_0.19to0.28_G0.csv"))
K200_20_G0b<-subset(K200_20_G0,subset = K200_20_G0$baby==0.28)

K200_10_G0<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_10/G0_no social learning/baby_0.19to0.28/K200_10_baby_0.19to0.28_G0.csv"))
K200_10_G0b<-subset(K200_10_G0,subset = K200_10_G0$baby==0.28)

K200_40_G0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_40/G0.5_no social learning/baby0.28/K200_40_G0.5.csv"))
K200_20_G0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_20/G0.5_no social learning/baby0.28/K200_20_G0.5.csv"))
K200_10_G0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_10/G0.5_no social learning/baby0.28/K200_10_G0.5.csv"))

K200_40_social0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_40/G0_social learning/baby0.28/K200_40_social0.5.csv"))
K200_20_social0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_20/G0_social learning/baby0.28/K200_20_social0.5.csv"))
K200_10_social0.5<-read_csv(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_10/G0_social learning/baby0.28/K200_10_social0.5.csv"))

K200_40_G0plot<-ggplot(K200_40_G0b,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle("Basic scenario")+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_40_G0.5plot<-ggplot(K200_40_G0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle("Generalization")+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_40_social0.5plot<-ggplot(K200_40_social0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle("Social learning")+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_20_G0plot<-ggplot(K200_20_G0b,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_20_G0.5plot<-ggplot(K200_20_G0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_20_social0.5plot<-ggplot(K200_20_social0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  xlab("")+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_10_G0plot<-ggplot(K200_10_G0b,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_10_G0.5plot<-ggplot(K200_10_G0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

K200_10_social0.5plot<-ggplot(K200_10_social0.5,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  ylab("")+
  geom_point(data = df,mapping = aes(x=P,y=F),size=3,shape=21,color="black",fill="white")

tiff("figure2.tiff",width = 10,height = 10,units = "in",res = 300)
gridExtra::grid.arrange(K200_40_G0plot,K200_40_G0.5plot,K200_40_social0.5plot,
             K200_20_G0plot,K200_20_G0.5plot,K200_20_social0.5plot,
             K200_10_G0plot,K200_10_G0.5plot,K200_10_social0.5plot,
             ncol=3)
dev.off()

 # Mapping empirical learning data onto the P-F space
# for G = 0
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.csv")
# For G = 0.5
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.5.csv")
# For G0_Eff0.5_Spr0.5
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.5.csv")


P_by_F$P<-rep(seq(0,0.5,length.out=10),each=20)
P_by_F$F<-rep(rep(seq(0,0.5,length.out=10),each=2),10)

P_by_F_babylow<-subset(P_by_F,subset = P_by_F$baby==0.19)
P_by_F_babyhigh<-subset(P_by_F,subset = P_by_F$baby==0.27)

df<-read_csv("~/Google Drive/Predator learning and warning signal diversity/P&F.csv")

plot1<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle(paste0("r0 = 0.17"))+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

plot2<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle(paste0("r0 = 0.28"))+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

grid.arrange(plot1,plot2,ncol=2)

# Social learning (annual prey)
library(tidyverse)
library(parallel)

source("learning_functions.R")

vector_P<-seq(0,0.5,length.out = 10)
vector_F<-seq(0,0.5,length.out = 10)

PbyF<-expand.grid(vector_P,vector_F)

v_PbyF<-PbyF %>% slice(rep(1:n(),each = 3))

baby<-rep(c(2,5,10),times=100)

parameter_space<-cbind(v_PbyF,baby)

no.Cores<-detectCores()

# For perennial prey, use mimic_social_fast
# For annual prey, use mimic_social_fast2

for (i in 1:10){
  
  dir<-paste0("D:/Desktop/",i)
  
  dir.create(dir)
  
  setwd(dir)
  
  clust<-makeCluster(no.Cores, setup_timeout = 0.5)
  
  clusterEvalQ(clust,library(tidyverse))
  clusterExport(clust,"learning")
  clusterExport(clust,"mimic_social_fast2")
  clusterExport(clust,"parameter_space")
  
  
  
  parApply(cl = clust,parameter_space,1,function(y){mimic_social_fast2(G = 0,P_OLD = y[1],Abundance = c(200,40),
                                                                       lambda_L = c(y[1],y[1]),alpha_F = y[2],baby = y[3],
                                                                       K = c(200,40),social_spread = 0,social_efficacy = 0,
                                                                       inter_learning = TRUE,
                                                                       GenTime = 100,timestep = 50000)})
  
  stopCluster(clust)
  
}

# Summarizing social learning results
library(tidyverse)
library(gridExtra)

rm(list=ls())

K<-c(200,20)

baby<-c(2,5,10)

for (j in 6:10){
  
  replication<-j
  
  # On Windows
  dir<-paste0("D:/Desktop/K",K[1],"_",K[2],"/",replication)
  
  # On Mac
  # dir<-paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_annualprey/K",K[1],"_",K[2],"/",replication)
  
  setwd(dir)
  
  files<-unique(list.files(dir))
  
  P_by_F<-matrix(0,nrow=length(files),ncol=6)

  
  for (i in 1:length(files)){
    
    # string<-files[i]
    # first<-str_extract(string = string,pattern = "G[0-1]+")
    # G<-as.numeric(substring(first,2:))
    # second<-str_extract(string = string,pattern = "NEW[0-9]+")
    # abnc_new_t0<-as.numeric(substring(second,4))
    
    data<-read_csv(files[i])
    abnc_OLD_final<-data[1,]$abnc_OLD
    abnc_new_final<-data[1,]$abnc_NEW
    
    P_by_F[i,1]<-as.numeric(str_split(str_split(as.character(files[i]),pattern = "_")[[1]][1],"P")[[1]][2])
    P_by_F[i,2]<-as.numeric(str_split(str_split(as.character(files[i]),pattern = "_")[[1]][2],"F")[[1]][2])
    P_by_F[i,3]<-as.numeric(str_split(as.character(files[i]),pattern = "_")[[1]][4])
    
    P_by_F[i,c(4:6)]<-c(abnc_OLD_final,abnc_new_final,-99)
    
  }
  
  P_by_F<-data.frame(P_by_F)
  colnames(P_by_F)<-c("P","F","baby","abnc_OLD_final","abnc_NEW_final","outcome")
  
  for (i in 1:nrow(P_by_F)){
    
    if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-0
    } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-1
    } else {P_by_F[i,]$outcome<-2}
    
  }
  
  P_by_F$selection_NEW<-(P_by_F$abnc_NEW_final/K[2])/(P_by_F$abnc_OLD_final/K[1])-1
  
  # P_by_F[is.na(P_by_F$selection_NEW),]$selection_NEW<--1
  
  # P_by_F_babylow<-subset(P_by_F,subset = P_by_F$baby==baby[1])
  # P_by_F_babyhigh<-subset(P_by_F,subset = P_by_F$baby==baby[2])
  # 
  # plot1<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[1])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot2<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=selection_NEW))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[1])))+
  #   theme_minimal()+  
  #   theme(legend.position = "none")+
  #   scale_fill_distiller(palette = "YlGnBu")
  # 
  # plot3<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=outcome))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[2])))+
  #   theme_minimal()+
  #   theme(legend.position = "none")+
  #   scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)
  # 
  # plot4<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=selection_NEW))+
  #   geom_tile()+
  #   ggtitle(paste0(paste0("r0 = ",baby[2])))+
  #   theme_minimal()+  
  #   theme(legend.position = "none")+
  #   scale_fill_distiller(palette = "YlGnBu")
  # 
  # library(gridExtra)
  # 
  # grid.arrange(plot1,plot3,plot2,plot4,ncol=2)
  
  dir2<-paste0("D:/Desktop/K200_",K[2])
  
  setwd(dir2)
  
  write_csv(P_by_F,paste0(replication,".csv"))
  
}

# Examine each replicate individually

rm(list=ls())

K<-c(200,10)

baby<-c(2,5,10)

setwd(paste0("~/Google Drive/Predator learning and warning signal diversity/social_annualprey/K200_",K[2]))

files<-list.files()

for (i in 1:length(files)){
  
  P_by_F<-read_csv(files[i])
  
  P_by_F_babylow<-P_by_F[P_by_F$baby==baby[1],]
  P_by_F_babymed<-P_by_F[P_by_F$baby==baby[2],]
  P_by_F_babyhigh<-P_by_F[P_by_F$baby==baby[3],]
  
  plot1<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=abnc_NEW_final))+
    geom_tile()+
    ggtitle(paste0(paste0("birth = ",baby[1])))+
    theme_minimal()+
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "Greys")+
    geom_text(aes(x = P,y = F,label=abnc_NEW_final))
  
  plot2<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=selection_NEW))+
    geom_tile()+
    ggtitle(paste0(paste0("r0 = ",baby[1])))+
    theme_minimal()+  
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "YlGnBu")
  
  plot3<-ggplot(P_by_F_babymed,aes(x = P,y = F,fill=abnc_NEW_final))+
    geom_tile()+
    ggtitle(paste0(paste0("birth = ",baby[2])))+
    theme_minimal()+
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "Greys")+
    geom_text(aes(x = P,y = F,label=abnc_NEW_final))
  
  plot4<-ggplot(P_by_F_babymed,aes(x = P,y = F,fill=selection_NEW))+
    geom_tile()+
    ggtitle(paste0(paste0("r0 = ",baby[2])))+
    theme_minimal()+  
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "YlGnBu")
  
  plot5<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=abnc_NEW_final))+
    geom_tile()+
    ggtitle(paste0(paste0("birth = ",baby[3])))+
    theme_minimal()+
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "Greys")+
    geom_text(aes(x = P,y = F,label=abnc_NEW_final))
  
  plot6<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=selection_NEW))+
    geom_tile()+
    ggtitle(paste0(paste0("r0 = ",baby[3])))+
    theme_minimal()+  
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "YlGnBu")
  
  library(gridExtra)
  
  grid.arrange(plot1,plot2,plot3,plot4,plot5,plot6,ncol=2)
  
}

# Summarizing results from different replications
rm(list=ls())

K<-c(200,20)

baby<-c(2,5,10)

##### Change the working directory as appropriate
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_annualprey/K200_",K[2]))
files<-list.files()

outcome<-matrix(-99,nrow = 300,ncol = 10)
selection_NEW<-matrix(-99,nrow = 300,ncol = 10)

for (i in 1:length(files)){
  
  df<-read_csv(files[i])
  
  outcome[,i]<-df$abnc_NEW_final
  selection_NEW[,i]<-df$selection_NEW
  
}

P_by_F<-data.frame(matrix(0,nrow = 300,ncol = 7))
colnames(P_by_F)<-c("P","F","baby","outcome","outcome_sd","selection_NEW","selection_NEW_sd")

P_by_F[,1:3]<-df[,1:3]

for (i in 1:300){
  
  # Change X in outcome[i,1:X] as new replications are done
  P_by_F[i,4]<-mean(outcome[i,1:length(files)])
  P_by_F[i,5]<-sd(outcome[i,1:length(files)])
  P_by_F[i,6]<-mean(selection_NEW[i,1:length(files)])
  P_by_F[i,7]<-sd(selection_NEW[i,1:length(files)])
  
}

# Heatmap of P x F based on outcome and selection against the new prey
# Create two lists to store ggplot objects in the loop
# Run this line to directly load summary datasets
rm(list=ls())

K<-c(200,40)

baby<-c(2,5,10)

##### Change the working directory as appropriate
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_annualprey/K200_",K[2]))

P_by_F<-read_csv("summary.csv")

plot1_list<-list(NA)
plot2_list<-list(NA)

baby<-c(2,5,10)

for (i in 1:length(baby)){
  
  df<-P_by_F[P_by_F$baby==baby[i],]
  
  plot1_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=outcome))+
    geom_tile()+
    ggtitle(paste0("birth = ",baby[i]))+
    theme_minimal()+
    theme(legend.position = "none")+
    scale_fill_distiller(palette = "Greys",limits = c(min(P_by_F$outcome),max(P_by_F$outcome)))+
    geom_text(aes(x = P,y = F,label=outcome))
  
  # Heatmap of P x F based on selection against the NEW prey
  
  plot2_list[[i]]<-ggplot(df,aes(x = P,y = F,fill=selection_NEW))+
    geom_tile()+
    ggtitle(paste0("birth = ",baby[i]))+
    theme_minimal()+  
    theme(legend.position = "right")+
    scale_fill_distiller(palette = "YlGnBu",limits = c(min(P_by_F$selection_NEW),max(P_by_F$selection_NEW)))
  
}

library(gridExtra)

grid.arrange(plot1_list[[1]],plot1_list[[2]],plot1_list[[3]],ncol=3)

grid.arrange(plot2_list[[1]],plot2_list[[2]],plot2_list[[3]],ncol=3)

write_csv(P_by_F,"summary.csv")


# Mapping empirical learning data onto the P-F space

# for G = 0
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.csv")
# For G = 0.5
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.5.csv")
# For G0_Eff0.5_Spr0.5
P_by_F<-read_csv("~/Google Drive/Predator learning and warning signal diversity/social learning/G0.5.csv")


P_by_F$P<-rep(seq(0,0.5,length.out=10),each=20)
P_by_F$F<-rep(rep(seq(0,0.5,length.out=10),each=2),10)

P_by_F_babylow<-subset(P_by_F,subset = P_by_F$baby==0.19)
P_by_F_babyhigh<-subset(P_by_F,subset = P_by_F$baby==0.27)

df<-read_csv("~/Google Drive/Predator learning and warning signal diversity/P&F.csv")

plot1<-ggplot(P_by_F_babylow,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle(paste0("r0 = 0.17"))+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

plot2<-ggplot(P_by_F_babyhigh,aes(x = P,y = F,fill=outcome))+
  geom_tile()+
  ggtitle(paste0("r0 = 0.28"))+
  theme_minimal()+
  theme(legend.position = "none")+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

gridExtra::grid.arrange(plot1,plot2,ncol=2)
