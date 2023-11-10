library(tidyverse)
library(parallel)

# 
# Separate K

source("learning_functions.R")

x2<-c(20,100,180)
K2<-c(20,100,180,200)
vector_P<-seq(0,0.5,length.out = 10)
vector_F<-seq(0,0.5,length.out = 10)
baby<-c(0.2,0.3)

all_combo<-matrix(unlist(expand.grid(x2,K2,vector_P,vector_F,baby)),ncol=5)
parameter_space<-all_combo[all_combo[,2]==all_combo[,1]|all_combo[,2]==200,]

G<-0
social_spread<-0.5
social_efficacy<-0.5

no.Cores<-detectCores()

# For perennial prey, use mimic_social_fast
# For annual prey, use mimic_social_fast2
# For perennial prey that share the same carrying capacity, use mimic_social_fast3

for (i in 1:10){
  
  dir<-paste0("D:/Desktop/separate K/",i)
  
  dir.create(dir)
  
  setwd(dir)
  
  clust<-makeCluster(no.Cores, setup_timeout = 0.5)
  
  clusterEvalQ(clust,library(tidyverse))
  clusterExport(clust,"learning")
  clusterExport(clust,"mimic_social_fast")
  clusterExport(clust,"parameter_space")
  clusterExport(clust,"G")
  clusterExport(clust,"social_spread")
  clusterExport(clust,"social_efficacy")
  
  
  
  parApply(cl = clust,parameter_space,1,function(y){mimic_social_fast(G = G,P_OLD = 0.5,Abundance = c(200,y[1]),
                                                                 lambda_L = c(y[3],y[3]),alpha_F = y[4],baby = y[5],
                                                                 K = c(200,y[2]),social_spread = social_spread,social_efficacy = social_efficacy,
                                                                 inter_learning = TRUE,
                                                                 GenTime = 100,timestep = 50000)})
  
  stopCluster(clust)
  
}

# perennial prey with shared carrying capacity
library(tidyverse)
library(parallel)

source("learning_functions.R")

x2<-c(20,100,180)
K<-c(220,300,380,400)
vector_P<-seq(0,0.5,length.out = 10)
vector_F<-seq(0,0.5,length.out = 10)
baby<-c(0.2,0.3)

all_combo<-matrix(unlist(expand.grid(x2,K,vector_P,vector_F,baby)),ncol=5)
parameter_space<-all_combo[all_combo[,2]-200==all_combo[,1]|all_combo[,2]-200==200,]

G<-0
social_spread<-0.5
social_efficacy<-0.5

no.Cores<-detectCores()

for (i in 6:10){
  
  dir<-paste0("D:/Desktop/shared K/",i)
  
  dir.create(dir)
  
  setwd(dir)
  
  clust<-makeCluster(no.Cores, setup_timeout = 0.5)
  
  clusterEvalQ(clust,library(tidyverse))
  clusterExport(clust,"learning")
  clusterExport(clust,"mimic_social_fast3")
  clusterExport(clust,"parameter_space")
  clusterExport(clust,"G")
  clusterExport(clust,"social_spread")
  clusterExport(clust,"social_efficacy")
  
  
  
  parApply(cl = clust,parameter_space,1,function(y){mimic_social_fast3(G = G,P_OLD = 0.5,Abundance = c(200,y[1]),
                                                                       lambda_L = c(y[3],y[3]),alpha_F = y[4],baby = y[5],
                                                                       K =y[2],social_spread = social_spread,social_efficacy = social_efficacy,
                                                                       inter_learning = TRUE,
                                                                       GenTime = 100,timestep = 50000)})
  
  stopCluster(clust)
  
}

## in case simulations are interrupted:
files<-list.files("D:/Desktop/shared K/5_incomplete")

combos<-NULL

for (i in 1:nrow(parameter_space)){
  
  combos[i]<-(paste0("rare",parameter_space[i,1],"_K",parameter_space[i,2],"_P",round(parameter_space[i,3],3),"_F",round(parameter_space[i,4],3),"_baby_",parameter_space[i,5],"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
  
  
}

missing<-setdiff(combos,files)

parameter_space_missing<-matrix(999,nrow = length(missing),ncol = 5)

for (i in 1:length(missing)){
  
  parameter_space_missing[i,1]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
  parameter_space_missing[i,2]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][2],pattern = "K")[[1]][2])
  parameter_space_missing[i,3]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][3],pattern = "P")[[1]][2]) 
  parameter_space_missing[i,4]<-as.numeric(str_split(str_split(missing[i],pattern = "_")[[1]][4],pattern = "F")[[1]][2])
  parameter_space_missing[i,5]<-as.numeric((str_split(missing[i],pattern = "_")[[1]][6]))
}

setwd("D:/Desktop/shared K/")

clust<-makeCluster(no.Cores, setup_timeout = 0.5)

clusterEvalQ(clust,library(tidyverse))
clusterExport(clust,"learning")
clusterExport(clust,"mimic_social_fast3")
clusterExport(clust,"parameter_space_missing")
clusterExport(clust,"G")
clusterExport(clust,"social_spread")
clusterExport(clust,"social_efficacy")


parApply(cl = clust,parameter_space_missing,1,function(y){mimic_social_fast3(G = G,P_OLD = 0.5,Abundance = c(200,y[1]),
                                                                             lambda_L = c(y[3],y[3]),alpha_F = y[4],baby = y[5],
                                                                             K =y[2],social_spread = social_spread,social_efficacy = social_efficacy,
                                                                             inter_learning = TRUE,
                                                                             GenTime = 100,timestep = 50000)})

stopCluster(clust)
## 

# Summarizing results from each replication
# separate K
library(tidyverse)
library(gridExtra)

rm(list=ls())

G<-0

for (j in 1:10){
  
  replication<-j
  
  # # G = 0
  # dir<-paste0("D:/Desktop/separate K/G0/",replication)
  # # G = 0.5
  # dir<-paste0("D:/Desktop/separate K/G0.5/",replication)
  # social learning
  dir<-paste0("D:/Desktop/separate K/social/",replication)
  
  setwd(dir)
  
  files<-unique(list.files(dir))
  
  P_by_F<-matrix(0,nrow=length(files),ncol=10)
  P_by_F<-data.frame(P_by_F)
  colnames(P_by_F)<-c("r","x2_0","K2","lambda","alphaF","G","abnc_OLD_final","abnc_NEW_final","outcome","selection_NEW")
  
  
  for (i in 1:length(files)){
    
    data<-read_csv(files[i])
    
    P_by_F[i,1]<-as.numeric(str_split(files[i],pattern = "_")[[1]][6])
    P_by_F[i,2]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
    P_by_F[i,3]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][2],pattern = "K2")[[1]][2])
    P_by_F[i,4]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][3],pattern = "P")[[1]][2])
    P_by_F[i,5]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][4],pattern = "F")[[1]][2])
    P_by_F[i,6]<-G
    
    
    abnc_OLD_final<-data[1,]$abnc_OLD
    abnc_new_final<-data[1,]$abnc_NEW
    
    P_by_F[i,c(7:10)]<-c(abnc_OLD_final,abnc_new_final,-99,-99)
    
  }
  
  
  for (i in 1:nrow(P_by_F)){
    
    if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-0
    } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-1
    } else {P_by_F[i,]$outcome<-2}
    
  }
  
  # dir2<-paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_no social learning/baby0.28/")
  dir2<-paste0("G:/我的雲端硬碟/Predator learning and warning signal diversity/social_perennialprey")
  
  dir.create(dir2)
  
  setwd(dir2)
  
  write_csv(P_by_F,paste0(replication,".csv"))
  
}

# Summarizing results from different replications
# separate K
rm(list=ls())

G<-0.5

##### Change the working directory as appropriate
# G = 0
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/G0"))
# G = 0.5
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/G0.5"))
# social learning
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/social"))


files<-list.files()

df<-read_csv(files[1])

for (i in 2:length(files)){
  
  dff<-read_csv(files[i])
  df<-rbind(df,dff)
  
}

df$id<-paste0(df$r,df$x2_0,df$K2,df$lambda,df$alphaF)

summary<-data.frame(matrix(0,nrow = 1200,ncol = 10))
colnames(summary)<-colnames(df)[1:10]

summary[,c(1:8)]<-df[1:1200,c(1:8)]

id_list<-unique(df$id)

for (i in 1:length(id_list)){
  
  summary[i,9]<-mean(df[df$id==id_list[i],]$outcome)
  summary[i,10]<-mean(df[df$id==id_list[i],]$selection_NEW)
  
  
}

# G = 0
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_G0.csv")

# G = 0.5
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_G0.5.csv")

# Social learning
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_social.csv")

# Heatmap of P x F based on outcome and selection against the new prey

# Run this line when directly reading the summary data frame
summary_G0<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_G0.csv"))
summary_G0.5<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_G0.5.csv"))
summary_social<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/summary_social.csv"))

summary_social$G<-99

df_plot<-rbind(summary_G0,summary_G0.5,summary_social)

df_plot$scenario<-paste0("x2 = ",df_plot$x2_0,", K2 = ",df_plot$K2)

df_plot$priority<-NA

for (i in 1:nrow(df_plot)){
  
  if (df_plot[i,]$scenario=="x2 = 20, K2 = 20"){df_plot[i,]$priority<-"A"
  } else if (df_plot[i,]$scenario=="x2 = 20, K2 = 200"){df_plot[i,]$priority<-"B"
  } else if (df_plot[i,]$scenario=="x2 = 100, K2 = 100"){df_plot[i,]$priority<-"C"
  } else if (df_plot[i,]$scenario=="x2 = 100, K2 = 200"){df_plot[i,]$priority<-"D"
  } else if (df_plot[i,]$scenario=="x2 = 180, K2 = 180"){df_plot[i,]$priority<-"E"
  } else if (df_plot[i,]$scenario=="x2 = 180, K2 = 200"){df_plot[i,]$priority<-"F"}
  
}

scenario_list<-c(A = "x2 = 20, K2 = 20",
                 B = "x2 = 20, K2 = 200",
                 C = "x2 = 100, K2 = 100",
                 D = "x2 = 100, K2 = 200",
                 E = "x2 = 180, K2 = 180",
                 F = "x2 = 180, K2 = 200")

df_plot$selection_NEW<-(df_plot$abnc_NEW_final/df_plot$x2_0)/(df_plot$abnc_OLD_final/200)-1
df_plot[is.na(df_plot$selection_NEW)==TRUE,]$selection_NEW<-0
df_plot[df_plot$selection_NEW==Inf,]$selection_NEW<-0

write_csv(df_plot,"sepK_summary_imb.csv")

# outcome
  sepK_outcome_G0_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==0 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_outcome_G0_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_outcome_G0.5_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==0.5 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_outcome_G0.5_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_outcome_social_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==99 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_outcome_social_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
# selection
  sepK_sel_G0_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.3)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_sel_G0_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.5)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_sel_G0.5_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.3)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_sel_G0.5_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.3)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

  sepK_sel_social_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.3)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  sepK_sel_social_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
    geom_raster()+
    # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
    theme_bw()+
    mdthemes::md_theme_bw()+
    theme(legend.position = "right")+
    labs(x = "\u019b",y = paste0("\u0251","F"))+
    scale_fill_distiller(type = "div",limits=c(-1.5,1.5),breaks=c(seq(-1.5,1.5,by = 0.3)),name = "S")+
    facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)
  
  
  gridExtra::grid.arrange(sepK_outcome_G0_r0.2,sepK_outcome_G0_r0.3,sepK_outcome_G0.5_r0.2,
                          sepK_outcome_G0.5_r0.3,sepK_outcome_social_r0.2,sepK_outcome_social_r0.3,ncol=3,as.table=FALSE)
  
  gridExtra::grid.arrange(sepK_sel_G0_r0.2,sepK_sel_G0_r0.3,sepK_sel_G0.5_r0.2,sepK_sel_G0.5_r0.3,
                          sepK_sel_social_r0.2,sepK_sel_social_r0.3,ncol=3,as.table=FALSE)

sepK_plot<-read_csv("sepK_summary_imb.csv")
  
hist_S_sepK_imb<-ggplot(data = sepK_plot[df_plot$outcome==2,])+
    geom_histogram(aes(x = selection_NEW))+
    mdthemes::md_theme_bw()+
    labs(x = "Selection on rare signal (S)", y = "Count")+
    ggtitle("(C) Separate carrying capacities")
  
# Mapping empirical learning data onto the P-F space
df_plot<-read_csv("sepK_summary_imb.csv")

df<-read_csv("P&F.csv")

# r = 0.3, G = 0
sepK_plot1<-ggplot(subset(df_plot,subset = r==0.3 & G == 0 & x2_0 == 20 & K2 == 20),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

# G = 0.5
sepK_plot2<-ggplot(subset(df_plot,subset = r==0.3 & G == 0.5 & x2_0 == 20 & K2 == 20),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

# Social learning
sepK_plot3<-ggplot(subset(df_plot,subset = r==0.3 & G == 99 & x2_0 == 20 & K2 == 20),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

gridExtra::grid.arrange(sepK_plo1,sepK_plot2,sepK_plot3,ncol=3)


# Summarizing results from each replication
# shared K
library(tidyverse)
library(gridExtra)

rm(list=ls())

G<-0

for (j in 1:10){
  
  replication<-j
  
  # # G = 0
  # dir<-paste0("D:/Desktop/shared K/G0/",replication)
  # # G = 0.5
  # dir<-paste0("D:/Desktop/shared K/G0.5/",replication)
  # social learning
  dir<-paste0("D:/Desktop/shared K/social/",replication)
  
  setwd(dir)
  
  files<-unique(list.files(dir))
  
  P_by_F<-matrix(0,nrow=length(files),ncol=10)
  P_by_F<-data.frame(P_by_F)
  colnames(P_by_F)<-c("r","x2_0","K","lambda","alphaF","G","abnc_OLD_final","abnc_NEW_final","outcome","selection_NEW")
  
  
  for (i in 1:length(files)){
    
    data<-read_csv(files[i])
    
    P_by_F[i,1]<-as.numeric(str_split(files[i],pattern = "_")[[1]][6])
    P_by_F[i,2]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][1],pattern = "rare")[[1]][2])
    P_by_F[i,3]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][2],pattern = "K")[[1]][2])
    P_by_F[i,4]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][3],pattern = "P")[[1]][2])
    P_by_F[i,5]<-as.numeric(str_split(str_split(files[i],pattern = "_")[[1]][4],pattern = "F")[[1]][2])
    P_by_F[i,6]<-G
    
    
    abnc_OLD_final<-data[1,]$abnc_OLD
    abnc_new_final<-data[1,]$abnc_NEW
    
    P_by_F[i,c(7:10)]<-c(abnc_OLD_final,abnc_new_final,-99,-99)
    
  }
  
  
  for (i in 1:nrow(P_by_F)){
    
    if (P_by_F[i,]$abnc_OLD_final==0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-0
    } else if (P_by_F[i,]$abnc_OLD_final!=0&P_by_F[i,]$abnc_NEW_final==0){
      P_by_F[i,]$outcome<-1
    } else {P_by_F[i,]$outcome<-2}
    
  }
  
  # dir2<-paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/Summary_K200_",K[2],"/G",G,"_no social learning/baby0.28/")
  dir2<-paste0("G:/我的雲端硬碟/Predator learning and warning signal diversity/social_perennialprey/shared K")
  
  dir.create(dir2)
  
  setwd(dir2)
  
  write_csv(P_by_F,paste0(replication,".csv"))
  
}

# Summarizing results from all replications
# separate K
rm(list=ls())

G<-0

##### Change the working directory as appropriate
# G = 0
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/G0"))
# G = 0.5
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/G0.5"))
# social learning
setwd(paste0("~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/social"))


files<-list.files()

df<-read_csv(files[1])

for (i in 2:length(files)){
  
  dff<-read_csv(files[i])
  df<-rbind(df,dff)
  
}

df$id<-paste0(df$r,df$x2_0,df$K,df$lambda,df$alphaF)

summary<-data.frame(matrix(0,nrow = 1200,ncol = 10))
colnames(summary)<-colnames(df)[1:10]

summary[,c(1:8)]<-df[1:1200,c(1:8)]

id_list<-unique(df$id)

for (i in 1:length(id_list)){
  
  summary[i,9]<-mean(df[df$id==id_list[i],]$outcome)
  
}

summary$selection_NEW<-(summary$abnc_NEW_final/summary$x2_0)/(summary$abnc_OLD_final/200)-1

# G = 0
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_G0.csv")

# G = 0.5
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_G0.5.csv")

# Social learning
write_csv(summary,"~/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_social.csv")

# Heatmap of P x F based on outcome and selection against the new prey

# Run this line when directly reading the summary data frame
summary_G0<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_G0.csv"))
summary_G0.5<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_G0.5.csv"))
summary_social<-read_csv(paste0("/Users/chi/Library/CloudStorage/GoogleDrive-kuochiyun@gmail.com/My Drive/Predator learning and warning signal diversity/social_perennialprey/shared K/summary_social.csv"))

summary_social$G<-99

df_plot<-rbind(summary_G0,summary_G0.5,summary_social)

df_plot$scenario<-paste0("x2 = ",df_plot$x2_0,", K = ",df_plot$K)

df_plot$priority<-NA

for (i in 1:nrow(df_plot)){
  
  if (df_plot[i,]$scenario=="x2 = 20, K = 220"){df_plot[i,]$priority<-"A"
  } else if (df_plot[i,]$scenario=="x2 = 20, K = 400"){df_plot[i,]$priority<-"B"
  } else if (df_plot[i,]$scenario=="x2 = 100, K = 300"){df_plot[i,]$priority<-"C"
  } else if (df_plot[i,]$scenario=="x2 = 100, K = 400"){df_plot[i,]$priority<-"D"
  } else if (df_plot[i,]$scenario=="x2 = 180, K = 380"){df_plot[i,]$priority<-"E"
  } else if (df_plot[i,]$scenario=="x2 = 180, K = 400"){df_plot[i,]$priority<-"F"}
  
}

scenario_list<-c(A = "x2 = 20, K = 20",
                 B = "x2 = 20, K = 220",
                 C = "x2 = 100, K = 300",
                 D = "x2 = 100, K = 400",
                 E = "x2 = 180, K = 380",
                 F = "x2 = 180, K = 400")

# df_plot$selection_NEW<-(df_plot$abnc_NEW_final/df_plot$x2_0)/(df_plot$abnc_OLD_final/200)-1
df_plot[is.na(df_plot$selection_NEW)==TRUE,]$selection_NEW<-0
df_plot[df_plot$selection_NEW==Inf,]$selection_NEW<-600

write_csv(df_plot,"sharedK_summary_imb.csv")

# outcome
sharedK_outcome_G0_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==0 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_G0_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_G0.5_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==0.5 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_G0.5_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_social_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G ==99 &priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_outcome_social_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = 0.2,", "G = ",G, ", H = 0.6"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

# selection
sharedK_sel_G0_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_G0_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==0 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_G0.5_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_G0.5_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==0.5 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_social_r0.2<-ggplot(subset(df_plot,subset = r==0.2 & G==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)

sharedK_sel_social_r0.3<-ggplot(subset(df_plot,subset = r==0.3 & G ==99 & priority %in% c("A","C","E")),aes(x = lambda,y = alphaF,fill=round(selection_NEW,3)))+
  geom_raster()+
  # ggtitle(paste0("Separate K, r0 = ",r0,", G = ",G, ", H = ",H))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "right")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_distiller(type = "div",limits=c(-10,10),breaks=c(seq(-10,10,by = 2)),name = "S")+
  facet_wrap(~priority,labeller = labeller(.cols = scenario_list),ncol=3)


gridExtra::grid.arrange(sharedK_outcome_G0_r0.2,sharedK_outcome_G0_r0.3,sharedK_outcome_G0.5_r0.2,
                        sharedK_outcome_G0.5_r0.3,sharedK_outcome_social_r0.2,sharedK_outcome_social_r0.3,ncol=3,as.table=FALSE)

gridExtra::grid.arrange(sharedK_sel_G0_r0.2,sharedK_sel_G0_r0.3,sharedK_sel_G0.5_r0.2,sharedK_sel_G0.5_r0.3,
                        sharedK_sel_social_r0.2,sharedK_sel_social_r0.3,ncol=3,as.table=FALSE)

sharedK_plot<-read_csv("sharedK_summary_imb.csv")

hist_S_sharedK_imb<-ggplot(data = sharedK_plot[df_plot$outcome==2,])+
  geom_histogram(aes(x = selection_NEW))+
  mdthemes::md_theme_bw()+
  labs(x = "Selection on rare signal (S)", y = "Count")+
  ggtitle("(D) Shared carrying capacities")

gridExtra::grid.arrange(hist_S_sepK_imb,hist_S_sharedK_imb,ncol=2)


# Mapping empirical learning data onto the P-F space
df_plot<-read_csv("sharedK_summary_imb.csv")

df<-read_csv("P&F.csv")

# r = 0.3, G = 0
sharedK_plot1<-ggplot(subset(df_plot,subset = r==0.3 & G == 0 & x2_0 == 20 & K == 220),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

# G = 0.5
sharedK_plot2<-ggplot(subset(df_plot,subset = r==0.3 & G == 0.5 & x2_0 == 20 & K == 220),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

# Social learning
sharedK_plot3<-ggplot(subset(df_plot,subset = r==0.3 & G == 99 & x2_0 == 20 & K == 220),aes(x = lambda,y = alphaF,fill=outcome))+
  geom_raster()+
  # ggtitle(paste0("r0 = 0.17"))+
  theme_bw()+
  mdthemes::md_theme_bw()+
  theme(legend.position = "none")+
  labs(x = "\u019b",y = paste0("\u0251","F"))+
  scale_fill_gradient2(low = "black",mid = "red",high = "blue",midpoint = 1)+
  geom_point(data = df,mapping = aes(x=P,y=F,size=3),shape=21,color="black",fill="white")

gridExtra::grid.arrange(sharedK_plot1,sharedK_plot2,sharedK_plot3,ncol=3)
