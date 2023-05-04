# Some proof-of-concept simulations for attack prob estimation

source("learning_functions.R")

# In an experiment, where 
mimic(G = 0,P_OLD = 0.4,Abundance = c(1000,1000),lambda_L = c(0.4,0.1),alpha_F = c(0,0),baby = 0,K = c(1,1000),GenTime = 100,timestep = 1000)

df<-read_csv("P0.4_F0.csv")

nrow(df[df$outcome=="OLD",])
nrow(df[df$outcome=="NEW",])


asymp_atk1<-0.4
asymp_atk2<-0.2

asymp_atk_estimated<-NULL

for (k in 1:1000){

outcome1<-NULL
outcome2<-NULL
outcome1[1]<-outcome2[1]<-"A"

i<-1
j<-1

while(i<50){

  outcome1[j]<-sample(x = c("A","P"),size = 1,replace = TRUE,prob = c(asymp_atk1,1-asymp_atk1))
  outcome2[j]<-sample(x = c("A","P"),size = 1,replace = TRUE,prob = c(asymp_atk2,1-asymp_atk2))
  
  i<-length(outcome1[outcome1=="A"])+length(outcome2[outcome2=="A"])
  j<-j+1
}

asymp_atk_estimated[k]<-length(outcome2[outcome2=="A"])/length(outcome1[outcome1=="A"])
}

hist(asymp_atk_estimated)

# Meta-analysis of the determinants of attack prob
library(tidyverse)
library(lme4)

source("HighstatLibV10.R")

df<-read_csv("attack_prob.csv")
df$study<-paste0(df$author,df$journal,df$volume)

# 1. Use data with us_concentration_cal (i.e., direct learning using bitters as US)
df2<-subset(df,subset = is.na(df$us_concentration_cal)==FALSE)

lm_all<-lm(attack_prob~us_concentration_cal*cs2*prey_presentation*species,data = df2)
# minus us_concentration
lm_red1<-lm(attack_prob~cs2*prey_presentation*species,data = df2)
# minus cs2
lm_red2<-lm(attack_prob~us_concentration_cal*prey_presentation*species,data = df2)
# minus prey_presentation
lm_red3<-lm(attack_prob~us_concentration_cal*cs2*species,data = df2)
# minus species
lm_red4<-lm(attack_prob~us_concentration_cal*cs2*prey_presentation,data = df2)

AIC(lm_red4,lm_all)
# lm_all is the best per AIC

anova(lm_all)

library(effectsize)
library(parameters)

# Check deviation from normality
op<-par(mfrow=c(2,3),mar = c(5,4,1,2))

plot(lm_all,add.smooth=FALSE,which=1)

E<-resid(lm_all)

hist(E,xlab="Residuals",main="")
plot(df2$us_concentration,E)
abline(h = 0,lty=2,col="gray")
plot(as.factor(df2$cs2),E)
plot(as.factor(df2$prey_presentation),E)
plot(as.factor(df2$species),E)

# Check colinearity among predictors
## us_concentration_cal vs. other predictors
par(mfrow=c(1,3))

plot(x = as.factor(df2$species),y = df2$us_concentration_cal)
plot(x = as.factor(df2$cs2),y = df2$us_concentration_cal)
plot(x = as.factor(df2$prey_presentation),y = df2$us_concentration_cal)

par(mfrow=c(1,3))

## Among ordinal predictors

plot(x = as.factor(df2$species),y = as.factor(df2$prey_presentation),width=2)
plot(x = as.factor(df2$species),y = as.factor(df2$cs2),col=c("red","navyblue","forestgreen","orange","yellow","gray"))
plot(x = as.factor(df2$prey_presentation),y = as.factor(df2$cs2),col=c("red","navyblue","forestgreen","orange","yellow","gray"))

# Effect sizes
lm_all_effectsize<-lm(attack_prob~us_concentration_cal+cs2+prey_presentation+species,data = df2)

standardize_parameters(lm_all_effectsize,method = "refit",two_sd = TRUE)
params<-parameters::model_parameters(lm_all_effectsize) 
effectsize<-t_to_r(params$t[-1],df_error = params$df_error[-1])
interpret_r(effectsize,rules = "gignac2016")

plot1<-ggplot(data=df2,aes(x = us_concentration_cal,y = attack_prob))+
  geom_point(aes(colour=factor(cs2)))+
  # geom_smooth(method = "lm",formula = y~x)+
  theme(legend.position = "top")


plot2<-ggplot(data=df2,aes(x = us_concentration_cal,y = attack_prob))+
  geom_point(aes(colour=factor(prey_presentation)))+
  # geom_smooth(method = "lm",formula = y~x)+
  theme(legend.position = "top")  


plot3<-ggplot(data=df2,aes(x = us_concentration_cal,y = attack_prob))+
  geom_point(aes(colour=factor(species)))+
  # geom_smooth(method = "lm",formula = y~x)+
  theme(legend.position = "top")


gridExtra::grid.arrange(plot1,plot2,plot3,ncol=3)

plot(x = as.factor(df2$cs2),y = df2$us_concentration_cal)
plot(x = as.factor(df2$prey_presentation),y = df2$us_concentration_cal)
plot(x = as.factor(df2$species),y = df2$us_concentration_cal)

# us_concentration main effect
plot(x = df2$us_concentration_cal,y = df2$attack_prob,pch=19)
lm.out<-lm(attack_prob~us_concentration_cal,data = df2)
newx<-seq(min(df2$us_concentration_cal),max(df2$us_concentration_cal),length.out=100)
conf_interval<-predict(lm.out,newdata = data.frame(us_concentration_cal = newx),interval = "confidence",level = 0.95)
abline(lm.out)
matlines(newx,conf_interval[,2:3],col="black",lty=2)


# plot(as.factor(df2$species),df2$attack_prob)
# means<-aggregate(attack_prob~species,data = df2,FUN = mean)
# points(means[,2],col="red")

# cs2 main effect
plot(as.factor(df2$cs2),df2$attack_prob)
means<-aggregate(attack_prob~cs2,data = df2,FUN = mean)
points(means[,2],col="red")

# prey presentation main effect
plot(as.factor(df2$prey_presentation),df2$attack_prob)
means<-aggregate(attack_prob~prey_presentation,data = df2,FUN = mean)
points(means[,2],col="red")

# us_concentration_cal*prey_presentation interaction
p_seq<-subset(df2,subset = df2$prey_presentation=="sequential")
p_seq<-p_seq[p_seq$cs2!="size",]
p_sim<-subset(df2,subset = df2$prey_presentation=="simultaneous")

plot(x = p_sim$us_concentration_cal,y = p_sim$attack_prob)
abline(lm(p_sim$attack_prob~p_sim$us_concentration_cal),lty=2)
points(x = p_seq$us_concentration_cal,y = p_seq$attack_prob,pch=19)
abline(lm(p_seq$attack_prob~p_seq$us_concentration_cal))

# cs2*prey_presentation interaction
par(mfrow=c(1,2))
plot(x = as.factor(p_sim$cs2),y = p_sim$attack_prob)
means<-aggregate(attack_prob~cs2,data = p_sim,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(p_seq$cs2),y = p_seq$attack_prob)
means<-aggregate(attack_prob~cs2,data = p_seq,FUN = mean)
lines(means[,2],type = "b",col="red")

# cs2*species interaction
species<-levels(as.factor(df2$species))

Cc<-df2[df2$species==species[1],]
Ggd<-df2[df2$species==species[2],]
Pm<-df2[df2$species==species[3],]
Sv<-df2[df2$species==species[4],]

par(mfrow=c(2,2))

plot(x = as.factor(Cc$cs2),y = Cc$attack_prob)
means<-aggregate(attack_prob~cs2,data = Cc,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Ggd$cs2),y = Ggd$attack_prob)
means<-aggregate(attack_prob~cs2,data = Ggd,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Pm$cs2),y = Pm$attack_prob)
means<-aggregate(attack_prob~cs2,data = Pm,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Sv$cs2),y = Sv$attack_prob)
means<-aggregate(attack_prob~cs2,data = Sv,FUN = mean)
lines(means[,2],type = "b",col="red")

# 2. Use all data (us_concentration_cal and species are no longer relevant in this analysis)

lm2_all<-lm(attack_prob~cs2*us2*prey_presentation*way_of_learning,data = df)

# minus cs2
lm2_red1<-lm(attack_prob~us2*prey_presentation*way_of_learning,data = df)

# minus us2
lm2_red2<-lm(attack_prob~cs2*prey_presentation*way_of_learning,data = df)

# minus prey_presentation
lm2_red3<-lm(attack_prob~cs2*us2*way_of_learning,data = df)

# minus way_of_learning
lm2_red4<-lm(attack_prob~cs2*us2*prey_presentation,data = df)

AIC(lm2_red4,lm2_red2)
# lm2_red2 is the best model per AIC, but let's use lm2_all anyway

# Check for normality
op<-par(mfrow=c(2,3),mar = c(5,4,1,2))

plot(lm2_all,add.smooth=FALSE,which=1)

E<-resid(lm2_all)

hist(E,xlab="Residuals",main="")
plot(as.factor(df$cs2),E)
plot(as.factor(df$us2),E)
plot(as.factor(df$prey_presentation),E)
plot(as.factor(df$way_of_learning),E)

# Check for collineariy

par(mfrow=c(2,3))

plot(x = as.factor(df$us2),y = as.factor(df$cs2))
plot(x = as.factor(df$us2),y = as.factor(df$prey_presentation))
plot(x = as.factor(df$us2),y = as.factor(df$way_of_learning))
plot(x = as.factor(df$cs2),y = as.factor(df$prey_presentation))
plot(x = as.factor(df$cs2),y = as.factor(df$way_of_learning))
plot(x = as.factor(df$prey_presentation),y = as.factor(df$way_of_learning))

# Effect size
lm2_all_effectsize<-lm(attack_prob~us2+cs2+prey_presentation+way_of_learning,data = df)

standardize_parameters(lm2_all_effectsize,method = "refit",two_sd = TRUE)

params<-parameters::model_parameters(lm2_all_effectsize) 

t_to_r(params$t[-1],df_error = params$df_error[-1])

par(mfrow=c(1,2),mar = c(5,4,1,2))

# cs2 main effect
plot(as.factor(df$cs2),df$attack_prob)
means<-aggregate(attack_prob~cs2,data = df,FUN = mean)
points(means[,2],col="red")

# us2 main effect
plot(as.factor(df$us2),df$attack_prob)
means<-aggregate(attack_prob~us2,data = df,FUN = mean)
points(means[,2],col="red")

# prey_presentation main effect
plot(as.factor(df$prey_presentation),df$attack_prob)
means<-aggregate(attack_prob~prey_presentation,data = df,FUN = mean)
points(means[,2],col="red")

# way_of_learning main effect
plot(as.factor(df$way_of_learning),df$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = df,FUN = mean)
points(means[,2],col="red")

# cs2*prey_presentation interaction
sim<-subset(df,subset = df$prey_presentation=="simultaneous")
seq<-subset(df,subset = df$prey_presentation=="sequential")

par(mfrow=c(1,2))

plot(x = as.factor(sim$cs2),y = sim$attack_prob)
means<-aggregate(attack_prob~cs2,data = sim,FUN = mean)
points(means,type="b",col="red")

plot(x = as.factor(seq$cs2),y = seq$attack_prob)
means<-aggregate(attack_prob~cs2,data = seq,FUN = mean)
points(means,type="b",col="red")

# 3. Use data that use bitters as US but exclude studies that use size and taste as CS
df3<-subset(df2,subset = !df2$cs2 %in% c("size","taste"))
df3b<-subset(df3,subset = df3$attack_prob>0)

# Stats
lm3<-gls(attack_prob~prey_presentation*us_concentration_cal,data = df3)
lm3.1<-gls(attack_prob~prey_presentation*us_concentration_cal,data = df3,weights = varIdent(form = ~1|prey_presentation))
lm3.2<-gls(attack_prob~us_concentration_cal,data = df3)
lm3.3<-lme(attack_prob~species*prey_presentation*us_concentration_cal,random = ~1|study,data = df3)
# Error in MEEM(object, conLin, control$niterEM) : 
#   Singularity in backsolve at level 0, block 1
lm3.4<-gls(attack_prob~species*prey_presentation*us_concentration_cal,data = df3)
# Error in glsEstimate(glsSt, control = glsEstControl) : 
#   computed "gls" fit is singular, rank 15
lm3.5<-glmmTMB::glmmTMB(attack_prob~prey_presentation*us_concentration_cal,family = Gamma(link="inverse"),data = df3b)

AIC(lm3,lm3.1,lm3.2,lm3.5)
# AIC showed that lm3.5 was the best model

# Model validation
library(DHARMa)

plot(E<-simulateResiduals(lm3.5))

par(mfrow=c(1,2))
plotResiduals(E,as.factor(df3b$us_concentration_cal),rank = FALSE)
plotResiduals(E,as.factor(df3b$prey_presentation),rank = FALSE)

#summary
summary(lm3.5)

par(mfrow=c(1,4))

# Graphic presentation of main effects
# prey presentation main effect
plot(as.factor(df3$prey_presentation),df3$attack_prob)
means<-aggregate(attack_prob~prey_presentation,data = df3,FUN = mean)
points(means[,2],col="red")


## prey_presentation:us_concentration_cal interaction
ggplot(df3b,aes(x = us_concentration_cal,y = attack_prob,color = prey_presentation))+
  scale_color_manual(values = c("orange","cornflowerblue"))+
  theme_bw()+
  geom_point(aes(color = prey_presentation),size=3,alpha = 0.5)+
  geom_smooth(method = "lm",fill = NA)+
  theme(axis.title = element_text(size=18),
        axis.text.x = element_text(size=16),
        axis.text.y = element_text(size=16),
        plot.title = element_text(size=18),
        legend.position = "")+
  xlab("Calibrated bitterness (weight %)")+
  ylab("Final attack probability")

# ## species main effect
# plot(as.factor(df3$species),df3$attack_prob)
# means<-aggregate(attack_prob~species,data = df3,FUN = mean)
# points(means[,2],col="red")

# # cs2 main effect
# plot(as.factor(df3$cs2),df3$attack_prob)
# means<-aggregate(attack_prob~cs2,data = df3,FUN = mean)
# points(means[,2],col="red")

# # cs2*prey_presentation interaction
# p_sim<-subset(df3,subset = df3$prey_presentation=="simultaneous")
# p_seq<-subset(df3,subset = df3$prey_presentation=="sequential")
# 
# par(mfrow=c(1,2))
# 
# plot(x = as.factor(p_seq$cs2),y = p_seq$attack_prob)
# means<-aggregate(attack_prob~cs2,data = p_seq,FUN = mean)
# lines(means[,2],type = "b",col="red")
# 
# plot(x = as.factor(p_sim$cs2),y = p_sim$attack_prob)
# means<-aggregate(attack_prob~cs2,data = p_sim,FUN = mean)
# lines(means[,2],type = "b",col="red")

# # cs2*species interaction
# species<-levels(as.factor(df3$species))
# 
# Cc<-df3[df3$species==species[1],]
# Ggd<-df3[df3$species==species[2],]
# Pm<-df3[df3$species==species[3],]
# Sv<-df3[df3$species==species[4],]
# 
# par(mfrow=c(2,2))
# 
# plot(x = as.factor(Cc$cs2),y = Cc$attack_prob)
# means<-aggregate(attack_prob~cs2,data = Cc,FUN = mean)
# lines(means[,2],type = "b",col="red")
# 
# plot(x = as.factor(Ggd$cs2),y = Ggd$attack_prob)
# means<-aggregate(attack_prob~cs2,data = Ggd,FUN = mean)
# lines(means[,2],type = "b",col="red")
# 
# plot(x = as.factor(Pm$cs2),y = Pm$attack_prob)
# means<-aggregate(attack_prob~cs2,data = Pm,FUN = mean)
# lines(means[,2],type = "b",col="red")
# 
# plot(x = as.factor(Sv$cs2),y = Sv$attack_prob)
# means<-aggregate(attack_prob~cs2,data = Sv,FUN = mean)
# lines(means[,2],type = "b",col="red")

# 4. Use all data but exclude studies that use color+pattern+gregariousness, size, taste as CS,
# as well as single-choice experiments

df4<-subset(df,subset = !df$cs2 %in% c("size","taste","color+pattern+gregariousness")&df$prey_presentation!="single choice")

lm4_all<-lm(attack_prob~species*us2*prey_presentation*cs2*way_of_learning,data = df4)

# Check for normality
op<-par(mfrow=c(2,4),mar = c(5,4,1,2))

plot(lm4_all,add.smooth=FALSE,which=1)

E<-resid(lm4_all)

hist(E,xlab="Residuals",main="")
plot(as.factor(df4$species),E)
plot(as.factor(df4$us2),E)
plot(as.factor(df4$cs2),E)
plot(as.factor(df4$prey_presentation),E)
plot(as.factor(df4$way_of_learning),E)

library(nlme)
var<-varIdent()
lm4_all_var<-gls(attack_prob~species*us2*prey_presentation*cs2,data = df4,weights = )


# Check for collinearity
## Among ordinal predictors
par(mfrow=c(2,5))

plot(x = as.factor(df4$species),y = as.factor(df4$us2))
plot(x = as.factor(df4$species),y = as.factor(df4$prey_presentation))
plot(x = as.factor(df4$species),y = as.factor(df4$cs2),col=c("gray","navyblue","forestgreen","orange"))
plot(x = as.factor(df4$species),y = as.factor(df4$way_of_learning))

plot(x = as.factor(df4$us2),y = as.factor(df4$cs2),col=c("gray","navyblue","forestgreen","orange"))
plot(x = as.factor(df4$us2),y = as.factor(df4$prey_presentation))
plot(x = as.factor(df4$us2),y = as.factor(df4$way_of_learning))

plot(x = as.factor(df4$cs2),y = as.factor(df4$prey_presentation))
plot(x = as.factor(df4$cs2),y = as.factor(df4$way_of_learning))

plot(x = as.factor(df4$prey_presentation),y = as.factor(df4$way_of_learning))

# Stats
anova(lm4_all)

# Effect size
lm4_effectsize<-lm(attack_prob~species+cs2+us2+prey_presentation+way_of_learning,data = df4)

effectsize::standardize_parameters(lm4_effectsize,method = "refit",two_sd = TRUE)

params<-parameters::model_parameters(lm4_effectsize) 

effectsize::t_to_r(params$t[-1],df_error = params$df_error[-1])

par(mfrow=c(1,2))

# Plots for main effects
## Species
plot(as.factor(df4$species),df4$attack_prob)
means<-aggregate(attack_prob~species,data = df4,FUN = mean)
points(means[,2],col="red")

## US2
plot(as.factor(df4$us2),df4$attack_prob)
means<-aggregate(attack_prob~us2,data = df4,FUN = mean)
points(means[,2],col="red")

## CS2
plot(as.factor(df4$cs2),df4$attack_prob)
means<-aggregate(attack_prob~cs2,data = df4,FUN = mean)
points(means[,2],col="red")

## prey_presentation
plot(as.factor(df4$prey_presentation),df4$attack_prob)
means<-aggregate(attack_prob~prey_presentation,data = df4,FUN = mean)
points(means[,2],col="red")

## way_of_learning
plot(as.factor(df4$way_of_learning),df4$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = df4,FUN = mean)
points(means[,2],col="red")

## species*cs2
species<-levels(as.factor(df4$species))

Cc<-df4[df4$species==species[1],]
Ggd<-df4[df4$species==species[2],]
Pm<-df4[df4$species==species[3],]
Sv<-df4[df4$species==species[4],]

par(mfrow=c(2,2))

plot(x = as.factor(Cc$cs2),y = Cc$attack_prob)
means<-aggregate(attack_prob~cs2,data = Cc,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Ggd$cs2),y = Ggd$attack_prob)
means<-aggregate(attack_prob~cs2,data = Ggd,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Pm$cs2),y = Pm$attack_prob)
means<-aggregate(attack_prob~cs2,data = Pm,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Sv$cs2),y = Sv$attack_prob)
means<-aggregate(attack_prob~cs2,data = Sv,FUN = mean)
lines(means[,2],type = "b",col="red")

## species*way_of_learning
par(mfrow=c(2,2))

plot(x = as.factor(Cc$way_of_learning),y = Cc$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = Cc,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Ggd$way_of_learning),y = Ggd$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = Ggd,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Pm$way_of_learning),y = Pm$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = Pm,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(Sv$way_of_learning),y = Sv$attack_prob)
means<-aggregate(attack_prob~way_of_learning,data = Sv,FUN = mean)
lines(means[,2],type = "b",col="red")

## prey_presentation*cs2
p_sim<-subset(df4,subset = df4$prey_presentation=="simultaneous")
p_seq<-subset(df4,subset = df4$prey_presentation=="sequential")

par(mfrow=c(1,2))

plot(x = as.factor(p_seq$cs2),y = p_seq$attack_prob)
means<-aggregate(attack_prob~cs2,data = p_seq,FUN = mean)
lines(means[,2],type = "b",col="red")

plot(x = as.factor(p_sim$cs2),y = p_sim$attack_prob)
means<-aggregate(attack_prob~cs2,data = p_sim,FUN = mean)
lines(means[,2],type = "b",col="red")

# 5. Use direct learning data only, exclude data with evasiveness as US and single choice as prey_presentation
df5<-subset(df,subset = df$way_of_learning=="D"&df$us2!="evasiveness"&df$cs2 %in% c("pattern","color+pattern","color","color+other")&df$prey_presentation!="single choice")

# Check for collinearity among predictors
df5b<-df5

df5b$prey_presentation2<-ifelse(df5b$prey_presentation=="sequential",1,2)

df5b$species2<-99
df5b[df5b$species=="Cyanistes caeruleus",]$species2<-1
df5b[df5b$species=="Gallus gallus domesticus",]$species2<-2
df5b[df5b$species=="Parus major",]$species2<-3
df5b[df5b$species=="Sturnus vulgaris",]$species2<-4

df5b$cs3<-99
df5b[df5b$cs2=="color",]$cs3<-1
df5b[df5b$cs2=="color+other",]$cs3<-2
df5b[df5b$cs2=="color+pattern",]$cs3<-3
df5b[df5b$cs2=="pattern",]$cs3<-4

df5b$us3<-ifelse(df5b$us2=="bitterness",1,2)

# Collinearity among predictors
pairs(df5b[,c(16:19)],lower.panel = panel.smooth2,upper.panel = panel.cor,diag.panel = panel.hist)
corvif(df5b[,c(16:19)]) 

# Variance inflation factors
# 
# GVIF
# prey_presentation2 1.666829
# species2           1.450296
# cs3                1.300346
# us3                1.144392

# Use species, us2, cs2, and prey_presentation as predictors
library(glmmTMB)

lm5_zi<-glmmTMB::glmmTMB(asin(attack_prob)~species+us2+cs2+prey_presentation,data = df5,ziformula = ~us2,family = ziGamma(link = "log"),dispformula = ~species+us2+cs2)
# lm5_zi<-glmmTMB::glmmTMB(asin(attack_prob)~species+us2+cs2+prey_presentation+(1|study),data = df5,ziformula = ~us2,family = ziGamma(link = "log"),dispformula = ~species+us2+cs2)
# Including study as a random-effect factor caused  overfitting

# Model validation
library(DHARMa)

plot(E<-simulateResiduals(lm5_zi))

par(mfrow=c(2,2))
plotResiduals(E,as.factor(df5$species),rank = FALSE)
plotResiduals(E,as.factor(df5$us2),rank = FALSE)
plotResiduals(E,as.factor(df5$cs2),rank = FALSE)
plotResiduals(E,as.factor(df5$prey_presentation),rank = FALSE)
# All good

summary(lm5_zi)

# Plots
## Prey presentation
prey_pres<-ggplot(data = df5b[df5b$attack_prob!=0,],aes(x = prey_presentation,y = attack_prob))+
  geom_boxplot(aes(fill=prey_presentation,alpha = 0.5),show.legend = FALSE)+
  geom_point(aes(x = prey_presentation,y = attack_prob),position = position_jitter(width = 0.05))+
  scale_fill_manual(values = c("steelblue4","darkorange2"))+
  theme_bw()+
  xlab("Prey presentation")+ylab("Final attack probability")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15))+
  ggtitle("(C)")

## Reason for unprofitability
unpro<-ggplot(data = df5b,aes(x = us2,y = attack_prob))+
  geom_boxplot(aes(fill=us2,alpha = 0.5),show.legend = FALSE)+
  geom_point(aes(x = us2,y = attack_prob),position = position_jitter(width = 0.05))+
  geom_point(data = df5b[df5b$attack_prob==0,],aes(x = us2,y = attack_prob),color = "red",position = position_jitter(height = 0.005))+
  scale_fill_manual(values = c("steelblue4","darkorange2"))+
  theme_bw()+
  xlab("Reason for unprofitability")+ylab("Final attack probability")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15))+
  ggtitle("(A)")

## Level of bitterness (using a different dataset)
bitter<-ggplot(df3b,aes(x = us_concentration_cal,y = attack_prob,color = prey_presentation))+
  scale_color_manual(values = c("orange","cornflowerblue"))+
  theme_bw()+
  geom_point(aes(color = prey_presentation),size=3,alpha = 0.5)+
  geom_smooth(method = "lm",fill = NA)+
  theme(axis.title = element_text(size=15),
        axis.text = element_text(size=15),
        legend.position = "")+
  xlab("Calibrated bitterness (weight %)")+
  ylab("")+
  ggtitle("(B)")

## Type of signal/cue
signal_type<-ggplot(data = df5b[df5b$attack_prob!=0,],aes(x = cs2,y = attack_prob))+
  geom_boxplot(aes(fill=cs2,alpha = 0.5),show.legend = FALSE)+
  geom_point(aes(x = cs2,y = attack_prob),position = position_jitter(width = 0.05))+
  scale_fill_manual(values = c("steelblue4","darkorange2","seagreen","tomato1"))+
  theme_bw()+
  xlab("Type of signal/cue")+ylab("")+
  theme(axis.title = element_text(size=15),axis.text = element_text(size=15))+
  ggtitle("(D)")

tiff("figure3.tiff",width = 12,height = 12,units = "in",res = 300)
gridExtra::grid.arrange(unpro,bitter,prey_pres,signal_type,ncol=2)
dev.off()
