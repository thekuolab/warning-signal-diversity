library(tidyverse)
# library(gmailr)
require(parallel)

learning<-function(G,P,Abundance,lambda_L,alpha_F,K){
  
  # G is the magnitude of generalization, G{0,1}
  # P is a vector of current attack probabilities for OLD and NEW, P{0,1}
  # Abundance is a vector of starting abundances of each prey type
  # lambda_L is a vector of asymptotic attack rates for each prey type
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW
  # K is the carrying capacity of OLD and NEW. We assume that there's a palatable prey present with abundance = K 
  prey_type<-c("OLD","NEW","N")
  
  # relative abundance of OLD
  RA_OLD<-Abundance[1]/(sum(Abundance)+K[1])
  RA_NEW<-Abundance[2]/(sum(Abundance)+K[1])
  
  predation<-sample(x = c(prey_type),size = 1,prob = c(RA_OLD*P[1],RA_NEW*P[2],1-RA_OLD*P[1]-RA_NEW*P[2]))
  
  alpha_L<-0.5+abs(lambda_L-0.5)
  deltaPL<-alpha_L*(lambda_L-P)
  deltaPF<-alpha_F*(0.5-P)
  
  if (predation=="OLD"){
    
    outcome<-predation
    
    # Generalization effect from OLD to NEW
    deltaPL[2]<-G*deltaPL[1]
    
    # since OLD was attacked, forgetting on OLD does not happen
    deltaPF[1]<-0
    
    # if G > 0, encountering OLD will also negate forgetting for NEW
    if (G > 0){ deltaPF[2]<-0}
    
    # overall changes in attack probability  = learning + forgetting
    deltaP<-deltaPL+deltaPF
    
  } else if (predation=="NEW") {
    
    outcome<-predation
    
    # Generalization effect from NEW to OLD
    deltaPL[1]<-G*deltaPL[2]
    
    # since NEW was attacked, forgetting on NEW does not happen
    deltaPF[2]<-0
    
    # if G > 0, encountering NEW will also negate forgetting for OLD
    if (G > 0){deltaPF[1]<-0}
    
    # since b was attacked, forgetting on b does not happen
    deltaPF[2]<-0
    
    # overall changes in attack probability  = learning + forgetting
    deltaP<-deltaPL+deltaPF
    
  } else {
    
    outcome<-predation
    
    deltaPL[1:2]<-c(0,0)
    deltaPF[1:2]<-alpha_F*(0.5-P)
  }
  c(outcome,deltaPL,deltaPF)
  
  # Activate the next two lines to check if the function learning is working properly
  # predation_outcome<<-list(outcome,c(deltaP,deltaPL,deltaPF))
  # print(predation_outcome)
}


mimic<-function(G,P_OLD,Abundance,lambda_L,alpha_F,baby,K,GenTime,timestep){
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_F of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type, {0,400} 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is # of offspring each butterfly produces during reproduction phase, set to 10
  
  # K is the carrying capacity of OLD and NEW, set to 400
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  attack_data_mimic<<-data.frame(matrix(0,ncol=13,nrow=timestep))
  colnames(attack_data_mimic)<-c("outcome","P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There's no attack at time step 1.
  attack_data_mimic[1,]<-c("N",prey_P_ini,rep(0,4),Abundance,rep(0,4))
  
  for (i in 2:timestep){
    
    # No reproduction when time step is not a multiple of GenTime
    # Changes in species abundance is governed entirely by predation
    if (i%%GenTime!=0){
      
      # I. Surviving predation
      # 1. Filling in the columns of attack_data with simulation results for each time step
      attack_data_mimic[i,c(1,4:7)]<-learning(G = 0, P = as.numeric(attack_data_mimic[i-1,2:3]),Abundance = as.numeric(attack_data_mimic[i-1,8:9]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      
      # 2. If OLD gets attacked 
      if (attack_data_mimic[i,]$outcome=="OLD"){
        
        # 2.1 Record the attack
        attack_data_mimic[i,]$attack_OLD<--1
        
        # 2.2 If OLD gets attacked, it will also change the deltaPL_NEW by a factor of G    
        attack_data_mimic[i,]$deltaPL_NEW<-as.numeric(attack_data_mimic[i,]$deltaPL_OLD)*G
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
        
        # 3. Same process if NEW gets attacked
      } else if (attack_data_mimic[i,]$outcome=="NEW") {
        
        # 3.1 Record the attack
        attack_data_mimic[i,]$attack_NEW<--1
        
        # 3.2 Generalization effect from New to OLD
        attack_data_mimic[i,]$deltaPL_OLD<-as.numeric(attack_data_mimic[i,]$deltaPL_NEW)*G
        
        # 3.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        # 3.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential production of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 3.5 Update attack probabilities for the next time step 
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
        
        # 4. if neither OLD nor NEW gets attacked
      } else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[i,]$attack_OLD<-0
        attack_data_mimic[i,]$attack_NEW<-0
        
        # 4.2 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        #new_attack_prob<-ifelse(new_attack_prob<prey_lambdaL[2],prey_lambdaL[2],new_attack_prob)  
        
        # 4.3 Update attack probabilities for the next time step 
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
      } # this is where the if else statement for attack outcomes ends
      
      # Update abundance for OLD and NEW
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)+as.numeric(attack_data_mimic[i,]$attack_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)+as.numeric(attack_data_mimic[i,]$attack_NEW)
    }# This is where the if (i%%100 != 0) statement ends
    
    # II Reproduction
    # Reproduction happens when time step is a multiple of GenTime
    # Species reproduce, and then all individuals from the present generation die
    # The next generation is reconstituted from the offspring
    else if (i%%GenTime==0){
      
      # No predation happens 
      attack_data_mimic[i,]$outcome<-"N"
      
      # Attack probabilities remain the same from the previous generation
      attack_data_mimic[i,]$P_OLD<-attack_data_mimic[i-1,]$P_OLD
      attack_data_mimic[i,]$P_NEW<-attack_data_mimic[i-1,]$P_NEW
      
      # Population growth of the two species based on the logistic growth function
      attack_data_mimic[i,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/K[1]))
      attack_data_mimic[i,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)/K[2]))
      
      # Population size of the next generation = current population size + growth
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)
    }
  } # This is where the i loop ends  
  
  write_csv(attack_data_mimic,paste0("P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),".csv"))    
  # write_csv(attack_data_mimic,paste0("OLD",Abundance[1],"_NEW",Abundance[2],".csv"))
  # write_csv(attack_data_mimic,paste0("G",G,"_NEW",Abundance[2],".csv"))
} # This is where the function "mimic" ends

mimic_sensitivity<-function(G,P_OLD,Abundance,lambda_L,alpha_F,baby,K,GenTime,timestep){
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_F of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type, {0,400} 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is # of offspring each butterfly produces during reproduction phase, set to 10
  
  # K is the carrying capacity of OLD and NEW, set to 400
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  attack_data_mimic<<-data.frame(matrix(0,ncol=13,nrow=timestep))
  colnames(attack_data_mimic)<-c("outcome","P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There's no attack at time step 1.
  attack_data_mimic[1,]<-c("N",prey_P_ini,rep(0,4),Abundance,rep(0,4))
  
  for (i in 2:timestep){
    
    # No reproduction when time step is not a multiple of GenTime
    # Changes in species abundance is governed entirely by predation
    if (i%%GenTime!=0){
      
      # I. Surviving predation
      # 1. Filling in the columns of attack_data with simulation results for each time step
      attack_data_mimic[i,c(1,4:7)]<-learning(G = 0, P = as.numeric(attack_data_mimic[i-1,2:3]),Abundance = as.numeric(attack_data_mimic[i-1,8:9]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      
      # 2. If OLD gets attacked 
      if (attack_data_mimic[i,]$outcome=="OLD"){
        
        # 2.1 Record the attack
        attack_data_mimic[i,]$attack_OLD<--1
        
        # 2.2 If OLD gets attacked, it will also change the deltaPL_NEW by a factor of G    
        attack_data_mimic[i,]$deltaPL_NEW<-as.numeric(attack_data_mimic[i,]$deltaPL_OLD)*G
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
        
        # 3. Same process if NEW gets attacked
      } else if (attack_data_mimic[i,]$outcome=="NEW") {
        
        # 3.1 Record the attack
        attack_data_mimic[i,]$attack_NEW<--1
        
        # 3.2 Generalization effect from New to OLD
        attack_data_mimic[i,]$deltaPL_OLD<-as.numeric(attack_data_mimic[i,]$deltaPL_NEW)*G
        
        # 3.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        # 3.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential production of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 3.5 Update attack probabilities for the next time step 
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
        
        # 4. if neither OLD nor NEW gets attacked
      } else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[i,]$attack_OLD<-0
        attack_data_mimic[i,]$attack_NEW<-0
        
        # 4.2 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_data_mimic[i-1,]$P_OLD)+as.numeric(attack_data_mimic[i,]$deltaPL_OLD)+as.numeric(attack_data_mimic[i,]$deltaPF_OLD)
        new_attack_prob_NEW<-as.numeric(attack_data_mimic[i-1,]$P_NEW)+as.numeric(attack_data_mimic[i,]$deltaPL_NEW)+as.numeric(attack_data_mimic[i,]$deltaPF_NEW)
        
        #new_attack_prob<-ifelse(new_attack_prob<prey_lambdaL[2],prey_lambdaL[2],new_attack_prob)  
        
        # 4.3 Update attack probabilities for the next time step 
        attack_data_mimic[i,]$P_OLD<-new_attack_prob_OLD
        attack_data_mimic[i,]$P_NEW<-new_attack_prob_NEW
      } # this is where the if else statement for attack outcomes ends
      
      # Update abundance for OLD and NEW
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)+as.numeric(attack_data_mimic[i,]$attack_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)+as.numeric(attack_data_mimic[i,]$attack_NEW)
    }# This is where the if (i%%100 != 0) statement ends
    
    # II Reproduction
    # Reproduction happens when time step is a multiple of GenTime
    # Species reproduce, and then all individuals from the present generation die
    # The next generation is reconstituted from the offspring
    else if (i%%GenTime==0){
      
      # No predation happens 
      attack_data_mimic[i,]$outcome<-"N"
      
      # Attack probabilities remain the same from the previous generation
      attack_data_mimic[i,]$P_OLD<-attack_data_mimic[i-1,]$P_OLD
      attack_data_mimic[i,]$P_NEW<-attack_data_mimic[i-1,]$P_NEW
      
      # Population growth of the two species based on the logistic growth function
      attack_data_mimic[i,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/K[1]))
      attack_data_mimic[i,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)/K[2]))
      
      # Population size of the next generation = current population size + growth
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)
    }
  } # This is where the i loop ends  
  
  write_csv(attack_data_mimic,paste0("P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),".csv"))    
  # write_csv(attack_data_mimic,paste0("OLD",Abundance[1],"_NEW",Abundance[2],".csv"))
  # write_csv(attack_data_mimic,paste0("G",G,"_NEW",Abundance[2],".csv"))
} # This is where the function "mimic_sensitivity" ends

mimic_social<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 100 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is the carrying capacity of OLD and NEW (should be set as equal to the initial abundance)
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=timestep))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=timestep))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  for (i in 2:timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[i,]$attacker<-attacker
      attack_data_mimic[i,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][i-1,1:2]),Abundance = as.numeric(attack_data_mimic[i-1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[i,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[i,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][i,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][i,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_OLD)
          attack_prob_predator[j][[1]][i,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][i,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_NEW)
          
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)

          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
          
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[i,]$attack_NEW<--1
        
        # 3.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 3.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 3.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 3.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 3.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][i,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][i,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_OLD)
          attack_prob_predator[j][[1]][i,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][i,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_NEW)
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)          
          
          attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[i,]$attack_OLD<-0
        attack_data_mimic[i,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 4.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 4.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 4.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 4.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)          
          
          attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)+as.numeric(attack_data_mimic[i,]$attack_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)+as.numeric(attack_data_mimic[i,]$attack_NEW)
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[i,]$outcome<-"N"
      attack_data_mimic[i,]$attacker<-"NA"
      attack_data_mimic[i,]$learner<-"NA"
      attack_data_mimic[i,]$attack_OLD<-0
      attack_data_mimic[i,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][i,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][i,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][i,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][i,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][i-1,]$P_OLD)+attack_prob_predator[n][[1]][i,]$deltaPL_OLD+attack_prob_predator[n][[1]][i,]$deltaPF_OLD
        attack_prob_predator[n][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][i-1,]$P_NEW)+attack_prob_predator[n][[1]][i,]$deltaPL_NEW+attack_prob_predator[n][[1]][i,]$deltaPF_NEW
      }
      
      
      # Reproduction
      # attack_data_mimic[i,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/K[1]))
      # attack_data_mimic[i,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)/K[2]))
      popsize<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i-1,]$abnc_NEW)
      
      birth_total<-round(baby*popsize*(1-popsize/K))
      
      if (popsize==0){
        
        attack_data_mimic[i,]$birth_OLD<-0
        attack_data_mimic[i,]$birth_NEW<-0
        
      } else if (popsize!=0){
      
      attack_data_mimic[i,]$birth_OLD<-round(birth_total*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/popsize)
      attack_data_mimic[i,]$birth_NEW<-birth_total-as.numeric(attack_data_mimic[i,]$birth_OLD)
      }
      
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)
      
    }
    
  write_csv(attack_data_mimic,paste0("P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))  
  }# end of timestep for loop
  
  # write_csv(attack_data_mimic,paste0("P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
  
} # end of the mimic_social function  


## Faster version of the simulation function
## Only retains data from the current and previous time steps so runs much faster
 

mimic_social_fast<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 10 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is a vector containing the carrying capacity of OLD and NEW (should be set as equal to the initial abundance)
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  # 3 rows: 1st row for initial condition
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=2))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=2))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  i<<-2
  
  while (i < timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[2,]$attacker<-attacker
      attack_data_mimic[2,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][1,1:2]),Abundance = as.numeric(attack_data_mimic[1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[2,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_NEW<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[2,]$attack_OLD<-0
        attack_data_mimic[2,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)+as.numeric(attack_data_mimic[2,]$attack_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)+as.numeric(attack_data_mimic[2,]$attack_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      # clean up the second row for the next step
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
      
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[2,]$outcome<-"N"
      attack_data_mimic[2,]$attacker<-"NA"
      attack_data_mimic[2,]$learner<-"NA"
      attack_data_mimic[2,]$attack_OLD<-0
      attack_data_mimic[2,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][2,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][1,]$P_OLD)+attack_prob_predator[n][[1]][2,]$deltaPL_OLD+attack_prob_predator[n][[1]][2,]$deltaPF_OLD
        attack_prob_predator[n][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][1,]$P_NEW)+attack_prob_predator[n][[1]][2,]$deltaPL_NEW+attack_prob_predator[n][[1]][2,]$deltaPF_NEW
        
        attack_prob_predator[n][[1]][1,]<-attack_prob_predator[n][[1]][2,]
      }
      
      
      # Reproduction
      # First check if both prey have become extinct. If yes, simulation terminates immediately
      popsize<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[1,]$abnc_NEW)
      if (popsize==0){i<-timestep} else {
        
        attack_data_mimic[2,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[1,]$abnc_OLD)/K[1]))
        attack_data_mimic[2,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[1,]$abnc_NEW)/K[2]))
        
        attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)
        attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)
        
        attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
        attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
        
        attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
      }
    }  
    
    # write_csv(attack_data_mimic,paste0("~/Desktop/P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))  
    i<<-i+1  
  }# end of timestep while loop
  write_csv(attack_data_mimic,paste0("rare",Abundance[2],"_K2",K[2],"_P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
} # end of the mimic_social_fast function  

# mimic_social_fast_pto
## 1. "pto" stands for predator turnover: 10% of the predators die each generation and are replaced by naive predators
## 2. Predator generation time set to either the same as prey generation time or 10X longer
## 3. Other setting same as mimic_social_fast

# mimic_social_fast2:
## 1. Prey are annual species, so a generation is reconstituted by new individual born after the previous generation
## based on Aubier and Sherratt 2015. Evolution

mimic_social_fast2<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 20 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is the carrying capacity of OLD and NEW (should be set as equal to the initial abundance)
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  # 3 rows: 1st row for initial condition
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=2))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=2))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  i<<-2
  
  while (i < timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[2,]$attacker<-attacker
      attack_data_mimic[2,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][1,1:2]),Abundance = as.numeric(attack_data_mimic[1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[2,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_NEW<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[2,]$attack_OLD<-0
        attack_data_mimic[2,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)+as.numeric(attack_data_mimic[2,]$attack_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)+as.numeric(attack_data_mimic[2,]$attack_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      # clean up the second row for the next step
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
      
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[2,]$outcome<-"N"
      attack_data_mimic[2,]$attacker<-"NA"
      attack_data_mimic[2,]$learner<-"NA"
      attack_data_mimic[2,]$attack_OLD<-0
      attack_data_mimic[2,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][2,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][1,]$P_OLD)+attack_prob_predator[n][[1]][2,]$deltaPL_OLD+attack_prob_predator[n][[1]][2,]$deltaPF_OLD
        attack_prob_predator[n][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][1,]$P_NEW)+attack_prob_predator[n][[1]][2,]$deltaPL_NEW+attack_prob_predator[n][[1]][2,]$deltaPF_NEW
        
        attack_prob_predator[n][[1]][1,]<-attack_prob_predator[n][[1]][2,]
      }
      
      
      # Reproduction
      attack_data_mimic[2,]$birth_OLD<-round(as.numeric(attack_data_mimic[1,]$abnc_OLD)*baby*(1/(1+(baby-1)*(as.numeric(attack_data_mimic[1,]$abnc_OLD)/K[1]))))
      attack_data_mimic[2,]$birth_NEW<-round(as.numeric(attack_data_mimic[1,]$abnc_NEW)*baby*(1/(1+(baby-1)*(as.numeric(attack_data_mimic[1,]$abnc_NEW)/K[2]))))
      
      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[2,]$birth_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[2,]$birth_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
    }
    
    
    
    i<<-i+1  
  }# end of timestep while loop
  
  write_csv(attack_data_mimic,paste0(Abundance[2],"_K2",K[2],"_P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
  
} # end of the mimic_social_fast2 function 

# mimic_social_fast3:
## 1. The two prey share the same carrying capacity that equals to N1+N2
## This scenario allows meaningful calculation of selection coefficient against the rare prey
mimic_social_fast3<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 10 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is the carrying capacity of OLD and NEW, which equals the sum of initial abundance of the two prey
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  # 3 rows: 1st row for initial condition
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=2))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=2))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  i<<-2
  
  while (i < timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[2,]$attacker<-attacker
      attack_data_mimic[2,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][1,1:2]),Abundance = as.numeric(attack_data_mimic[1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[2,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_NEW<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[2,]$attack_OLD<-0
        attack_data_mimic[2,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)+as.numeric(attack_data_mimic[2,]$attack_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)+as.numeric(attack_data_mimic[2,]$attack_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      # clean up the second row for the next step
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
      
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[2,]$outcome<-"N"
      attack_data_mimic[2,]$attacker<-"NA"
      attack_data_mimic[2,]$learner<-"NA"
      attack_data_mimic[2,]$attack_OLD<-0
      attack_data_mimic[2,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][2,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][1,]$P_OLD)+attack_prob_predator[n][[1]][2,]$deltaPL_OLD+attack_prob_predator[n][[1]][2,]$deltaPF_OLD
        attack_prob_predator[n][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][1,]$P_NEW)+attack_prob_predator[n][[1]][2,]$deltaPL_NEW+attack_prob_predator[n][[1]][2,]$deltaPF_NEW
        
        attack_prob_predator[n][[1]][1,]<-attack_prob_predator[n][[1]][2,]
      }
      
      
      # Reproduction
      # First check if both prey have become extinct. If yes, simulation terminates immediately
      popsize<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[1,]$abnc_NEW)
      birth_total<-round(baby*popsize*(1-popsize/K))
      
      if (popsize==0){
        
      attack_data_mimic[2,]$birth_OLD<-0
      attack_data_mimic[2,]$birth_NEW<-0
      
      i<-timestep
      
      } else {
      
      attack_data_mimic[2,]$birth_OLD<-round(birth_total*as.numeric(attack_data_mimic[1,]$abnc_OLD)/popsize)
      attack_data_mimic[2,]$birth_NEW<-birth_total-as.numeric(attack_data_mimic[2,]$birth_OLD)
      
      }

      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
    }
    
    # write_csv(attack_data_mimic,paste0("~/Desktop/P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
    i<<-i+1
  }# end of timestep while loop
  
  write_csv(attack_data_mimic,paste0("rare",Abundance[2],"_K",K,"_P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
  
} # end of the mimic_social_fast3 function 

mimic_social_pto<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,GenTime_P,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 100 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is the carrying capacity of OLD and NEW (should be set as equal to the initial abundance)
  
  # GenTime is the generation time for OLD and NEW
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=timestep))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=timestep))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  for (i in 2:timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[i,]$attacker<-attacker
      attack_data_mimic[i,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][i-1,1:2]),Abundance = as.numeric(attack_data_mimic[i-1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[i,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[i,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][i,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][i,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_OLD)
          attack_prob_predator[j][[1]][i,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][i,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_NEW)
          
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)
          
          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
          
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[i,]$attack_NEW<--1
        
        # 3.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 3.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 3.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 3.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 3.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][i,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][i,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_OLD)
          attack_prob_predator[j][[1]][i,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][i,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][i,]$P_NEW)
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)          
          
          attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[i,]$attack_OLD<-0
        attack_data_mimic[i,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][i,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][i,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][i,]$deltaPF_NEW<-deltaPF_NEW
        
        # 4.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][i-1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 4.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 4.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][i,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][i,]$P_NEW<-new_attack_prob_NEW
        
        # 4.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][i,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][i-1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_OLD)+attack_prob_predator[j][[1]][i,]$deltaPL_OLD+attack_prob_predator[j][[1]][i,]$deltaPF_OLD
          attack_prob_predator[j][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][i-1,]$P_NEW)+attack_prob_predator[j][[1]][i,]$deltaPL_NEW+attack_prob_predator[j][[1]][i,]$deltaPF_NEW
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          # attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_OLD)
          # attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][i-1,]$P_NEW)          
          
          attack_prob_predator[k][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[k][[1]][i,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_OLD)+attack_prob_predator[k][[1]][i,]$deltaPL_OLD+attack_prob_predator[k][[1]][i,]$deltaPF_OLD
          attack_prob_predator[k][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][i-1,]$P_NEW)+attack_prob_predator[k][[1]][i,]$deltaPL_NEW+attack_prob_predator[k][[1]][i,]$deltaPF_NEW
          
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)+as.numeric(attack_data_mimic[i,]$attack_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)+as.numeric(attack_data_mimic[i,]$attack_NEW)
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[i,]$outcome<-"N"
      attack_data_mimic[i,]$attacker<-"NA"
      attack_data_mimic[i,]$learner<-"NA"
      attack_data_mimic[i,]$attack_OLD<-0
      attack_data_mimic[i,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][i,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][i,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][i,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][i,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][i,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][i-1,]$P_OLD)+attack_prob_predator[n][[1]][i,]$deltaPL_OLD+attack_prob_predator[n][[1]][i,]$deltaPF_OLD
        attack_prob_predator[n][[1]][i,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][i-1,]$P_NEW)+attack_prob_predator[n][[1]][i,]$deltaPL_NEW+attack_prob_predator[n][[1]][i,]$deltaPF_NEW
      }
      
      
      # Reproduction
      # attack_data_mimic[i,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/K[1]))
      # attack_data_mimic[i,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)/K[2]))
      
      attack_data_mimic[i,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)/K[1]))
      attack_data_mimic[i,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[i-1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)/K[2]))
      
      attack_data_mimic[i,]$abnc_OLD<-as.numeric(attack_data_mimic[i-1,]$abnc_OLD)+as.numeric(attack_data_mimic[i,]$birth_OLD)
      attack_data_mimic[i,]$abnc_NEW<-as.numeric(attack_data_mimic[i-1,]$abnc_NEW)+as.numeric(attack_data_mimic[i,]$birth_NEW)
      
      if (i%%GenTime_P==0){
        
        predator_died<-sample(c(1:20),size = 2,replace = FALSE)
        predator_survived<-c(1:20)[-predator_died]
        
        attack_data_mimic[i,]$attacker<-paste(predator_died,collapse = ".")
        
        for (x in predator_died){
          
          attack_prob_predator[x][[1]][i,]$P_OLD<-0.5
          attack_prob_predator[x][[1]][i,]$P_NEW<-0.5
          attack_prob_predator[x][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[x][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[x][[1]][i,]$deltaPL_NEW<-0
          attack_prob_predator[x][[1]][i,]$deltaPF_NEW<-0
          
        }
        
        for (y in predator_survived){
          
          attack_prob_predator[y][[1]][i,]$P_OLD<-attack_prob_predator[y][[1]][i-1,]$P_OLD
          attack_prob_predator[y][[1]][i,]$P_NEW<-attack_prob_predator[y][[1]][i-1,]$P_NEW
          attack_prob_predator[y][[1]][i,]$deltaPL_OLD<-0
          attack_prob_predator[y][[1]][i,]$deltaPF_OLD<-0
          attack_prob_predator[y][[1]][i,]$deltaPL_NEW<-0
          attack_prob_predator[y][[1]][i,]$deltaPF_NEW<-0
          
        }
      }  # end of if (i%%GenTime_P==0) statement
      
      
    } # end of if (i%%GenTime==0) statement
    
    
  }# end of timestep for loop
  for (z in 1:20){
    
    write.csv(attack_prob_predator[z][[1]],paste0(z,".csv"))
    # Outputting changes in attack prob for each predator for validation
  }
  write_csv(attack_data_mimic,paste0("P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))  
} # end of the mimic_social_pto_test function  

mimic_social_fast_pto<-function(G,P_OLD,Abundance,lambda_L,social_spread,social_efficacy,inter_learning,alpha_F,baby,K,GenTime,GenTime_P,turnover,timestep){
  
  # This mimic function incorporates social learning. There are two predator species, each of which has 10 individuals
  
  # G is the magnitude of generalization, G{0,1}
  
  # P_OLD is the attack probability for the local prey at t0, set equal to lambda_L of the OLD prey 
  # if predators have already learned to avoid them at the beginning of the simulation
  
  # Abundance is a vector of abundances of each prey type 
  
  # lambda_L is a vector of asymptotic attack rates for each prey type
  
  # social_spread is input as a percentage {0,1}, indicating the probability of other predators observing the foraging outcome
  
  # social_efficacy is input as a percentage {0,1}, indicating the "discount" in attack probability relative to 1st-hand experience
  
  # inter_learning: if TRUE, social learning occurs across species, if FALSE, social learning occurs only within species
  
  # alpha_F is the forgetting rate, set to be the same for OLD and NEW prey
  
  # baby is the intrinsic rate of increase during reproduction
  
  # K is the carrying capacity of OLD and NEW (should be set as equal to the initial abundance)
  
  # GenTime is the generation time for OLD and NEW
  
  # GenTime_P is the generation time for predators, necessary for modelling predator turnover
  
  # timestep is the # of total time steps for the simulation
  
  # Initial attack probability toward OLD and NEW
  prey_P_ini<-c(P_OLD,0.5-(0.5-P_OLD)*G)
  
  # 3 rows: 1st row for initial condition
  attack_data_mimic<<-data.frame(matrix(0,ncol=9,nrow=2))
  colnames(attack_data_mimic)<-c("outcome","attacker","learner","abnc_OLD","abnc_NEW","birth_OLD","birth_NEW","attack_OLD","attack_NEW")
  
  # Initial condition. There is no attack at time step 1
  attack_data_mimic[1,]<-c("N","NA","NA",Abundance[1],Abundance[2],rep(0,4))
  
  attack_prob_predator<<-as.list(NULL)
  
  for (b in 1:20){
    
    attack_prob_predator[[b]]<-data.frame(matrix(0,ncol=6,nrow=2))
    attack_prob_predator[[b]][1,]<-c(P_OLD,0.5,rep(0,4))
    colnames(attack_prob_predator[[b]])<-c("P_OLD","P_NEW","deltaPL_OLD","deltaPL_NEW","deltaPF_OLD","deltaPF_NEW")
    
  }
  
  i<<-2
  
  while (i < timestep){
    
    if (i%%GenTime!=0){
      
      # No reproduction when time step is not a multiple of GenTime
      # Changes in species abundance is governed entirely by predation
      
      
      # I. Surviving predation
      # Randomly choose the attacker
      attacker<-sample(x = c(1:20),size = 1)
      
      # Randomly choose the learner based on social_spread
      ## If social learning can occur across species
      if (inter_learning==TRUE){
        
        # Create a potential learner pool
        learner_pool<-c(1:20)[-attacker]
        # Choose the birds that will observe the foraging outcome
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        
        # If social learning can only occur within species  
        # If the attacker belongs to species 1 (birds 1-10)
      } else if (inter_learning==FALSE & attacker<=10){  
        learner_pool<-c(1:20)[-attacker][1:9]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
        # If the attacker belongs to species 2 (birds 11-20)   
      } else if (inter_learning==FALSE & attacker>10) {
        learner_pool<-c(1:20)[-attacker][11:19]
        learners<-sample(x = learner_pool,size = round(social_spread*length(learner_pool)),replace = FALSE)
      }
      
      # Record attacker and learner in attack_data_mimic
      attack_data_mimic[2,]$attacker<-attacker
      attack_data_mimic[2,]$learner<-paste(learners,collapse = ".")
      
      # 1. Predation outcome and the associated change in attack prob
      result<-learning(G = G, P = as.numeric(attack_prob_predator[attacker][[1]][1,1:2]),Abundance = as.numeric(attack_data_mimic[1,4:5]),lambda_L = lambda_L,alpha_F = alpha_F,K = K)
      outcome<-result[1]
      attack_data_mimic[2,]$outcome<-result[1]
      deltaPL_OLD<-as.numeric(result[2])
      deltaPL_NEW<-as.numeric(result[3])
      deltaPF_OLD<-as.numeric(result[4])
      deltaPF_NEW<-as.numeric(result[5])
      
      # 2. If OLD gets attacked 
      if (outcome=="OLD"){
        
        # 2.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_OLD<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_OLD*G
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-deltaPL_OLD*social_efficacy
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-0
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)          
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][1,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # Alternative: others do not forget when it is not their turns to feed
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        # 3. Same process if "NEW" is attacked  
      } else if (outcome=="NEW"){ 
        
        # 3.1 Record the attack in attack_data_mimic
        attack_data_mimic[2,]$attack_NEW<--1
        
        # 2.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW*G
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-deltaPL_NEW*social_efficacy
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-0         
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          # This step is to safeguard the potential result of negative attack prob due to social learning
          attack_prob_predator[j][[1]][2,]$P_OLD<-ifelse(attack_prob_predator[j][[1]][2,]$P_OLD<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_OLD)
          attack_prob_predator[j][[1]][2,]$P_NEW<-ifelse(attack_prob_predator[j][[1]][2,]$P_NEW<lambda_L[1],lambda_L[1],attack_prob_predator[j][[1]][2,]$P_NEW)
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
          
        }# end of the others loop
        
        
      }  else {
        
        # 4.1 No attack is recorded
        attack_data_mimic[2,]$attack_OLD<-0
        attack_data_mimic[2,]$attack_NEW<-0
        
        # 4.2 Updating numbers in the chosen predator's data frame
        attack_prob_predator[attacker][[1]][2,]$deltaPL_OLD<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPL_NEW<-deltaPL_NEW
        attack_prob_predator[attacker][[1]][2,]$deltaPF_OLD<-deltaPF_OLD
        attack_prob_predator[attacker][[1]][2,]$deltaPF_NEW<-deltaPF_NEW
        
        # 2.3 Calculate the attack probabilities = attack probability from the previous time step + learning + forgetting
        new_attack_prob_OLD<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_OLD)+deltaPL_OLD+deltaPF_OLD
        new_attack_prob_NEW<-as.numeric(attack_prob_predator[attacker][[1]][1,]$P_NEW)+deltaPL_NEW+deltaPF_NEW
        
        # 2.4 Accept the new attack probability if it is >= lambda_L, otherwise the new attack probability is set to the asymptotic value
        # This step is to safeguard the potential result of negative attack probabilities due to generalization effects
        new_attack_prob_OLD<-ifelse(new_attack_prob_OLD<lambda_L[1],lambda_L[1],new_attack_prob_OLD)
        new_attack_prob_NEW<-ifelse(new_attack_prob_NEW<lambda_L[2],lambda_L[2],new_attack_prob_NEW)
        
        # 2.5 Update attack probabilities for the next time step for the attacker
        attack_prob_predator[attacker][[1]][2,]$P_OLD<-new_attack_prob_OLD
        attack_prob_predator[attacker][[1]][2,]$P_NEW<-new_attack_prob_NEW
        
        attack_prob_predator[attacker][[1]][1,1:2]<-attack_prob_predator[attacker][[1]][2,1:2]
        attack_prob_predator[attacker][[1]][2,]<-0
        
        # 2.6 Social learning:
        # Modify learners' attack prob based on social learning efficacy
        ## Updating numbers for learners
        for (j in learners){
          
          # deltaPL_OLD changes based on social learning efficacy
          attack_prob_predator[j][[1]][2,]$deltaPL_OLD<-0
          # deltaPL_NEW is 0, as generalization is assumed to not occur via social learning
          attack_prob_predator[j][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD is 0 b/c social learning
          attack_prob_predator[j][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_OLD) 
          # deltaPF_NEW is calculated based on the equation
          attack_prob_predator[j][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[j][[1]][1,]$P_NEW)        
          
          # Updating P_OLD and P_NEW
          attack_prob_predator[j][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[j][[1]][1,]$P_OLD)+attack_prob_predator[j][[1]][2,]$deltaPL_OLD+attack_prob_predator[j][[1]][2,]$deltaPF_OLD
          attack_prob_predator[j][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[j][[1]][1,]$P_NEW)+attack_prob_predator[j][[1]][2,]$deltaPL_NEW+attack_prob_predator[j][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[j][[1]][1,1:2]<-attack_prob_predator[j][[1]][2,1:2]
          attack_prob_predator[j][[1]][2,]<-0
          
        }# end of the learners loop
        
        others<-c(1:20)[-c(attacker,learners)]
        ## P_OLD and P_NEW do not change for all other predators
        for (k in others){
          
          # deltaPL_OLD and NEW are 0
          attack_prob_predator[k][[1]][2,]$deltaPL_OLD<-0
          attack_prob_predator[k][[1]][2,]$deltaPL_NEW<-0
          # deltaPF_OLD and NEW are calculated based on the equation
          attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_OLD)
          attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-alpha_F*(0.5-attack_prob_predator[k][[1]][1,]$P_NEW)
          
          # attack_prob_predator[k][[1]][2,]$deltaPF_OLD<-0
          # attack_prob_predator[k][[1]][2,]$deltaPF_NEW<-0   
          
          # # Updating P_OLD and P_NEW
          attack_prob_predator[k][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[k][[1]][1,]$P_OLD)+attack_prob_predator[k][[1]][2,]$deltaPL_OLD+attack_prob_predator[k][[1]][2,]$deltaPF_OLD
          attack_prob_predator[k][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[k][[1]][1,]$P_NEW)+attack_prob_predator[k][[1]][2,]$deltaPL_NEW+attack_prob_predator[k][[1]][2,]$deltaPF_NEW
          
          attack_prob_predator[k][[1]][1,1:2]<-attack_prob_predator[k][[1]][2,1:2]
          attack_prob_predator[k][[1]][2,]<-0
        }# end of the others loop
        
      } # end of the if else statement for attack outcomes
      
      # Update abundance for OLD and NEW
      attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)+as.numeric(attack_data_mimic[2,]$attack_OLD)
      attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)+as.numeric(attack_data_mimic[2,]$attack_NEW)
      
      attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
      attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
      
      # clean up the second row for the next step
      attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
      
    } # end of the i%%Gentime!=0 if statement
    else if (i%%GenTime==0){
      
      # No predation, no attacker and learner chosen
      attack_data_mimic[2,]$outcome<-"N"
      attack_data_mimic[2,]$attacker<-"NA"
      attack_data_mimic[2,]$learner<-"NA"
      attack_data_mimic[2,]$attack_OLD<-0
      attack_data_mimic[2,]$attack_NEW<-0
      
      # Updating numbers for all predators
      for (n in 1:20){
        
        attack_prob_predator[n][[1]][2,]$deltaPL_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPL_NEW<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_OLD<-0
        attack_prob_predator[n][[1]][2,]$deltaPF_NEW<-0
        
        attack_prob_predator[n][[1]][2,]$P_OLD<-as.numeric(attack_prob_predator[n][[1]][1,]$P_OLD)+attack_prob_predator[n][[1]][2,]$deltaPL_OLD+attack_prob_predator[n][[1]][2,]$deltaPF_OLD
        attack_prob_predator[n][[1]][2,]$P_NEW<-as.numeric(attack_prob_predator[n][[1]][1,]$P_NEW)+attack_prob_predator[n][[1]][2,]$deltaPL_NEW+attack_prob_predator[n][[1]][2,]$deltaPF_NEW
        
        attack_prob_predator[n][[1]][1,]<-attack_prob_predator[n][[1]][2,]
      }
      
      
      # Reproduction
      # First check if both prey have become extinct. If yes, simulation terminates immediately
      popsize<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[1,]$abnc_NEW)
      if (popsize==0){i<-timestep} else {
        
        attack_data_mimic[2,]$birth_OLD<-round(baby*as.numeric(attack_data_mimic[1,]$abnc_OLD)*(1-as.numeric(attack_data_mimic[1,]$abnc_OLD)/K[1]))
        attack_data_mimic[2,]$birth_NEW<-round(baby*as.numeric(attack_data_mimic[1,]$abnc_NEW)*(1-as.numeric(attack_data_mimic[1,]$abnc_NEW)/K[2]))
        
        attack_data_mimic[2,]$abnc_OLD<-as.numeric(attack_data_mimic[1,]$abnc_OLD)+as.numeric(attack_data_mimic[2,]$birth_OLD)
        attack_data_mimic[2,]$abnc_NEW<-as.numeric(attack_data_mimic[1,]$abnc_NEW)+as.numeric(attack_data_mimic[2,]$birth_NEW)
        
        attack_data_mimic[1,]$abnc_OLD<-attack_data_mimic[2,]$abnc_OLD
        attack_data_mimic[1,]$abnc_NEW<-attack_data_mimic[2,]$abnc_NEW
        
        attack_data_mimic[2,]<-c("N",rep("NA",2),rep(0,6))
        
        if (i%%GenTime_P==0){
          
          predator_died<-sample(c(1:20),size = turnover,replace = FALSE)
          predator_survived<-c(1:20)[-predator_died]
          
          # Resetting the attack prob of "dead" predators to 0.5 for both prey
          for (x in predator_died){
            
            attack_prob_predator[x][[1]][1,]$P_OLD<-0.5
            attack_prob_predator[x][[1]][1,]$P_NEW<-0.5
            attack_prob_predator[x][[1]][2,]$deltaPL_OLD<-0
            attack_prob_predator[x][[1]][2,]$deltaPF_OLD<-0
            attack_prob_predator[x][[1]][2,]$deltaPL_NEW<-0
            attack_prob_predator[x][[1]][2,]$deltaPF_NEW<-0
            
          }
          
          # No learning occurs for surviving predators during this step
          for (y in predator_survived){
            
            attack_prob_predator[y][[1]][2,]$P_OLD<-attack_prob_predator[y][[1]][1,]$P_OLD
            attack_prob_predator[y][[1]][2,]$P_NEW<-attack_prob_predator[y][[1]][1,]$P_NEW
            attack_prob_predator[y][[1]][2,]$deltaPL_OLD<-0
            attack_prob_predator[y][[1]][2,]$deltaPF_OLD<-0
            attack_prob_predator[y][[1]][2,]$deltaPL_NEW<-0
            attack_prob_predator[y][[1]][2,]$deltaPF_NEW<-0
            
          }
          
          
        }  # end of if (i%%GenTime_P==0) statement
      } # end of if (popsize==0){i<-timestep} else statement
    } # end of  else if (i%%GenTime==0) statement
    # write_csv(attack_data_mimic,paste0("~/Desktop/P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))  
    i<<-i+1  
  }# end of timestep while loop
  write_csv(attack_data_mimic,paste0("rare",Abundance[2],"_K2",K[2],"_P",round(lambda_L[1],3),"_F",round(alpha_F[1],3),"_baby_",round(baby,3),"_sspr_",social_spread,"_seff_",social_efficacy,".csv"))
} # end of the mimic_social_fast function  