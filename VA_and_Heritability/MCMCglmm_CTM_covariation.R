

library(MCMCglmm)

##see DS_glmmMCMC_survival.R code for how to make family pedigree file

##make pedigree file
setwd("~/UCDavis/FCCL/Spawning_2020/CTMs")
setwd("~/UCDavis/Smelt_heritability")
ped <- read.delim2("DS_ped_family_offspring.txt", header=T)
ped$animal<-as.factor(ped$animal)
ped$FATHER<-as.factor(ped$FATHER)
ped$MOTHER<-as.factor(ped$MOTHER)
colnames(ped) <- c("Family", "MOTHER", "FATHER")
ped2 <- ped

ped$Rep <- rep("1",nrow(ped))
ped$animal <- paste(ped$Family, ped$Rep, sep = "_")

ped2$Rep <- rep("2",nrow(ped2))
ped2$animal <- paste(ped2$Family, ped2$Rep, sep = "_")

ped_offspring <- rbind(ped, ped2)
ped_offspring <- ped_offspring[c("animal", "MOTHER", "FATHER")]
#write.table(ped_offspring, file = "DS_ped_family_offspring_covariation.txt", sep = "\t", row.names=FALSE, quote = FALSE)
##manually add parents back to file from family ped file

Ped <- read.delim2("DS_ped_family_parents_offspring_covariation.txt", header = T)

#read in phenotype data
ctm<- read.delim2("ctm_DI_parents_r2.txt", header = T)
ctm$CTM <- as.numeric(ctm$CTM)
ctm_1_15 <- subset(ctm, System=="1-15")
ctm_2_15 <- subset(ctm, System=="2-15")
ctm_1_18 <- subset(ctm, System=="1-18")
ctm_2_18 <- subset(ctm, System=="2-18")

mean_ctm_1_15 <- aggregate(ctm_1_15$CTM, list(ctm_1_15$Family), mean)
colnames(mean_ctm_1_15) <- c("Family", "CTM_15")
mean_ctm_1_15$System <- rep("1-15",nrow(mean_ctm_1_15))
mean_ctm_1_15$Rear_temp <- rep("15",nrow(mean_ctm_1_15))
mean_ctm_1_15$Rep <- rep("1",nrow(mean_ctm_1_15))
mean_ctm_1_15$animal <- paste(mean_ctm_1_15$Family, mean_ctm_1_15$Rep, sep = "_")

mean_ctm_1_18 <- aggregate(ctm_1_18$CTM, list(ctm_1_18$Family), mean)
colnames(mean_ctm_1_18) <- c("Family", "CTM_18")
mean_ctm_1_18$System <- rep("1-18",nrow(mean_ctm_1_18))
mean_ctm_1_18$Rear_temp <- rep("18",nrow(mean_ctm_1_18))
mean_ctm_1_18$Rep <- rep("1",nrow(mean_ctm_1_18))
mean_ctm_1_18$animal <- paste(mean_ctm_1_18$Family, mean_ctm_1_18$Rep, sep = "_")

mean_ctm_rep1 <- merge(mean_ctm_1_15, mean_ctm_1_18, by="animal")
mean_ctm_rep1 <- mean_ctm_rep1[c("animal", "Family.x", "CTM_15", "CTM_18", "Rep.x")]

mean_ctm_2_15 <- aggregate(ctm_2_15$CTM, list(ctm_2_15$Family), mean)
colnames(mean_ctm_2_15) <- c("Family", "CTM_15")
mean_ctm_2_15$System <- rep("2-15",nrow(mean_ctm_2_15))
mean_ctm_2_15$Rear_temp <- rep("15",nrow(mean_ctm_2_15))
mean_ctm_2_15$Rep <- rep("2",nrow(mean_ctm_2_15))
mean_ctm_2_15$animal <- paste(mean_ctm_2_15$Family, mean_ctm_2_15$Rep, sep = "_")

mean_ctm_2_18 <- aggregate(ctm_2_18$CTM, list(ctm_2_18$Family), mean)
colnames(mean_ctm_2_18) <- c("Family", "CTM_18")
mean_ctm_2_18$System <- rep("2-18",nrow(mean_ctm_2_18))
mean_ctm_2_18$Rear_temp <- rep("18",nrow(mean_ctm_2_18))
mean_ctm_2_18$Rep <- rep("2",nrow(mean_ctm_2_18))
mean_ctm_2_18$animal <- paste(mean_ctm_2_18$Family, mean_ctm_2_18$Rep, sep = "_")

mean_ctm_rep2 <- merge(mean_ctm_2_15, mean_ctm_2_18, by="animal")
mean_ctm_rep2 <- mean_ctm_rep2[c("animal", "Family.x", "CTM_15", "CTM_18", "Rep.x")]

mean_ctm2 <- rbind(mean_ctm_rep1, mean_ctm_rep2)
colnames(mean_ctm2) <- c("animal", "Family", "CTM_15", "CTM_18", "Rep")

Fish_DI <- read.delim2("Fish_DI.txt", header=T)
ped_fam <- Fish_DI[c("Family", "Dam", "Sire")]
mean_ctm <- merge(mean_ctm2, ped_fam, by="Family")
mean_ctm_noNA <- na.omit(mean_ctm)



prior1.1<-list(G=list(G1=list(V=diag(2)*(0.002/1.002),n=1.002), #number in diag(X) is the number of response variables (ie temperature is 2 because its either 15 or 18C)
                      G2=list(V=diag(2)*(0.002/1.002),n=1.002),
                      G3=list(V=diag(2)*(0.002/1.002),n=1.002),
                      G4=list(V=diag(2)*(0.002/1.002),n=1.002)),
               R=list(V=diag(2)*(0.002/1.002),n=1.002)) #the number in this R list depends on how many variables in rcov

cov_mcmc.model <- MCMCglmm(fixed = cbind(CTM_15,CTM_18) ~  trait, 
                           random=~us(trait):animal + idh(trait):Dam + idh(trait):Rep +idh(trait):Family, 
                           nitt=1300000, thin=50, burnin=300000,
                           rcov = ~us(trait):units,
                           pedigree=Ped, data=mean_ctm_noNA, prior=prior1.1, 
                           family=c("gaussian", "gaussian"), verbose=TRUE)

save(cov_mcmc.model, file = "DS_hert_CoVmodels_.RData")
load("DS_hert_CoVmodels_1.RData")
load("DS_hert_CoVmodels_2.RData")

cov_mcmc.model <- cov_mcmc.model1
cov_mcmc.model <- cov_mcmc.model2

posterior.mode(cov_mcmc.model$VCV)
HPDinterval(cov_mcmc.model$VCV[,"traitCTM_15:traitCTM_18.animal"],0.95)
HPDinterval(cov_mcmc.model$VCV[,"traitCTM_15:traitCTM_15.animal"],0.95)
HPDinterval(cov_mcmc.model$VCV[,"traitCTM_18:traitCTM_18.animal"],0.95)
