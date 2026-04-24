#Helpful links:
#https://wam-tutorial.github.io/en/  ##check out this one first, definitely most helpful and what most of my code is based off
#https://github.com/tmalsburg/MCMCglmm-intro
#https://github.com/matthewwolak/wolakR

library(MCMCglmm)
library(dplyr)
library(ggplot2)
library("cowplot")

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign")
setwd("~/UCDavis/Smelt_heritability")

##############
# Data set up
##############


ctm <- read.table("2021_CTM_data_trial1-120.txt", header = T)
ctm$CTM <- as.numeric(ctm$CTM)
ctm$FL <- as.numeric(ctm$FL)
ctm$DI<- as.factor(ctm$DI)
ctm$Rear_temp<- as.factor(ctm$Rear_temp)

pedigree <- read.table("all_pedigree_meta_DI.txt", header=T)

ctm_pedigree <- merge(pedigree, ctm, by="Fish_ID", all.y = T)

ctm_pedigree <- na.omit(ctm_pedigree)
ctm_pedigree$animal <- ctm_pedigree$Fish_ID
ctm_pedigree$assigned_DI <- as.factor(ctm_pedigree$assigned_DI)
ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)

## start pedigree, manually add in parents, then read in final pedigree file. Run only once
#ped <- ctm_pedigree[,c(1,7,9)]
#write.table(ped, file = "MCMCglmm_Pedigree.txt", sep="\t", row.names = F)

ped <- read.table("MCMCglmm_Pedigree_withParents.txt", header=T, stringsAsFactors = FALSE)
ped <- ped[,c(1,3,2)] #I think dam needs to be second column

#############################
##check for missing parents
missing_dams <- data.frame(ped$dam[!(ped$dam %in% ped$animal)])
missing_dams <- unique(missing_dams)

missing_sires <- data.frame(ped$sire[!(ped$sire %in% ped$animal)])
missing_sires <- unique(missing_sires)

##manually add missing dams and sires to pedigree to the top of the pedigree!
##find duplicates in pedigree:
n_occur <- data.frame(table(ped$animal))
###########################


for (x in 1:3) ped[, x] <- as.factor(ped[, x]) #makes all columns factors
ped <- as.matrix(ped)


############
#PRIOR SETUP
############

##choose one!!
p.var<-var(ctm_pedigree$CTM,na.rm=TRUE)
p.var<-var(ctm_pedigree$FL,na.rm=TRUE)

prior.ctm.G1 <- list(G=list(G1=list(V=matrix(p.var/2), nu=1)),
                        R=list(V=matrix(p.var/2), nu=1))

prior.ctm.G2 <- list(G=list(G1=list(V=matrix(p.var/3), nu=1),
                               G2=list(V=matrix(p.var/3), nu=1)),
                        R=list(V=matrix(p.var/3), nu=1))

prior.ctm.G3 <- list(
  G = list(
    G1 = list(V = diag(2) * (p.var / 3), nu = 2)  # For idh(Rear_temp):animal
  ),
  R = list(
    V = diag(2) * (p.var / 3),  # For idh(Rear_temp):units
    nu = 2
  )
)


prior.ctm.G3.5 <- list(
  G = list(
    G1 = list(V = diag(2) * (p.var / 4), nu = 2),  # idh(Rear_temp):animal → 2 levels
    G2 = list(V = diag(3) * (p.var / 4), nu = 3)   # idh(assigned_DI):animal → 3 levels
  ),
  R = list(
    R1 = list(V = diag(2) * (p.var / 4), nu = 2),  # idh(Rear_temp):units → 2 levels
    R2 = list(V = diag(3) * (p.var / 4), nu = 3)   # idh(assigned_DI):units → 3 levels
  )
)




prior.ctm.G4 <- list(
  G = list(
    G1 = list(V = diag(2) * (p.var / 5), nu = 2),  # idh(Rear_temp):animal → 2 levels
    G2 = list(V = diag(3) * (p.var / 5), nu = 3),  # idh(assigned_DI):animal → 3 levels
    G3 = list(V = matrix(p.var / 5), nu = 1)       # dam (scalar)
  ),
  R = list(
    R1 = list(V = diag(2) * (p.var / 5), nu = 2),  # idh(Rear_temp):units → 2 levels
    R2 = list(V = diag(3) * (p.var / 5), nu = 3)   # idh(assigned_DI):units → 3 levels
  )
)

prior.ctm.G5 <- list(
  G = list(
    G1 = list(V = diag(2) * (p.var / 5), nu = 2),  # idh(Rear_temp):animal → 2 levels
    G2 = list(V = matrix(p.var / 5), nu = 1)       # dam (scalar)
  ),
  R = list(
    R1 = list(V = diag(2) * (p.var / 5), nu = 2)  # idh(Rear_temp):units → 2 levels
  )
)

prior.cov <- list(
  G = list(G1 = list(V = diag(2), nu = 1.002)),
  R = list(V = diag(2), nu = 1.002)
)

#this prior was too weak to run covmodel2
prior.cov2 <- list(
  G = list(
    G1 = list(V = diag(2), nu = 1.002),
    G2 = list(V = diag(2), nu = 1.002)),
  R = list(V = diag(2), nu = 1.002)
)


# compute phenotypic covariance for scaling
pheno <- ctm_pedigree[, c("CTM","FL")]
pheno <- ctm_pedigree[, c("CTM","Offspring_DI")]

S <- cov(pheno, use = "pairwise.complete.obs")
S
# scale V to S (here we use S / 2 as a reasonable start)
Vscale <- S / 2

prior.cov2 <- list(
  G = list(
    G1 = list(V = Vscale, nu = 3),   # animal
    G2 = list(V = Vscale, nu = 3) # Rear_temp
  ),
  R = list(V = Vscale, nu = 3)
)

prior.cov3 <- list(
  G = list(
    G1 = list(V = Vscale, nu = 3),   # for at.level(15):trait:animal
    G2 = list(V = Vscale, nu = 3)    # for at.level(18):trait:animal
  ),
  R = list(
    R1 = list(V = Vscale, nu = 3),   # for at.level(15):trait residual
    R2 = list(V = Vscale, nu = 3)    # for at.level(18):trait residual
  )
)


#################
## Models CTM
################


model1 <- MCMCglmm(CTM ~ System,
                   random = ~ idh(Rear_temp):animal,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model1, file="model1.RData")

model1.5 <- MCMCglmm(CTM ~ System + FL,
                   random = ~ idh(Rear_temp):animal,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model1.5, file="model1.FL.RData")


model2 <- MCMCglmm(CTM ~ System,
                   random = ~ idh(Rear_temp):animal + idh(assigned_DI):animal,
                   rcov = ~ idh(Rear_temp):units + idh(assigned_DI):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G3.5, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model2, file="model2.RData")


model2.5 <- MCMCglmm(CTM ~ System,
                   random = ~ idh(Rear_temp):animal + Dam.1,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G5, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model2.5, file="model2.5.RData")



model3 <- MCMCglmm(CTM ~ System,
                      random = ~ idh(Rear_temp):animal + idh(DI):animal + dam,
                      rcov = ~ idh(Rear_temp):units + idh(DI):units,
                      family = "gaussian",
                      pedigree = ped, data = ctm_pedigree,
                      nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                      prior = prior.ctm.G4, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model3, file="model3.RData")



Low <- subset(ctm_pedigree, assigned_DI=="L")
model4.L <- MCMCglmm(CTM ~ System,
                   random = ~ idh(Rear_temp):animal,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = Low,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.L, file="model4.L.RData")


Med <- subset(ctm_pedigree, assigned_DI=="M")
model4.M <- MCMCglmm(CTM ~ System,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = Med,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.M, file="model4.M.RData")


High <- subset(ctm_pedigree, assigned_DI=="H")
model4.H <- MCMCglmm(CTM ~ System,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = High,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.H, file="model4.H.RData")

#### with fl in models
Low <- subset(ctm_pedigree, assigned_DI=="L")
model4.L <- MCMCglmm(CTM ~ System + FL,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = Low,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.L, file="model4.L.FL.RData")


Med <- subset(ctm_pedigree, assigned_DI=="M")
model4.M <- MCMCglmm(CTM ~ System +FL,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = Med,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.M, file="model4.M.FL.RData")


High <- subset(ctm_pedigree, assigned_DI=="H")
model4.H <- MCMCglmm(CTM ~ System + FL,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = High,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.H, file="model4.H.FL.RData")


###############
## Covar model
###############
cov_mcmc.model1 <- MCMCglmm(fixed = cbind(CTM,FL) ~ trait -1 + trait:System, 
                            random=~us(trait):animal, 
                            rcov = ~us(trait):units,
                            nitt = 500000, burnin = 1000, thin = 500,
                            pedigree=ped, data=ctm_pedigree, prior=prior.cov, 
                            family=c("gaussian", "gaussian"), verbose=T)
save(cov_mcmc.model1, file = "DS_hert_CoVmodels1_ctm_fl.RData")


cov_mcmc.model2 <- MCMCglmm(fixed = cbind(CTM,FL) ~ trait -1 + trait:System, 
                            random=~us(trait):animal + us(trait):Rear_temp, 
                            rcov = ~us(trait):units,
                            nitt = 500000, burnin = 1000, thin = 500,
                            pedigree=ped, data=ctm_pedigree, prior=prior.cov2, 
                            family=c("gaussian", "gaussian"), verbose=T)
save(cov_mcmc.model2, file = "DS_hert_CoVmodels2_ctm_fl.RData")


cov_mcmc.model3 <- MCMCglmm(fixed = cbind(CTM,FL) ~ trait -1 + trait:System, 
                            random= ~ us(at.level(Rear_temp, "15"):trait):animal + us(at.level(Rear_temp, "18"):trait):animal, 
                            rcov = ~ us(at.level(Rear_temp, "15"):trait):units + us(at.level(Rear_temp, "18"):trait):units,
                            nitt = 500000, burnin = 1000, thin = 500,
                            pedigree=ped, data=ctm_pedigree, prior=prior.cov3, 
                            family=c("gaussian", "gaussian"), verbose=T)
save(cov_mcmc.model3, file = "DS_hert_CoVmodels3_ctm_fl.RData")


cov_mcmc.model4 <- MCMCglmm(fixed = cbind(CTM, Offspring_DI) ~ trait -1 + trait:System, 
                            random= ~ us(at.level(Rear_temp, "15"):trait):animal + us(at.level(Rear_temp, "18"):trait):animal, 
                            rcov = ~ us(at.level(Rear_temp, "15"):trait):units + us(at.level(Rear_temp, "18"):trait):units,
                            nitt = 500000, burnin = 1000, thin = 500,
                            pedigree=ped, data=ctm_pedigree, prior=prior.cov3, 
                            family=c("gaussian", "gaussian"), verbose=T)
save(cov_mcmc.model4, file = "DS_hert_CoVmodels4_ctm_DI.RData")

##############
# VALIDITY CHECKS
###############

load("model1.RData") #final model 
load("model1.FL.RData")
load("DS_hert_CoVmodels1_ctm_fl.RData")
load("DS_hert_CoVmodels2_ctm_fl.RData")
load("DS_hert_CoVmodels3_ctm_fl.RData")
load("DS_hert_CoVmodels4_ctm_DI.RData")
load("model4.L.RData")
load("model4.M.RData")
load("model4.H.RData")



# Assess significance of fixed and random effects
#check fixed effects interval doesn't cross 0
posterior.mode(model1$Sol)
HPDinterval(model1$Sol,0.95)

posterior.mode(model1.5$Sol)
HPDinterval(model1.5$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(model2$Sol)
HPDinterval(model2$Sol,0.95)

posterior.mode(model3$Sol)
HPDinterval(model3$Sol,0.95)

posterior.mode(model4.L$Sol)
HPDinterval(model4.L$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(model4.M$Sol)
HPDinterval(model4.M$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(model4.H$Sol)
HPDinterval(model4.H$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(cov_mcmc.model1$Sol)
HPDinterval(cov_mcmc.model1$Sol,0.95)

posterior.mode(cov_mcmc.model2$Sol)
HPDinterval(cov_mcmc.model2$Sol,0.95)

posterior.mode(cov_mcmc.model3$Sol)
HPDinterval(cov_mcmc.model3$Sol,0.95)

posterior.mode(cov_mcmc.model4$Sol)
HPDinterval(cov_mcmc.model4$Sol,0.95)


model1$DIC #11367.52
model1.5$DIC #11369.86
model2$DIC
model2.5$DIC #11400.87
model3$DIC
model4.L$DIC
model4.M$DIC
model4.H$DIC

# autocorrelation plot
autocorr(model1$Sol)
autocorr(model1$VCV)
plot(model1$VCV)

autocorr(model2$Sol)
autocorr(model2$VCV)
plot(model2$VCV)

autocorr(model3$Sol)
autocorr(model3$VCV)
plot(model3$VCV)

autocorr(model4.L$Sol)
autocorr(model4.L$VCV)
plot(model4.L$VCV)

autocorr(model4.M$Sol)
autocorr(model4.M$VCV)
plot(model4.M$VCV)

autocorr(model4.H$Sol)
autocorr(model4.H$VCV)
plot(model4.H$VCV)

autocorr(cov_mcmc.model1$Sol)
autocorr(cov_mcmc.model1$VCV)
plot(cov_mcmc.model1$VCV)

autocorr(cov_mcmc.model2$Sol)
autocorr(cov_mcmc.model2$VCV)
plot(cov_mcmc.model2$VCV)

autocorr(cov_mcmc.model3$Sol)
autocorr(cov_mcmc.model3$VCV)
plot(cov_mcmc.model3$VCV)

autocorr(cov_mcmc.model4$Sol)
autocorr(cov_mcmc.model4$VCV)
plot(cov_mcmc.model4$VCV)

#################
# Phenotypic and genetic variances
################

#Va and other variances
posterior.mode(model1$VCV)
HPDinterval(model1$VCV,0.95)
VP_posterior <- rowSums(model1$VCV) # Calculate the posterior distribution of VP:
mean(VP_posterior) # Get the point estimate (Posterior Mean) of VP:
VP_mcmc_object <- as.mcmc(VP_posterior)
HPDinterval(VP_mcmc_object,0.95)

posterior.mode(model1.5$VCV)
HPDinterval(model1.5$VCV,0.95)

posterior.mode(model2$VCV)
HPDinterval(model2$VCV,0.95)

posterior.mode(model3$VCV)
HPDinterval(model3$VCV,0.95)

posterior.mode(model4.L$VCV)
HPDinterval(model4.L$VCV,0.95)

posterior.mode(model4.M$VCV)
HPDinterval(model4.M$VCV,0.95)

posterior.mode(model4.H$VCV)
HPDinterval(model4.H$VCV,0.95)

posterior.mode(cov_mcmc.model1$VCV)
HPDinterval(cov_mcmc.model1$VCV,0.95)

posterior.mode(cov_mcmc.model2$VCV)
HPDinterval(cov_mcmc.model2$VCV,0.95)

posterior.mode(cov_mcmc.model3$VCV)
HPDinterval(cov_mcmc.model3$VCV,0.95)

posterior.mode(cov_mcmc.model4$VCV)
HPDinterval(cov_mcmc.model4$VCV,0.95)


#################
##Heritability
#################

#################### model 1
posterior.heritability1.18 <- model1$VCV[, "Rear_temp18.animal"] /
  (model1$VCV[, "Rear_temp18.animal"] + model1$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability1.18)
HPDinterval(posterior.heritability1.18, 0.95)


posterior.heritability1.15 <- model1$VCV[, "Rear_temp15.animal"] /
  (model1$VCV[, "Rear_temp15.animal"] + model1$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability1.15)
HPDinterval(posterior.heritability1.15, 0.95)


# Mean for the reference group (Intercept)
mu_S1_15 <- model1$Sol[, "(Intercept)"]

# Means for the other groups (Intercept + Contrast)
mu_S1_18 <- model1$Sol[, "(Intercept)"] + model1$Sol[, "System1-18"]
mu_S2_15 <- model1$Sol[, "(Intercept)"] + model1$Sol[, "System2-15"]
mu_S2_18 <- model1$Sol[, "(Intercept)"] + model1$Sol[, "System2-18"]

# System 1-15
IA_S1_15 <- model1$VCV[, "Rear_temp15.animal"] / (mu_S1_15^2)
posterior.mode(IA_S1_15)
HPDinterval(IA_S1_15, 0.95)

#  System 1-18
IA_S1_18 <- model1$VCV[, "Rear_temp18.animal"] / (mu_S1_18^2)
posterior.mode(IA_S1_18)
HPDinterval(IA_S1_18, 0.95)

# System 2-15
IA_S2_15 <- model1$VCV[, "Rear_temp15.animal"] / (mu_S2_15^2)
posterior.mode(IA_S2_15)
HPDinterval(IA_S2_15, 0.95)

#  System 2-18
IA_S2_18 <- model1$VCV[, "Rear_temp18.animal"] / (mu_S2_18^2)
posterior.mode(IA_S2_18)
HPDinterval(IA_S2_18, 0.95)




# Combine them (Simple Average of the posteriors)
mu_15_combined <- (mu_S1_15 + mu_S2_15) / 2
# Calculate combined Evolvability using the 15 degree Variance
IA_1_combined <- model1$VCV[, "Rear_temp15.animal"] / (mu_15_combined^2)
# Combine them (Simple Average of the posteriors)
mu_18_combined <- (mu_S1_18 + mu_S2_18) / 2
# Calculate combined Evolvability using the 15 degree Variance
IA_18_combined <- model1$VCV[, "Rear_temp18.animal"] / (mu_18_combined^2)

# System 2-15
IA_15 <- model1$VCV[, "Rear_temp15.animal"] / (mu_15_combined^2)
posterior.mode(IA_15)
HPDinterval(IA_15, 0.95)

#  System 2-18
IA_18 <- model1$VCV[, "Rear_temp18.animal"] / (mu_18_combined^2)
posterior.mode(IA_18)
HPDinterval(IA_18, 0.95)






#################### model 1.5
posterior.heritability1.18 <- model1.5$VCV[, "Rear_temp18.animal"] /
  (model1.5$VCV[, "Rear_temp18.animal"] + model1.5$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability1.18)
HPDinterval(posterior.heritability1.18, 0.95)


posterior.heritability1.15 <- model1.5$VCV[, "Rear_temp15.animal"] /
  (model1.5$VCV[, "Rear_temp15.animal"] + model1.5$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability1.15)
HPDinterval(posterior.heritability1.15, 0.95)



#################### model2
posterior.heritability2.18 <- model2$VCV[, "Rear_temp18.animal"] /
  (model2$VCV[, "Rear_temp18.animal"] + model2$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability2.18)
HPDinterval(posterior.heritability2.18, 0.95)


posterior.heritability2.15 <- model2$VCV[, "Rear_temp15.animal"] /
  (model2$VCV[, "Rear_temp15.animal"] + model2$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability2.15)
HPDinterval(posterior.heritability2.15, 0.95)




posterior.heritability2.L <- model2$VCV[, "assigned_DIL.animal"] /
  (model2$VCV[, "assigned_DIL.animal"] + model2$VCV[, "assigned_DIL.units"])

posterior.mode(posterior.heritability2.L)
HPDinterval(posterior.heritability2.L, 0.95)


posterior.heritability2.M <- model2$VCV[, "assigned_DIM.animal"] /
  (model2$VCV[, "assigned_DIM.animal"] + model2$VCV[, "assigned_DIM.units"])

posterior.mode(posterior.heritability2.M)
HPDinterval(posterior.heritability2.M, 0.95)



posterior.heritability2.H <- model2$VCV[, "assigned_DIH.animal"] /
  (model2$VCV[, "assigned_DIH.animal"] + model2$VCV[, "assigned_DIH.units"])

posterior.mode(posterior.heritability2.H)
HPDinterval(posterior.heritability2.H, 0.95)


#################### model3
posterior.heritability3.18 <- model3$VCV[, "Rear_temp18.animal"] /
  (model3$VCV[, "Rear_temp18.animal"] + model3$VCV[, "Rear_temp18.units"] + model3$VCV[,"Dam"])

posterior.mode(posterior.heritability3.18)
HPDinterval(posterior.heritability3.18, 0.95)


posterior.heritability3.15 <- model3$VCV[, "Rear_temp15.animal"] /
  (model3$VCV[, "Rear_temp15.animal"] + model3$VCV[, "Rear_temp15.units"] + model3$VCV[,"Dam"])

posterior.mode(posterior.heritability3.15)
HPDinterval(posterior.heritability3.15, 0.95)




posterior.heritability3.L <- model3$VCV[, "assigned_DIL.animal"] /
  (model3$VCV[, "assigned_DIL.animal"] + model3$VCV[, "assigned_DIL.units"] + model3$VCV[,"Dam"])

posterior.mode(posterior.heritability3.L)
HPDinterval(posterior.heritability3.L, 0.95)


posterior.heritability3.M <- model3$VCV[, "assigned_DIM.animal"] /
  (model3$VCV[, "assigned_DIM.animal"] + model3$VCV[, "assigned_DIM.units"] + model3$VCV[,"Dam"])

posterior.mode(posterior.heritability3.M)
HPDinterval(posterior.heritability3.M, 0.95)



posterior.heritability3.H <- model3$VCV[, "assigned_DIH.animal"] /
  (model3$VCV[, "assigned_DIH.animal"] + model3$VCV[, "assigned_DIH.units"] + model3$VCV[,"Dam"])

posterior.mode(posterior.heritability3.H)
HPDinterval(posterior.heritability3.H, 0.95)


#################### model4

##LOW
posterior.heritability4.L.18 <- model4.L$VCV[, "Rear_temp18.animal"] /
  (model4.L$VCV[, "Rear_temp18.animal"] + model4.L$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.L.18)
HPDinterval(posterior.heritability4.L.18, 0.95)


posterior.heritability4.L.15 <- model4.L$VCV[, "Rear_temp15.animal"] /
  (model4.L$VCV[, "Rear_temp15.animal"] + model4.L$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.L.15)
HPDinterval(posterior.heritability4.L.15, 0.95)


## MED
posterior.heritability4.M.18 <- model4.M$VCV[, "Rear_temp18.animal"] /
  (model4.M$VCV[, "Rear_temp18.animal"] + model4.M$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.M.18)
HPDinterval(posterior.heritability4.M.18, 0.95)


posterior.heritability4.M.15 <- model4.M$VCV[, "Rear_temp15.animal"] /
  (model4.M$VCV[, "Rear_temp15.animal"] + model4.M$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.M.15)
HPDinterval(posterior.heritability4.M.15, 0.95)


## HIGH
posterior.heritability4.H.18 <- model4.H$VCV[, "Rear_temp18.animal"] /
  (model4.H$VCV[, "Rear_temp18.animal"] + model4.H$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.H.18)
HPDinterval(posterior.heritability4.H.18, 0.95)


posterior.heritability4.H.15 <- model4.H$VCV[, "Rear_temp15.animal"] /
  (model4.H$VCV[, "Rear_temp15.animal"] + model4.H$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.H.15)
HPDinterval(posterior.heritability4.H.15, 0.95)




################## genetic correlation
genetic.correlation <- cov_mcmc.model1$VCV[, "traitCTM:traitFL.animal"] / sqrt(cov_mcmc.model1$VCV[, "traitCTM:traitCTM.animal"] * cov_mcmc.model1$VCV[, "traitFL:traitFL.animal"])
posterior.mode(genetic.correlation)
corr_interval <- HPDinterval(genetic.correlation,0.95)


genetic.correlation <- cov_mcmc.model2$VCV[, "traitCTM:traitFL.animal"] / sqrt(cov_mcmc.model2$VCV[, "traitCTM:traitCTM.animal"] * cov_mcmc.model2$VCV[, "traitFL:traitFL.animal"])
posterior.mode(genetic.correlation)
corr_interval <- HPDinterval(genetic.correlation,0.95)
corr_interval




num15 <- cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"15\"):traitCTM:at.level(Rear_temp, \"15\"):traitFL.animal"]
den15 <- sqrt(
  cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"15\"):traitCTM:at.level(Rear_temp, \"15\"):traitCTM.animal"] *
    cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"15\"):traitFL:at.level(Rear_temp, \"15\"):traitFL.animal"]
)
genetic.correlation15 <- num15 / den15
posterior.mode(genetic.correlation15)
HPDinterval(genetic.correlation15)

num18 <- cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"18\"):traitCTM:at.level(Rear_temp, \"18\"):traitFL.animal"]
den18 <- sqrt(
  cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"18\"):traitCTM:at.level(Rear_temp, \"18\"):traitCTM.animal"] *
    cov_mcmc.model3$VCV[, "at.level(Rear_temp, \"18\"):traitFL:at.level(Rear_temp, \"18\"):traitFL.animal"]
)
genetic.correlation18 <- num18 / den18
posterior.mode(genetic.correlation18)
HPDinterval(genetic.correlation18)




num15 <- cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"15\"):traitCTM:at.level(Rear_temp, \"15\"):traitOffspring_DI.animal"]
den15 <- sqrt(
  cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"15\"):traitCTM:at.level(Rear_temp, \"15\"):traitCTM.animal"] *
    cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"15\"):traitOffspring_DI:at.level(Rear_temp, \"15\"):traitOffspring_DI.animal"]
)
genetic.correlation15 <- num15 / den15
posterior.mode(genetic.correlation15)
HPDinterval(genetic.correlation15)

num18 <- cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"18\"):traitCTM:at.level(Rear_temp, \"18\"):traitOffspring_DI.animal"]
den18 <- sqrt(
  cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"18\"):traitCTM:at.level(Rear_temp, \"18\"):traitCTM.animal"] *
    cov_mcmc.model4$VCV[, "at.level(Rear_temp, \"18\"):traitOffspring_DI:at.level(Rear_temp, \"18\"):traitOffspring_DI.animal"]
)
genetic.correlation18 <- num18 / den18
posterior.mode(genetic.correlation18)
HPDinterval(genetic.correlation18)

#####################
## Figues
###################

##model without fl
temp <- factor(c("15", "18", "15", "18", "15", "18", "15", "18"), levels = c("15", "18"))
DI <- factor(c("All", "All", "Low", "Low", "Medium", "Medium", "High", "High"), levels=c("All", "Low","Medium", "High"))
va <- c("1.91", "0.53", "1.8", "0.64", "1.46", "0.87", "1.32", "0.31")
lower <- c("1.02", "0.29", "0.39", "0.19", "0.67", "0.36", "0.32", "0.17")
upper <- c("2.70", "0.87", "6.41", "1.69", "3.04", "1.67", "2.98", "0.81")



# Combine vectors into a data frame
ad_gen_var <- data.frame(va, lower, upper, temp, DI)
ad_gen_var$va <- as.numeric(ad_gen_var$va)
ad_gen_var$lower <- as.numeric(ad_gen_var$lower)
ad_gen_var$upper <- as.numeric(ad_gen_var$upper)

windows()
plot1 <- ggplot(data=ad_gen_var, aes(x=DI, y=va, color= temp, group = temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=ad_gen_var, aes(ymin=lower, ymax=upper), width=0, linewidth = 1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Additive Genetic Variation"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 7, by = 1))



## heritability, model without fl
temp <- c("15", "18", "15", "18", "15", "18", "15", "18")
DI <- factor(c("All", "All", "Low", "Low", "Medium", "Medium", "High", "High"), levels=c("All", "Low","Medium", "High"))
h2 <- c(0.2594873, 0.1561065 , 0.1855789, 0.1479088, 0.2431994, 0.2732441, 0.2672351, 0.1307828 )
lower <- c(0.1667586 , 0.1086305 , 0.06420441, 0.05785591, 0.1279629, 0.155862, 0.1003419, 0.07240127)
upper <- c(0.397811, 0.2962406, 0.6137362, 0.4282, 0.4589953, 0.5352722, 0.498974, 0.3005448)



# Combine vectors into a data frame
heritability <- data.frame(h2, lower, upper, temp, DI)

windows()
plot2 <- ggplot(data=heritability, aes(x=DI, y=h2, color= temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=heritability, aes(ymin=lower, ymax=upper), width=0, linewidth = 1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Proportion of variation"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) 
coord_cartesian(ylim = 0:0.5)

## Evolvability, model without fl
temp <- c("15", "18", "15", "18", "15", "18", "15", "18")
DI <- factor(c("All", "All", "L", "L", "M", "M", "H", "H"), levels=c("All", "L","M", "H"))
h2 <- c(0.2596702, 0.06512294, 0.1672319, 0.08204916, 0.2266543, 0.1095768, 0.1515815, 0.03811622)
lower <- c(0.1369523, 0.03812755, 0.05412295, 0.02450536, 0.08839724, 0.04628809, 0.04602467, 0.02071461)
upper <- c(0.3556641, 0.1108644, 0.8647676, 0.2170645, 0.4065484, 0.2106581, 0.3749894, 0.09912414)



# Combine vectors into a data frame
heritability <- data.frame(h2, lower, upper, temp, DI)

windows()
plot3 <- ggplot(data=heritability, aes(x=DI, y=h2, color= temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=heritability, aes(ymin=lower, ymax=upper), width=0, linewidth = 1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Coeff. of Genetic Variation (%)"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.position = "none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) 
coord_cartesian(ylim = 0:0.5)


#################
## Models FL
################


model1.FL.h2 <- MCMCglmm(FL ~ System,
                   random = ~ idh(Rear_temp):animal,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model1.FL.h2, file="model1.FL.h2.RData")

model3.FL.h2 <- MCMCglmm(FL ~ System,
                   random = ~ idh(Rear_temp):animal + Dam.1,
                   rcov = ~ idh(Rear_temp):units,
                   family = "gaussian",
                   pedigree = ped, data = ctm_pedigree,
                   nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                   prior = prior.ctm.G5, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model3.FL.h2, file="model3.FL.h2.RData")


Low <- subset(ctm_pedigree, assigned_DI=="L")
model4.L.FL.h2 <- MCMCglmm(FL ~ System,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = Low,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.L.FL.h2, file="model4.L.FL.h2.RData")


Med <- subset(ctm_pedigree, assigned_DI=="M")
model4.M.FL.h2 <- MCMCglmm(FL ~ System,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = Med,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.M.FL.h2, file="model4.M.FL.h2.RData")


High <- subset(ctm_pedigree, assigned_DI=="H")
model4.H.FL.h2 <- MCMCglmm(FL ~ System,
                     random = ~ idh(Rear_temp):animal,
                     rcov = ~ idh(Rear_temp):units,
                     family = "gaussian",
                     pedigree = ped, data = High,
                     nitt = 500000, burnin = 1000, thin = 500, ##for actual results
                     prior = prior.ctm.G3, verbose = TRUE, slice = FALSE, DIC = TRUE)
save(model4.H.FL.h2, file="model4.H.FL.h2.RData")


##############
# VALIDITY CHECKS
###############

load("model1.FL.h2.RData")

load("model3.FL.h2.RData")
load("model4.L.FL.h2.RData")
load("model4.M.FL.h2.RData")
load("model4.H.FL.h2.RData")

model1 <- model1.FL.h2
model3 <- model3.FL.h2
model4.L <- model4.L.FL.h2
model4.M <- model4.M.FL.h2
model4.H <- model4.H.FL.h2

# Assess significance of fixed and random effects
#check fixed effects interval doesn't cross 0
posterior.mode(model1$Sol)
HPDinterval(model1$Sol,0.95)
            
posterior.mode(model3$Sol)
HPDinterval(model3$Sol,0.95)

posterior.mode(model4.L$Sol)
HPDinterval(model4.L$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(model4.M$Sol)
HPDinterval(model4.M$Sol,0.95) ##FL crosses 0 so not useful addition to the model

posterior.mode(model4.H$Sol)
HPDinterval(model4.H$Sol,0.95) ##FL crosses 0 so not useful addition to the model

model1$DIC #11080.52
model3$DIC #11299
model4.L$DIC
model4.M$DIC
model4.H$DIC

# autocorrelation plot
autocorr(model1$Sol)
autocorr(model1$VCV)
plot(model1$VCV)

autocorr(model3.FL.h2$Sol)
autocorr(model3.FL.h2$VCV)
plot(model3$VCV)

autocorr(model4.L$Sol)
autocorr(model4.L$VCV)
plot(model4.L$VCV)

autocorr(model4.M$Sol)
autocorr(model4.M$VCV)
plot(model4.M$VCV)

autocorr(model4.H$Sol)
autocorr(model4.H$VCV)
plot(model4.H$VCV)

#################
# Phenotypic and genetic variances
################

#Va and other variances
posterior.mode(model1$VCV)
HPDinterval(model1$VCV,0.95)
VP_posterior <- rowSums(model1$VCV) # Calculate the posterior distribution of VP:
mean(VP_posterior) # Get the point estimate (Posterior Mean) of VP:
VP_mcmc_object <- as.mcmc(VP_posterior)
HPDinterval(VP_mcmc_object,0.95)

posterior.mode(model3.FL.h2$VCV)
HPDinterval(model3.FL.h2$VCV,0.95)

posterior.mode(model4.L.FL.h2$VCV)
HPDinterval(model4.L.FL.h2$VCV,0.95)

posterior.mode(model4.M.FL.h2$VCV)
HPDinterval(model4.M.FL.h2$VCV,0.95)

posterior.mode(model4.H.FL.h2$VCV)
HPDinterval(model4.H.FL.h2$VCV,0.95)


#################
##Heritability
#################

#################### model 1
posterior.heritability1.18 <- model1$VCV[, "Rear_temp18.animal"] /
  (model1$VCV[, "Rear_temp18.animal"] + model1$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability1.18)
HPDinterval(posterior.heritability1.18, 0.95)


posterior.heritability1.15 <- model1$VCV[, "Rear_temp15.animal"] /
  (model1$VCV[, "Rear_temp15.animal"] + model1$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability1.15)
HPDinterval(posterior.heritability1.15, 0.95)



#################### model 3
posterior.heritability3.18 <- model3$VCV[, "Rear_temp18.animal"] /
  (model3$VCV[, "Rear_temp18.animal"] + model3$VCV[, "Rear_temp18.units"] + model3$VCV[, "Dam.1"])

posterior.mode(posterior.heritability3.18)
HPDinterval(posterior.heritability3.18, 0.95)


posterior.heritability3.15 <- model3$VCV[, "Rear_temp15.animal"] /
  (model3$VCV[, "Rear_temp15.animal"] + model3$VCV[, "Rear_temp15.units"] + model3$VCV[, "Dam.1"])

posterior.mode(posterior.heritability3.15)
HPDinterval(posterior.heritability3.15, 0.95)




#################### model4

##LOW
posterior.heritability4.L.18 <- model4.L$VCV[, "Rear_temp18.animal"] /
  (model4.L$VCV[, "Rear_temp18.animal"] + model4.L$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.L.18)
HPDinterval(posterior.heritability4.L.18, 0.95)


posterior.heritability4.L.15 <- model4.L$VCV[, "Rear_temp15.animal"] /
  (model4.L$VCV[, "Rear_temp15.animal"] + model4.L$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.L.15)
HPDinterval(posterior.heritability4.L.15, 0.95)


## MED
posterior.heritability4.M.18 <- model4.M$VCV[, "Rear_temp18.animal"] /
  (model4.M$VCV[, "Rear_temp18.animal"] + model4.M$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.M.18)
HPDinterval(posterior.heritability4.M.18, 0.95)


posterior.heritability4.M.15 <- model4.M$VCV[, "Rear_temp15.animal"] /
  (model4.M$VCV[, "Rear_temp15.animal"] + model4.M$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.M.15)
HPDinterval(posterior.heritability4.M.15, 0.95)


## HIGH
posterior.heritability4.H.18 <- model4.H$VCV[, "Rear_temp18.animal"] /
  (model4.H$VCV[, "Rear_temp18.animal"] + model4.H$VCV[, "Rear_temp18.units"])

posterior.mode(posterior.heritability4.H.18)
HPDinterval(posterior.heritability4.H.18, 0.95)


posterior.heritability4.H.15 <- model4.H$VCV[, "Rear_temp15.animal"] /
  (model4.H$VCV[, "Rear_temp15.animal"] + model4.H$VCV[, "Rear_temp15.units"])

posterior.mode(posterior.heritability4.H.15)
HPDinterval(posterior.heritability4.H.15, 0.95)

##
temp <- c("15", "18","15", "18", "15", "18", "15", "18", "15", "18")
DI <- factor(c("All", "All", "All.Dam", "All.Dam", "Low", "Low", "Medium", "Medium", "High", "High"), levels=c("All", "All.Dam", "Low","Medium", "High"))
va <- c("2.22", "3.48", "1.23", "2.42", "1.34", "4.46", "3.39", "2.91", "1.34", "3.50")
lower <- c("1.47", "1.99", "0.80", "1.38", "0.64", "1.46", "1.93", "1.52", "0.68", "1.59")
upper <- c("3.34", "5.16", "2.26", "3.87", "2.30", "6.49", "4.78", "6.27", "3.01", "8.00")




# Combine vectors into a data frame
ad_gen_var2 <- data.frame(va, lower, upper, temp, DI)
ad_gen_var2$va <- as.numeric(ad_gen_var2$va)
ad_gen_var2$lower <- as.numeric(ad_gen_var2$lower)
ad_gen_var2$upper <- as.numeric(ad_gen_var2$upper)

windows()
plot4 <- ggplot(data=ad_gen_var2, aes(x=DI, y=va, color= temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=ad_gen_var2, aes(ymin=lower, ymax=upper), width=0, linewidth=1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Additive Genetic Variation"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), legend.position = "none") +
  scale_y_continuous(breaks = seq(0, 8, by = 1))



temp <- c("15", "18","15", "18", "15", "18", "15", "18", "15", "18")
DI <- factor(c("All", "All", "All.Dam", "All.Dam", "Low", "Low", "Medium", "Medium", "High", "High"), levels=c("All", "All.Dam", "Low","Medium", "High"))
h2 <- c(0.5809397, 0.4970823, 0.3134261, 0.3447556, 0.651319, 0.7162546, 0.6787631, 0.4288928, 0.4143741, 0.5669099)
lower <- c(0.4106893, 0.3278506, 0.1987495, 0.2175055, 0.3438686, 0.344096, 0.4662616, 0.2483182, 0.2329916, 0.2588781)
upper <- c(0.7448825, 0.6771461, 0.5111179, 0.5131906, 0.7996708, 0.8992364, 0.8681884, 0.7737474, 0.6945326, 0.8314059)




# Combine vectors into a data frame
heritability2 <- data.frame(h2, lower, upper, temp, DI)

windows()
plot5 <- ggplot(data=heritability2, aes(x=DI, y=h2, color= temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=heritability2, aes(ymin=lower, ymax=upper), width=0, linewidth=1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Proportion of variation"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), legend.position = "none") 
#coord_cartesian(ylim = 0:0.5)

## Evolvability
temp <- c("15", "18","15", "18", "15", "18", "15", "18", "15", "18")
DI <- factor(c("All", "All","All.Dam", "All.Dam", "L", "L", "M", "M", "H", "H"), levels=c("All", "All.Dam", "L","M", "H"))
h2 <- c(0.7816576, 0.7441232, 0.4858595, 0.6525496, 0.4260335, 0.9832069, 1.162496, 0.7575702, 0.5978339, 0.9429375)
lower <- c(0.4980565, 0.5729288, 0.2684277, 0.3714167, 0.234068, 0.4466539, 0.5993086, 0.3818246, 0.2679136, 0.4368071)
upper <- c(1.160729, 1.421189, 0.794027, 1.048118, 0.831465, 1.916819, 1.566284, 1.59959, 1.137799, 2.074198)



# Combine vectors into a data frame
heritability <- data.frame(h2, lower, upper, temp, DI)

windows()
plot6 <- ggplot(data=heritability, aes(x=DI, y=h2, color= temp)) +
  geom_point(stat="identity", position=position_dodge(0.5), size = 3) +
  geom_errorbar(data=heritability, aes(ymin=lower, ymax=upper), width=0, linewidth = 1, position=position_dodge(0.5)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  labs(y=expression("Coeff. of Genetic Variation (%)"), x="", fill="") +
  #geom_text(data = generate_label_df(posthoc1, 'lev'), aes(x = plot.labels, y = V1, label = labels)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), legend.position = "none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) 
coord_cartesian(ylim = 0:0.5)




windows()
plot_grid(plot1, plot2, plot3, plot4,plot5, plot6, labels = "AUTO")
