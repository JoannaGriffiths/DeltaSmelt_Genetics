
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library("multcomp")
library("bestNormalize")
library("outliers")
library("cowplot")
library("performance")

setwd("~/UCDavis/FCCL/Spawning_2021/CTMs")

ctm <- read.table("2021_CTM_data_trial1-120.txt", header = T)
ctm$CTM <- as.numeric(ctm$CTM)
ctm$FL <- as.numeric(ctm$FL)
ctm$DI<- as.factor(ctm$DI)
ctm$Rear_temp<- as.factor(ctm$Rear_temp)
shapiro.test(ctm$CTM) #p-value < 2.2e-16
plot(hist(ctm$CTM))
leveneTest(ctm$CTM, group = ctm$Rear_temp) #< 2.2e-16 ***
leveneTest(ctm$CTM, group = ctm$DI) #5.079e-10 ***


#bestNormalize(ctm$CTM) #suggests orderNorm is best transformation
ctm$orderNorm_CTM <- predict(orderNorm(ctm$CTM))
shapiro.test(ctm$orderNorm_CTM) #p-value = 0.7411
plot(hist(ctm$orderNorm_CTM))
leveneTest(ctm$orderNorm_CTM, group = ctm$Rear_temp) #< 2.2e-16 ***; significant, but also clear from data that variance is lower and maybe less genetic variation
leveneTest(ctm$orderNorm_CTM, group = ctm$DI) #0.0001012 ***

##normlize by FL?
ctm$FL_CTM <- ctm$CTM / ctm$FL
shapiro.test(ctm$FL_CTM) #p-value < 2.2e-16
ctm$orderNorm_FL_CTM <- predict(orderNorm(ctm$FL_CTM))
shapiro.test(ctm$orderNorm_FL_CTM) #p-value =1
plot(hist(ctm$orderNorm_FL_CTM))
leveneTest(ctm$orderNorm_FL_CTM, group = ctm$Rear_temp) #< 1.319e-06 ****; significant, but also clear from data that variance is lower and maybe less genetic variation
leveneTest(ctm$orderNorm_FL_CTM, group = ctm$DI) #2.123e-13 ***


#####normalize FL data
#omit outliers:
grubbs.test(ctm$FL, type = 10)
ctm<-subset(ctm, Fish_ID!="076S")
ctm<-subset(ctm, Fish_ID!="108M")



shapiro.test(ctm$FL) #p-value = 1.115e-12
plot(hist(ctm$FL))
leveneTest(ctm$FL, group = ctm$Rear_temp)
leveneTest(ctm$FL, group = ctm$DI)

bestNormalize(ctm$FL) #suggests center_scale(x) is best transformation
ctm$orderNorm_FL <- predict(orderNorm(ctm$FL))
ctm$center_FL <- predict(center_scale(ctm$FL))
plot(hist(ctm$center_FL))
shapiro.test(ctm$orderNorm_FL) #p-value = 8.854e-07
shapiro.test(ctm$center_FL) #p-value = 1.115e-12

##performs worse:
#ctm.FL <- ctm[!is.na(ctm$FL),] # remove NAs for fixed effects
#ctm.FL$sqrt_FL <- sqrt(ctm.FL$FL)
#ctm.FL$sqrt_FL <- as.numeric(ctm.FL$sqrt_FL)
#shapiro.test(ctm.FL$sqrt_FL)

plot(hist(ctm$orderNorm_FL))
leveneTest(ctm$orderNorm_CTM, group = ctm$Rear_temp) #significant, but also clear from data that variance is lower and maybe less genetic variation
leveneTest(ctm$orderNorm_CTM, group = ctm$DI)

##trying rank normalization
#ctm$rank_ctm <- rank(ctm$CTM)
#shapiro.test(ctm$rank_ctm)
#plot(hist(ctm$rank_ctm))
#leveneTest(ctm$rank_ctm, group = ctm$Rear_temp)

##remove Mix for now since this will be sorted into L,M,H when I have parentage assignments
#ctm<-subset(ctm, DI!="Mix")

##Add in DI calculations
pedigree <- read.table("all_pedigree_meta_DI.txt", header=T)

ctm_pedigree <- merge(pedigree, ctm, by="Fish_ID", all.y = T)


########what are DI range for L, M, and h
ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)

##remove individuals where tank DI doesn't match correct crosses for each DI group. Identified below.
##remove lows with too high of DI: 113V,065C,116G,001S,001H
##remove mediums with too high or too low of DI: "035T","069I","071U"
##remove highs with too low of DI: 117L, 117G,011K, 093B
##remove mixed with wrong crosses: 
ctm_pedigree <- ctm_pedigree[-c(2934,1667,2997,19,8,904,1777,1841,3028,3023,271,2394),]


Low <- subset(ctm_pedigree, DI=="L")
Low <- na.omit(Low)
mean(Low$Offspring_DI) #6.512888
range(Low$Offspring_DI) #5.1423 10.6000


Medium <- subset(ctm_pedigree, DI=="M")
Medium <- na.omit(Medium)
mean(Medium$Offspring_DI) #8.644792
range(Medium$Offspring_DI) #5.85 10.55

High <- subset(ctm_pedigree, DI=="H")
High <- na.omit(High)
mean(High$Offspring_DI) #9.753929
range(High$Offspring_DI) #5.15 10.35

Mixed <- subset(ctm_pedigree, DI=="Mix")
Mixed <- na.omit(Mixed)
mean(Mixed$Offspring_DI) #9.66
range(Mixed$Offspring_DI) #7.90 11.125
##There are lots of O39 individuals from this cross. Male 39 is PC 213 from the high crosses. Male 57 from the mixed crosses is also PC 213, so likely is not an incorrect genetic assignment. Female O is PC 173, but there are no similar PCs from the mixed crosses, but almost all mixed for PC213 are in the "high DI" range. So do I leave those crosses in or remove them?

LowHigh <- subset(ctm_pedigree, assigned_DI=="L" | assigned_DI=="H")
##############

##only remove for RNAseq or epigenetics paper
#ctm<-subset(ctm, DI!="M")

ctm_Mix_1_15<-subset(ctm, Source_tank=="Mix-1-15")
ctm_Mix_2_15<-subset(ctm, Source_tank=="Mix-2-15")
ctm_Mix_1_18<-subset(ctm, Source_tank=="Mix-1-18")
ctm_Mix_2_18<-subset(ctm, Source_tank=="Mix-2-18")

ctm_L_1_15<-subset(ctm, Source_tank=="L-1-15")
ctm_L_2_15<-subset(ctm, Source_tank=="L-2-15")
ctm_L_1_18<-subset(ctm, Source_tank=="L-1-18")
ctm_L_2_18<-subset(ctm, Source_tank=="L-2-18")

ctm_H_1_15<-subset(ctm, Source_tank=="H-1-15")
ctm_H_2_15<-subset(ctm, Source_tank=="H-2-15")
ctm_H_1_18<-subset(ctm, Source_tank=="H-1-18")
ctm_H_2_18<-subset(ctm, Source_tank=="H-2-18")

ctm_M_1_15<-subset(ctm, Source_tank=="M-1-15")
ctm_M_2_15<-subset(ctm, Source_tank=="M-2-15")
ctm_M_1_18<-subset(ctm, Source_tank=="M-1-18")
ctm_M_2_18<-subset(ctm, Source_tank=="M-2-18")

#omit outliers:
grubbs.test(ctm_pedigree$CTM, type = 10)
ctm_pedigree<-subset(ctm_pedigree, Fish_ID!="001H")
ctm_pedigree<-subset(ctm_pedigree, Fish_ID!="001R")

grubbs.test(ctm$CTM, type = 10)
ctm<-subset(ctm, Fish_ID!="001H")
ctm<-subset(ctm, Fish_ID!="001R")


###################################
Random Variable Models
###################################

ctm.model_size = lmer(orderNorm_CTM ~ FL + (1|System:DI), data = ctm)
Anova(ctm.model_size,test.statistic = "F") #significant 2.454e-12! expected since correlated with temp, include as random effect in model 

ctm.model_weight = lmer(orderNorm_CTM ~ Mass + (1|System:DI), data = ctm)
Anova(ctm.model_weight,test.statistic = "F") #significant at 2.2e-16

ctm.model_age = lmer(orderNorm_CTM ~ Age + (1|System:DI), data = ctm)
Anova(ctm.model_age,test.statistic = "F") #significant 0.0004! Also confounded with treatment a little bit though

ctm.model_obs = lmer(orderNorm_CTM ~ Observer + (1|System:DI), data = ctm)
Anova(ctm.model_obs,test.statistic = "F") #significant! 5.354e-12
ctm_obs.emm <- emmeans(ctm.model_obs, ~ Observer, adjust = "tukey")
multcomp::cld(ctm_obs.emm, alpha = 0.05, Letters = LETTERS)

ctm.model_rep = lmer(orderNorm_CTM ~ Rep +(1|System:DI), data = ctm)
Anova(ctm.model_rep,test.statistic = "F") # not significant! 0.6699



ctm.model_size = lmer(orderNorm_CTM ~ FL + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model_size,test.statistic = "F") #significant < 2.2e-16 ***! expected since correlated with temp, include as random effect in model 

ctm.model_weight = lmer(orderNorm_CTM ~ Mass + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model_weight,test.statistic = "F") #significant at < 2.2e-16 ***

ctm.model_age = lmer(orderNorm_CTM ~ Age + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model_age,test.statistic = "F") #significant 0.0003211! Also confounded with treatment a little bit though

ctm.model_obs = lmer(orderNorm_CTM ~ Observer + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model_obs,test.statistic = "F") #significant! 1.03e-11 ***
ctm_obs.emm <- emmeans(ctm.model_obs, ~ Observer, adjust = "tukey")
multcomp::cld(ctm_obs.emm, alpha = 0.05, Letters = LETTERS)

ctm.model_rep = lmer(orderNorm_CTM ~ Rep +(1|System:DI), data = ctm_pedigree)
Anova(ctm.model_rep,test.statistic = "F") # not significant! 0.6699

###################################
#CTM Models Categorical
###################################

ctm.model1 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI), data = ctm) #mixes and 2 outliers removed
Anova(ctm.model1, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI           24.1764  2      6 0.001345 **
Rear_temp    18.2378  1      6 0.005260 **
DI:Rear_temp  2.9171  2      6 0.130329 
'''

ctm.model1.2 = lmer(orderNorm_CTM ~ assigned_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model1.2, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
                           F Df Df.res  Pr(>F)  
assigned_DI           3.4258  2 525.19 0.03325 *
Rear_temp             1.6759  1  13.08 0.21786  
assigned_DI:Rear_temp 1.3112  2 507.84 0.27040
'''
AIC(ctm.model1.2) #8384.512

ctm.model1.L = lmer(orderNorm_CTM ~ Rear_temp + (1|System:DI), data = Low)
Anova(ctm.model1.L, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
Rear_temp 5.1632  1 1.9997  0.151
'''

ctm.model1.M = lmer(orderNorm_CTM ~ Rear_temp + (1|System:DI), data = Medium)
Anova(ctm.model1.M, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
Rear_temp 57.488  1 1.9333 0.01863 *
'''

ctm.model1.H = lmer(orderNorm_CTM ~ Rear_temp + (1|System:DI), data = High)
Anova(ctm.model1.H, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
Rear_temp 0.0844  1 2.0259 0.7984
'''

ctm.model1.Mix = lmer(orderNorm_CTM ~ Rear_temp + (1|System:DI), data = Mixed)
Anova(ctm.model1.Mix, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
Rear_temp 0.8945  1 1.9972 0.4442
'''

ctm.model1.2 = lmer(orderNorm_CTM ~ Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model1.2, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
Rear_temp 1.5544  1 14.139 0.2327
'''

#same as above so sticking with System:DI
ctm.model1.5 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|Source_tank), data = ctm)
Anova(ctm.model1.5, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI           26.1857  2 5.9865 0.001096 **
Rear_temp    19.2617  1 5.9866 0.004649 **
DI:Rear_temp  3.0254  2 5.9866 0.123587
'''


##Models with fork length below
ctm.model2 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI) + (1|FL), data = ctm)
Anova(ctm.model2, test.statistic = "F")
'''
                   F Df Df.res    Pr(>F)    
DI           36.1765  2 6.1897 0.0003846 ***
Rear_temp    56.8859  1 8.8804 3.801e-05 ***
DI:Rear_temp  1.3211  2 6.0851 0.3338124
'''

ctm.model2.1 = lmer(orderNorm_CTM ~ assigned_DI*Rear_temp + (1|System:DI) + (1|FL), data = ctm_pedigree)
Anova(ctm.model2.1, test.statistic = "F")
'''
assigned_DI           10.8079  2 94.172 5.949e-05 ***
Rear_temp             17.5672  1 13.656 0.0009535 ***
assigned_DI:Rear_temp  0.5431  2 92.008 0.5827895 
'''
AIC(ctm.model2.1) #8141.213

ctm.model2.2 = lmer(orderNorm_CTM ~ assigned_DI*Rear_temp + FL + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model2.2, test.statistic = "F")
'''
assigned_DI            7.4941  2  237.85 0.0006979 ***
Rear_temp             12.3404  1   13.62 0.0035789 ** 
FL                    85.1630  1 2864.72 < 2.2e-16 ***
assigned_DI:Rear_temp  1.2050  2  233.79 0.3015324
'''
AIC(ctm.model2.2) #8311.203



##for epigenetics manuscript. two outliers removed
LowHigh <- subset(ctm, DI=="L" | DI=="H")
ctm.model2.4 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI) + (1|FL), data = LowHigh)
Anova(ctm.model2.4, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI           60.3299  1 4.0957 0.001344 **
Rear_temp    26.1668  1 6.8479 0.001471 **
DI:Rear_temp  0.0192  1 4.0302 0.896482 
'''



ctm.model3 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI) + (1|Observer), data = ctm)
Anova(ctm.model3, test.statistic = "F")
'''
DI           40.7683  2 5.9515 0.000336 ***
Rear_temp    26.4629  1 5.9423 0.002190 ** 
DI:Rear_temp  5.2339  2 6.0029 0.048341 *
'''

ctm.model4 = lmer(CTM ~ DI*Rear_temp + (1|System:DI) + (1|FL) + (1|Observer), data = ctm)
Anova(ctm.model4, test.statistic = "F")
'''
DI           34.8562  3  7.6432 8.163e-05 ***
Rear_temp    71.9271  1 12.1939 1.841e-06 ***
DI:Rear_temp  4.2994  3  7.7969   0.04528 *
'''

ctm.model4 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI) + (1|FL) + (1|Observer), data = ctm)
Anova(ctm.model4, test.statistic = "F")
'''
DI           48.4390  2 6.1664 0.0001695 ***
Rear_temp    68.0215  1 9.4643 1.272e-05 ***
DI:Rear_temp  2.2482  2 6.1277 0.1852494    
'''

AIC(ctm.model1, ctm.model2, ctm.model3, ctm.model4) #model 4 has a lower AIC, so a better fit for the data
'''
           df      AIC
ctm.model1  8 6440.649
ctm.model2  9 6304.330
ctm.model3  9 6403.341
ctm.model4 10 6243.829
'''

emmeans(ctm.model1, list(pairwise ~ DI*Rear_temp), adjust = "tukey")
ctm.emm <- emmeans(ctm.model1, ~ DI*Rear_temp, adjust = "tukey")
multcomp::cld(ctm.emm, alpha = 0.05, Letters = LETTERS)

emmeans(ctm.model2, list(pairwise ~ DI*Rear_temp), adjust = "tukey")
ctm.emm <- emmeans(ctm.model2, ~ DI*Rear_temp, adjust = "tukey")
multcomp::cld(ctm.emm, alpha = 0.05, Letters = LETTERS)

emmeans(ctm.model3, list(pairwise ~ DI*Rear_temp), adjust = "tukey")
ctm.emm <- emmeans(ctm.model3, ~ DI*Rear_temp, adjust = "tukey")
multcomp::cld(ctm.emm, alpha = 0.05, Letters = LETTERS)

emmeans(ctm.model4, list(pairwise ~ DI*Rear_temp), adjust = "tukey")
ctm.emm <- emmeans(ctm.model4, ~ DI*Rear_temp, adjust = "tukey")
multcomp::cld(ctm.emm, alpha = 0.05, Letters = LETTERS)





###################################
#CTM Models Continuous
###################################

ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)
ctm.model5 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(ctm.model5, test.statistic = "F")
##run with individuals assigned to wrong crosses above removed, except I didn't remove wrong mixed groups
'''
                            F Df Df.res Pr(>F)
Offspring_DI           0.0573  1 238.440 0.8110
Rear_temp              0.9971  1  13.359 0.3358
Offspring_DI:Rear_temp 2.2472  1 235.393 0.1352
'''

#ctm.model5.5 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|Source_tank), data = ctm_pedigree)
#Anova(ctm.model5.5, test.statistic = "F")
##run with individuals assigned to wrong crosses above removed, except I didn't remove wrong mixed groups
'''
                            F Df Df.res Pr(>F)
Offspring_DI           0.0023  1 240.576 0.9619
Rear_temp              0.6271  1  13.181 0.4425
Offspring_DI:Rear_temp 1.8021  1 226.735 0.1808
'''

ctm.model6 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|System:DI) + (1|AAFam) , data = ctm_pedigree)
Anova(ctm.model6, test.statistic = "F") #takes a while to run, trying with type 3, made significance worse
'''
Offspring_DI           1.9071  1 43.266 0.1744
Rear_temp              1.6129  1 11.855 0.2284
Offspring_DI:Rear_temp 0.8647  1 85.493 0.3551
'''

ctm.model7 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|System:DI) + (1|FL) , data = ctm_pedigree)
Anova(ctm.model7, test.statistic = "F") #takes a while to run, trying with type 3 made significance worse
'''
Offspring_DI           27.1814  1 28.486 1.473e-05 ***
Rear_temp              18.7371  1 14.649 0.0006282 ***
Offspring_DI:Rear_temp  0.0885  1 28.080 0.7682599
'''
AIC(ctm.model7) #7672.4

ctm.model8 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + FL + (1|System:DI) , data = ctm_pedigree)
Anova(ctm.model8, test.statistic = "F") #takes a while to run, trying with type 3 made significance worse
'''
Offspring_DI            9.4243  1   73.33  0.003002 ** 
Rear_temp               9.1661  1   13.70  0.009228 ** 
FL                     78.5466  1 2630.91 < 2.2e-16 ***
Offspring_DI:Rear_temp  0.8195  1   73.66  0.368286  
'''
AIC(ctm.model8) #7827.534

temp15 <- subset(ctm_pedigree, Rear_temp=="15")
temp18 <- subset(ctm_pedigree, Rear_temp=="18")
temp15 <- na.omit(temp15$CTM)
temp18 <- na.omit(temp18$CTM)
mean(temp15) #27.7
mean(temp18) #28.3

ctm.model9 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp*FL + (1|System:DI) , data = ctm_pedigree)
Anova(ctm.model9, test.statistic = "F") #takes a while to run, trying with type 3 made significance worse
'''
Offspring_DI              29.6168  1   28.46 7.877e-06 ***
Rear_temp                 25.1113  1   14.59 0.0001684 ***
FL                        89.3378  1 2007.89 < 2.2e-16 ***
Offspring_DI:Rear_temp     0.0094  1   38.27 0.9234315    
Offspring_DI:FL            0.3538  1 2808.44 0.5520328    
Rear_temp:FL              88.4622  1 1546.68 < 2.2e-16 ***
Offspring_DI:Rear_temp:FL  0.0392  1 2800.26 0.8431426
'''
AIC(ctm.model9) #7773.479

ctm.model10 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp*FL + (1|System:DI) + (1|Trial), data = ctm_pedigree)
Anova(ctm.model10, test.statistic = "F") 
'''
Offspring_DI               23.8933  1   42.24 1.512e-05 ***
Rear_temp                  26.2028  1   17.40 7.968e-05 ***
FL                        119.6595  1 2109.33 < 2.2e-16 ***
Offspring_DI:Rear_temp      0.1308  1   55.44    0.7190    
Offspring_DI:FL             1.5896  1 2781.59    0.2075    
Rear_temp:FL              131.5707  1 1682.78 < 2.2e-16 ***
Offspring_DI:Rear_temp:FL   0.4925  1 2789.54    0.4829
'''

library(effectsize)
eta_squared(ctm.model9, partial = TRUE) #effectsize package

##################
install.packages("partR2")
install.packages("r2glmm")
library(partR2)
library(r2glmm)

# Using partR2, doesn't work
 result_partr2 <- partR2(ctm.model9, partially_dependent_on = Rear_temp, R2_type = "marginal")

# Using r2glmm, doesn't work
r2_partial <- r2beta(ctm.model9, method = "partial", data = ctm_pedigree)
###############


ctm.model10 = lmer(orderNorm_FL ~ Offspring_DI*Rear_temp + (1|System:DI) , data = ctm_pedigree)
Anova(ctm.model10, test.statistic = "F") 
'''
Offspring_DI           27.4046  1 1114.99 1.971e-07 ***
Rear_temp               7.2130  1   14.03   0.01772 *  
Offspring_DI:Rear_temp  0.0423  1 1067.37   0.83711
'''


ctm_pedigree_largefish <- subset(ctm_pedigree, FL > 15)
ctm.model11 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree_largefish)
Anova(ctm.model11, test.statistic = "F")
'''
Offspring_DI           6.2333  1 61.060 0.01525 *
Rear_temp              6.7398  1 12.585 0.02265 *
Offspring_DI:Rear_temp 0.7075  1 61.048 0.40356 
'''

ctm_pedigree_largefish <- subset(ctm_pedigree, FL > 15)
ctm.model12 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + FL + (1|System:DI), data = ctm_pedigree_largefish)
Anova(ctm.model12, test.statistic = "F")
'''
Offspring_DI           7.3482  1   55.36 0.008922 **
Rear_temp              7.9396  1   13.66 0.013965 * 
FL                     0.8543  1 2406.69 0.355435   
Offspring_DI:Rear_temp 0.6702  1   56.33 0.416442  
'''

###################################
#FL Models Categorical
###################################

FL.model1 = lmer(orderNorm_FL ~ DI*Rear_temp + (1|System:DI), data = ctm)
Anova(FL.model1, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI            1.3281  3 8.1474 0.3304131    
Rear_temp    37.3181  1 8.1462 0.0002665 ***
DI:Rear_temp  1.0883  3 8.1439 0.4071232
'''
AIC(FL.model1) #7635.195

FL.model1.2 = lmer(orderNorm_FL ~ assigned_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(FL.model1.2, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI            1.3281  3 8.1474 0.3304131    
Rear_temp    37.3181  1 8.1462 0.0002665 ***
DI:Rear_temp  1.0883  3 8.1439 0.4071232
'''



###################################
#FL Models Continuos
###################################

Anova(FL.model2, test.statistic = "F")

'''
F Df Df.res   Pr(>F)   
Offspring_DI           49.2556  1 1257.7 3.664e-12 ***
Rear_temp              17.4829  1   14.1 0.0009108 ***
Offspring_DI:Rear_temp  1.4322  1 1216.3 0.2316339  
'''
#with 2 outliers removed:
'''
Offspring_DI           54.1827  1 1198.27 3.384e-13 ***
Rear_temp              15.1568  1   13.67  0.001692 ** 
Offspring_DI:Rear_temp  0.4942  1 1186.08  0.482220
'''
summary(FL.model2)
'''
                         Estimate Std. Error t value
(Intercept)              -1.82665    0.31728  -5.757
Offspring_DI              0.14793    0.02957   5.003
Rear_temp18               0.69478    0.46428   1.496
Offspring_DI:Rear_temp18  0.03108    0.04371   0.711
'''

r2(FL.model2) #from performance package
'''
Conditional R2: 0.468
Marginal R2: 0.267
'''



FL.model3 = lmer(center_FL ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(FL.model3, test.statistic = "F")
'''
Offspring_DI           45.1053  1 1113.82 2.971e-11 ***
Rear_temp              17.5895  1   13.65 0.0009498 ***
Offspring_DI:Rear_temp  2.0437  1 1064.23 0.1531295
'''

FL.model4 = lmer(FL ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(FL.model4, test.statistic = "F")
'''
Offspring_DI           45.1053  1 1113.82 2.971e-11 ***
Rear_temp              17.5894  1   13.65 0.0009498 ***
Offspring_DI:Rear_temp  2.0437  1 1064.23 0.1531295 
'''



AIC(FL.model2) #6929.319

###################################
#Random effect plots
###################################

###plotting CTM by age
windows()
ggplot(ctm, aes(x=Age, y=CTM)) + geom_point() + geom_smooth(method=lm)

###plotting CTM by FL
windows()
ggplot(ctm, aes(x=FL, y=CTM)) + geom_point() + geom_smooth(method=lm)
ggplot(ctm, aes(x=FL, y=orderNorm_CTM)) + geom_point() + geom_smooth(method=lm)


###################################
#Main effect plots Categorical
###################################
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      se = (sd(x[[col]]))/sqrt(length(x[[col]])))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df1 <- data_summary(ctm, varname="CTM", 
                    groupnames=c("DI", "Rear_temp"))
head(df1)

df1.1 <- data_summary(ctm_pedigree, varname="CTM", 
                    groupnames=c("Offspring_DI", "Rear_temp"))
head(df1.1)

df1.2 <- data_summary(ctm_pedigree, varname="CTM", 
                      groupnames=c("Rear_temp"))
head(df1.2)

df2 <- data_summary(ctm, varname="FL", 
                    groupnames=c("DI", "Rear_temp"))
head(df2)

df_LH <- data_summary(LowHigh, varname="CTM", 
                    groupnames=c("DI", "Rear_temp"))
head(df_LH)

df_LH <- data_summary(LowHigh, varname="CTM", 
                      groupnames=c("Rear_temp"))
head(df_LH)


df1$Rear_temp <- factor(df1$Rear_temp,levels = c("15", "18"))
df1$DI <- factor(df1$DI,levels = c("W", "L", "M", "H"))

windows()
ggplot(data=ctm, aes(x=DI, y=CTM, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  geom_smooth(method=lm, aes(fill=Rear_temp))+
  labs(y=expression("CTM (℃)"), x="Domestication Index", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
 # scale_y_continuous(breaks=seq(20,30,1))


ctm$DI <- factor(ctm$DI,levels = c("L", "M", "H", "Mix"))
windows()
ggplot(data=ctm, aes(x=DI, y=CTM, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))

######################
#Epigenetics Manuscript
######################
LowHigh$DI <- factor(LowHigh$DI,levels = c("L", "H"))
LowHigh$Rear_temp <- factor(LowHigh$Rear_temp,levels = c("15", "18"))
windows()
ggplot(data=LowHigh, aes(x=DI, y=CTM, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))
################


ctm_pedigree$assigned_DI <- factor(ctm_pedigree$assigned_DI,levels = c("L", "M", "H"))
ctm_pedigree_plot <- ctm_pedigree[!is.na(ctm_pedigree$assigned_DI),] # remove NAs for DI
windows()
ggplot(data=ctm_pedigree_plot, aes(x=assigned_DI, y=CTM, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))


#CTM normalized by FL

windows()
ggplot(data=ctm_pedigree, aes(x=Offspring_DI, y=FL_CTM, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  geom_smooth(method=lm, aes(fill=Rear_temp))+  labs(y=expression("CTM / FL"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
  #scale_y_continuous(breaks=seq(15,34,1))

#CTM normalized by FL CONTINUOUS
ctm_pedigree$assigned_DI <- factor(ctm_pedigree$assigned_DI,levels = c("L", "M", "H"))
ctm_pedigree_plot <- ctm_pedigree[!is.na(ctm_pedigree$assigned_DI),]
ctm_pedigree_plot$CTMFL <- ctm_pedigree_plot$CTM / ctm_pedigree_plot$FL# remove NAs for DI
windows()
ggplot(data=ctm_pedigree_plot, aes(x=assigned_DI, y=CTMFL, color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM / FL"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
#scale_y_continuous(breaks=seq(15,34,1))


###################################
#Main effect plots Continuous
###################################

##create equation
windows()
# m < lm(put the important variables here)
# X2hr_acclim X CTM
# this part of the code creates the equation of the line and the R2 value
lm_eqn <- function(df){
  m <- lm(CTM ~ Offspring_DI, ctm_pedigree);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)
windows()
plot1 <- ggplot(data=ctm_pedigree, aes(x=Offspring_DI, y=CTM, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  geom_smooth(method=lm, aes(colour = Rear_temp, group = Rear_temp)) +
  #geom_text(x = 17, y = 33, label = lm_eqn(ctm_pedigree), parse = TRUE) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="DI", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), 
        legend.text = element_text(size=14), legend.position = "none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))


##CTM by FL and rearing temp
ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)
windows()
plot2 <- ggplot(data=ctm_pedigree, aes(x=FL, y=CTM, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  geom_smooth(method=lm, aes(colour = Rear_temp, group = Rear_temp)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="FL (mm)", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), 
        legend.text = element_text(size=14), legend.position = "none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))

##FL by DI and rearing temp
ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)
windows()
plot3 <- ggplot(data=ctm_pedigree, aes(x=Offspring_DI, y=FL, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  geom_smooth(method=lm, aes(colour = Rear_temp, group = Rear_temp)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("FL (mm)"), x="DI", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(10,34,1))



windows()
plot_grid(plot1, plot2, plot3, labels = "AUTO")

#####################
boxplot for Observer Bias
#####################

windows()
ggplot(data=ctm, aes(x=Observer, y=CTM,  color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  labs(y=expression("CTM (?C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))


windows()
ggplot(data=ctm_pedigree, aes(x=Observer, y=CTM,  color=Rear_temp)) +
  geom_boxplot() +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  labs(y=expression("CTM (?C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
  