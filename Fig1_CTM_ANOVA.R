
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

#### FINAL MODEL
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




###################################
#FL Models Continuos
###################################


FL.model3 = lmer(center_FL ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(FL.model3, test.statistic = "F")
'''
Offspring_DI           45.1053  1 1113.82 2.971e-11 ***
Rear_temp              17.5895  1   13.65 0.0009498 ***
Offspring_DI:Rear_temp  2.0437  1 1064.23 0.1531295
'''

##### FINAL MODEL
FL.model4 = lmer(FL ~ Offspring_DI*Rear_temp + (1|System:DI), data = ctm_pedigree)
Anova(FL.model4, test.statistic = "F")
'''
Offspring_DI           45.1053  1 1113.82 2.971e-11 ***
Rear_temp              17.5894  1   13.65 0.0009498 ***
Offspring_DI:Rear_temp  2.0437  1 1064.23 0.1531295 
'''


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
## boxplot for Observer Bias
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

  
