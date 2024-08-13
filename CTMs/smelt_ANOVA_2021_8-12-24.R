
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
library("multcomp")
library("bestNormalize")
library("outliers")

setwd("~/UCDavis/FCCL/Spawning_2021/CTMs")

ctm <- read.table("2021_CTM_data_trial1-120.txt", header = T)
ctm$CTM <- as.numeric(ctm$CTM)
ctm$FL <- as.numeric(ctm$FL)
ctm$DI<- as.factor(ctm$DI)
ctm$Rear_temp<- as.factor(ctm$Rear_temp)
shapiro.test(ctm$CTM)
plot(hist(ctm$CTM))
leveneTest(ctm$CTM, group = ctm$Rear_temp)
leveneTest(ctm$CTM, group = ctm$DI)


#bestNormalize(ctm$CTM) #suggests orderNorm is best transformation
ctm$orderNorm_CTM <- predict(orderNorm(ctm$CTM))
shapiro.test(ctm$orderNorm_CTM)
plot(hist(ctm$orderNorm_CTM))
leveneTest(ctm$orderNorm_CTM, group = ctm$Rear_temp) #significant, but also clear from data that variance is lower and maybe less genetic variation
leveneTest(ctm$orderNorm_CTM, group = ctm$DI)

##trying rank normalization
ctm$rank_ctm <- rank(ctm$CTM)
shapiro.test(ctm$rank_ctm)
plot(hist(ctm$rank_ctm))
leveneTest(ctm$rank_ctm, group = ctm$Rear_temp)

##remove Mix for now since this will be sorted into L,M,H when I have parentage assignments
#ctm<-subset(ctm, DI!="Mix")

##Add in DI calculations
pedigree <- read.table("../../../Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign/all_pedigree_meta_DI.txt", header=T)
colnames(pedigree) <- c("AAFam", "Dam_Seq_ID", "Sire_Seq_ID", "Seq_ID", "Fish_ID", "Treatment", "Sire.1", "Sire.2", "Dam.1", "Dam.2", "Offspring_DI")

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


###################################
CTM Models Categorical
###################################

ctm.model1 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI), data = ctm)
Anova(ctm.model1, test.statistic = "F")
'''
                   F Df Df.res   Pr(>F)   
DI           26.7755  2 5.9967 0.001025 **
Rear_temp    20.4109  1 5.9967 0.004032 **
DI:Rear_temp  3.0415  2 5.9968 0.122483
'''

ctm.model2 = lmer(orderNorm_CTM ~ DI*Rear_temp + (1|System:DI) + (1|FL), data = ctm)
Anova(ctm.model2, test.statistic = "F")
'''
                   F Df Df.res    Pr(>F)    
DI           36.1765  2 6.1897 0.0003846 ***
Rear_temp    56.8859  1 8.8804 3.801e-05 ***
DI:Rear_temp  1.3211  2 6.0851 0.3338124
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
Offspring_DI           0.3518  1 283.367 0.5536
Rear_temp              0.7888  1  13.456 0.3901
Offspring_DI:Rear_temp 2.7032  1 278.501 0.1013
'''

ctm.model6 = lmer(orderNorm_CTM ~ Offspring_DI*Rear_temp + (1|System:DI) + (1|AAFam) , data = ctm_pedigree)
Anova(ctm.model6, test.statistic = "F") #takes a while to run, trying with type 3, made significance worse
'''
Offspring_DI           1.9071  1 43.266 0.1744
Rear_temp              1.6129  1 11.855 0.2284
Offspring_DI:Rear_temp 0.8647  1 85.493 0.3551
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


ctm$DI <- factor(ctm$DI,levels = c("W", "L", "M", "H"))
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



###################################
#Main effect plots Continuous
###################################

ctm_pedigree$Offspring_DI <- as.numeric(ctm_pedigree$Offspring_DI)
windows()
ggplot(data=ctm_pedigree, aes(x=Offspring_DI, y=CTM, color=Rear_temp)) +
  geom_point(aes(fill = Rear_temp), size = 1, shape = 21) +
  geom_smooth(method=lm, aes(colour = Rear_temp, group = Rear_temp)) +
  scale_color_manual(values=c("lightsteelblue4", "lightsalmon4")) +
  scale_fill_manual(values=c("lightsteelblue2", "lightsalmon")) +
  labs(y=expression("CTM (◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))+
  scale_y_continuous(breaks=seq(15,34,1))




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
  