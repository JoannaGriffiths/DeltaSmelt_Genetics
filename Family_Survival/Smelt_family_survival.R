
library("dplyr")
library("ggplot2")
library(tidyr)
library("bestNormalize")

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign")

load("smelt_survival1.RData")
load("smelt_survival2.RData")

ctm <- read.delim2("all_pedigree_meta_DI.txt", header = T)
ctm <- na.omit(ctm)

eggs <- read.delim2("egg_counts.txt", header=T)

family_count <- ctm %>% count(AAFam) #started with 121 families, left with 104 at end
family_count$AAFam <- as.character(family_count$AAFam)
colnames(family_count) <- c("Family", "Count")
family_count_eggs <- merge(family_count, eggs, by="Family", all.x = T)
family_count_eggs$survival <- family_count_eggs$Count/family_count_eggs$Total

Dam_count <- ctm %>% count(Dam.2) #started with 121 families, left with 90 at end
Dam_count$Dam.2 <- as.character(Dam_count$Dam.2)
colnames(Dam_count) <- c("Dam", "Count")

Sire_count <- ctm %>% count(Sire.2) #started with 121 families, left with 90 at end
Sire_count$Sire.2 <- as.character(Sire_count$Sire.2)
colnames(Sire_count) <- c("Sire", "Count")


#write.table(family_count_eggs, file="Family_survival.txt",sep = "\t", row.names = F)

##Graph family counts
windows()
ggplot(data=family_count, aes(x=Family, y=Count)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    #geom_errorbar(data=dat, aes(ymin=CTM-ster, ymax=CTM+ster), width=.2, position=position_dodge(.7)) +
    #scale_fill_manual(values=c("#3399FF", "#FF9900")) +
    labs(y=expression("Offspring per family"), x="", fill="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))


##Graph family survival
windows()
ggplot(data=family_count_eggs, aes(x=Family, y=survival)) +
  geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
  #geom_errorbar(data=dat, aes(ymin=CTM-ster, ymax=CTM+ster), width=.2, position=position_dodge(.7)) +
  #scale_fill_manual(values=c("#3399FF", "#FF9900")) +
  labs(y=expression("Recovered offspring per family (%)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))

##Graph Dam counts
windows()
ggplot(data=Dam_count, aes(x=Dam, y=Count)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    #geom_errorbar(data=dat, aes(ymin=CTM-ster, ymax=CTM+ster), width=.2, position=position_dodge(.7)) +
    #scale_fill_manual(values=c("#3399FF", "#FF9900")) +
    labs(y=expression("Offspring per Dam"), x="", fill="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))

  
##Graph Sire counts
Sire_count$Sire <- factor(Sire_count$Sire, levels=c("1","2","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59","60","61","62","63","64","65"))
  windows()
  ggplot(data=Sire_count, aes(x=Sire, y=Count)) +
    geom_bar(stat="identity", position=position_dodge(), colour="black", width = 0.7) +
    #geom_errorbar(data=dat, aes(ymin=CTM-ster, ymax=CTM+ster), width=.2, position=position_dodge(.7)) +
    #scale_fill_manual(values=c("#3399FF", "#FF9900")) +
    labs(y=expression("Offspring per Sire"), x="", fill="") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
          axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
          strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))
  

  
#################### Models

library(lme4)
library(car)
library(emmeans)
library(ggplot2)
  
ctm <- ctm %>% separate_wider_delim(Treatment, "-", names = c("DI", "Replicate", "Temp"))

##choose one or the other below
treatment_family_counts <- ctm %>%
  group_by(Temp, Replicate) %>%
  count(AAFam)

#treatment_family_counts <- ctm %>%
#  group_by(Temp) %>%
#  count(AAFam)

#ctm$System <- paste(ctm$Replicate, ctm$Temp, sep = "-")
#treatment_family_counts <- ctm %>%
#  group_by(Temp, System) %>%
#  count(AAFam)


##only do below if ran first two options above
treatment_family_counts$Temp_by_Fam <- paste(treatment_family_counts$Temp, treatment_family_counts$AAFam, sep="_")

#get means for each family by temperature replicate
treatment_family_counts$n <- as.numeric(treatment_family_counts$n)
treatment_family_means <- treatment_family_counts %>%
  group_by(Temp_by_Fam) %>%
  summarise(mean_count = mean(n, na.rm=TRUE))


##manipulate egg counts file so you can merge it with treatment_family_counts file
eggs$mean_15 <- ((eggs$R1.15 + eggs$R2.15)/2)
eggs$mean_18 <- ((eggs$R1.18 + eggs$R2.18)/2)
eggs$total_15 <- (eggs$R1.15 + eggs$R2.15)
eggs$total_18 <- (eggs$R1.18 + eggs$R2.18)

eggs_15 <- eggs[,c(1,7,9)]
colnames(eggs_15) <- c("Family", "Mean_eggs", "Total")

eggs_18 <- eggs[,c(1,8,10)]
colnames(eggs_18) <- c("Family", "Mean_eggs", "Total")

eggs_15$Temp_by_Fam <- paste("15_", eggs_15$Family, sep="")
eggs_18$Temp_by_Fam <- paste("18_", eggs_18$Family, sep="")

eggs_mean <- rbind(eggs_15, eggs_18)

treatment_family_percent <- merge(treatment_family_means, eggs_mean, by="Temp_by_Fam")
treatment_family_percent$Mean_survival <- treatment_family_percent$mean_count / treatment_family_percent$Mean_eggs
colnames(treatment_family_percent) <- c("Temp_by_Fam", "mean_count", "AAFam", "Mean_eggs", "Total", "Mean_survival")


treatment_family_percent2 <- treatment_family_percent %>% 
  mutate(
    Temp_by_Fam = Temp_by_Fam
    ) %>%
  separate(Temp_by_Fam, into=c("Temp", "Fam"), sep="_")

family_DI <- read.delim2("Family_DIs.txt", header=T)

treatment_family_percent_DI <- merge(treatment_family_percent2, family_DI, by="AAFam")

save(ctm, Dam_count, eggs, eggs_15, eggs_18, eggs_mean, family_count, family_count_eggs, family_DI, Sire_count, treatment_family_counts, treatment_family_means, treatment_family_percent, treatment_family_percent_DI, treatment_family_percent2,
     file="smelt_survival1.RData")



## Low vs High Comparison with fam by replicate treatment ###

eggs_15_R1 <- eggs[,c(1,3,9)]
colnames(eggs_15_R1) <- c("Family", "Mean_eggs", "Total")

eggs_15_R2 <- eggs[,c(1,4,9)]
colnames(eggs_15_R2) <- c("Family", "Mean_eggs", "Total")

eggs_18_R1 <- eggs[,c(1,5,10)]
colnames(eggs_18_R1) <- c("Family", "Mean_eggs", "Total")

eggs_18_R2 <- eggs[,c(1,6,10)]
colnames(eggs_18_R2) <- c("Family", "Mean_eggs", "Total")

eggs_15_R1$treatment_Fam <- paste(eggs_15_R1$Family, "_1_15", sep="")
eggs_15_R2$treatment_Fam <- paste(eggs_15_R2$Family, "_2_15", sep="")
eggs_18_R1$treatment_Fam <- paste(eggs_18_R1$Family, "_1_18", sep="")
eggs_18_R2$treatment_Fam <- paste(eggs_18_R2$Family, "_2_18", sep="")

eggs_mean_R <- rbind(eggs_15_R1, eggs_15_R2, eggs_18_R1, eggs_18_R2)

treatment_family_percent_R <- merge(treatment_family_counts, eggs_mean_R, by="treatment_Fam")
treatment_family_percent_R$survival <- treatment_family_percent_R$n / treatment_family_percent_R$Mean_eggs
treatment_family_percent_R$Temp_Rep <- paste(treatment_family_percent_R$Temp, treatment_family_percent_R$Replicate, sep="_")

slopes_Rep <- treatment_family_percent_R %>%
  group_by(AAFam) %>%
  filter(n_distinct(Temp_Rep) == 4)

slopes_Rep_filtered <- subset(slopes_Rep, survival != "Inf")

slopes <- slopes_Rep_filtered %>%
  group_by(AAFam) %>%
  summarize(slope = coef(lm(survival ~ Temp))[2]) #2 keeps only slopes

slope_DI <- merge(slopes, family_DI, by="AAFam")
slope_DI$Offspring_DI <- as.numeric(slope_DI$Offspring_DI)

slope_DI_L <- subset(slope_DI, Offspring_DI < 7)
slope_DI_H <- subset(slope_DI, Offspring_DI > 9)
slope_DI_L$DI <- 'L'
slope_DI_H$DI <- 'H'
slope_DI_LH <- rbind(slope_DI_L, slope_DI_H)
slope_DI_LH$DI <- factor(slope_DI_LH$DI,levels = c("L", "H"))


#### Survival Slope Models

slopes.model1 = lm(slope ~ DI, data = slope_DI_LH)
Anova(slopes.model1, test.statistic = "F")
'''
Response: slope
           Sum Sq Df F value Pr(>F)
DI        0.00090  1  0.1009 0.7521
Residuals 0.42996 48  
'''

windows()
ggplot(data=slope_DI_LH, aes(x=DI, y=slope, color=DI)) +
  geom_boxplot() +
  geom_point(aes(fill = DI), size = 1, shape = 21, position = position_jitterdodge()) +
  scale_color_manual(values=c("grey", "black")) +
  scale_fill_manual(values=c("grey", "black")) +
  labs(y=expression("plasticity (family slope 15◦C-18◦C)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.position="none",
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))


#######################
##Graph family survival
#########################

windows()
ggplot(data=treatment_family_percent_DI, aes(x=AAFam, y=Mean_survival, fill=Temp)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  labs(y=expression("Recovered offspring per family (%)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))


####Graph family survival by rear temp and DI
windows()
ggplot(data=treatment_family_percent_DI, aes(x=Offspring_DI, y=Mean_survival, color=Temp)) +
  geom_point(aes(fill = Temp), size = 1, shape = 21) +
  scale_color_manual(values=c("lightsteelblue3", "lightsalmon")) +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  geom_smooth(method=lm, aes(fill=Temp))+
  labs(y=expression("Percent Recovered Offspring Per Family (%)"), x="Domestication Index", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
# scale_y_continuous(breaks=seq(20,30,1))


##Graph family plasticity
treatment_family_percent_DI$Temp <- factor(treatment_family_percent_DI$Temp,levels = c("15", "18"))
windows()
ggplot(data=treatment_family_percent_DI, aes(x=Temp, y=Mean_survival, group = AAFam)) +
  geom_point(size = 1, shape = 21) +
  geom_line() +
  labs(y=expression("Percent Recovered Offspring Per Family (%)"), x="Rearing Temperature", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
# scale_y_continuous(breaks=seq(20,30,1))


########################
##models
########################

#survival.model1 = lm(n ~ Temp*AAFam, data = treatment_family_counts)  
survival.model1 = lm(n ~ Temp, data = treatment_family_counts)

survival.model2 = lmer(percent_survival ~ Offspring_DI*Rear_temp + (1|System:DI), data = family_survival)

Anova(survival.model1, test.statistic = "F")
'''
Response: n
Sum Sq  Df F value Pr(>F)
Temp        12.9   1  0.5186 0.4719
Residuals 9479.3 381
'''

egg_counts_model <- read.delim2("egg_counts_for_model.txt", header=T)
egg_counts_model$treatment_Fam <- paste(egg_counts_model$Family, egg_counts_model$Replicate, egg_counts_model$Temp, sep="_")
treatment_family_counts$treatment_Fam <- paste(treatment_family_counts$AAFam, treatment_family_counts$Replicate, treatment_family_counts$Temp, sep="_")

treatment_family_model <- merge(treatment_family_counts, egg_counts_model, by="treatment_Fam")
treatment_family_model$survival <- treatment_family_model$n / treatment_family_model$Egg_counts

treatment_family_model_DI <- merge(treatment_family_model, family_DI, by="AAFam")

treatment_family_model_DI$survival <- gsub("Inf","0",as.character(treatment_family_model_DI$survival))
treatment_family_model_DI$survival <- as.numeric(treatment_family_model_DI$survival)
treatment_family_model_DI$Offspring_DI <- as.numeric(treatment_family_model_DI$Offspring_DI)

treatment_family_model_DI$System <- paste(treatment_family_model_DI$Temp.x, treatment_family_model_DI$Replicate.x, sep="_")
family_ID <- read.delim2("Family_IDs.txt", header=T)
family_ID <- family_ID[,c(1,8,9)]
colnames(family_ID) <- c("AAFam", "DI", "assigned_DI")

treatment_family_model_DI2 <- merge(treatment_family_model_DI, family_ID, by="AAFam")

save(egg_counts_model, treatment_family_counts, treatment_family_model, treatment_family_model_DI, treatment_family_model_DI2,
     file="smelt_survival2.RData")

shapiro.test(treatment_family_model_DI2$survival) #p-value = 9.326e-11
plot(hist(treatment_family_model_DI2$survival))
#bestNormalize(treatment_family_model_DI2$survival) #best was Yeo-Johnson transformation (Box-Cox transformation)
treatment_family_model_DI2$survival_norm <- predict(yeojohnson(treatment_family_model_DI2$survival, standardize = T))
shapiro.test(treatment_family_model_DI2$survival_norm) #p-value = 0.01578
plot(hist(treatment_family_model_DI2$survival_norm))

survival.model3 = lmer(survival_norm ~ Offspring_DI*Temp.x + (1|System:DI), data = treatment_family_model_DI2)
Anova(survival.model3, test.statistic = "F")
'''
                         F Df Df.res Pr(>F)
Offspring_DI        0.0846  1 19.139 0.7743
Temp.x              0.6695  1 11.666 0.4296
Offspring_DI:Temp.x 1.0442  1 19.146 0.3196
'''

treatment_family_model_DI2_LH$survival_norm <- predict(yeojohnson(treatment_family_model_DI2_LH$survival, standardize = T))
survival.model3.LH = lmer(survival_norm ~ assigned_DI*Temp.x + (1|System:DI), data = treatment_family_model_DI2_LH)
Anova(survival.model3.LH, test.statistic = "F")
'''
                         F Df Df.res Pr(>F)
assigned_DI        4.1388  1 4.0996 0.1099
Temp.x             0.0160  1 3.9116 0.9056
assigned_DI:Temp.x 1.9556  1 4.0995 0.2329
'''

##Use this model to look at effects of DI and Temp
survival.model4 = lmer(survival_norm ~ Offspring_DI*Temp.x + (1|System:DI) + (1|AAFam), data = treatment_family_model_DI2)
Anova(survival.model4, test.statistic = "F")
'''
                         F Df Df.res Pr(>F)
Offspring_DI        0.1855  1 57.088 0.6683
Temp.x              0.9212  1 10.320 0.3591
Offspring_DI:Temp.x 1.6042  1 22.054 0.2185
'''

##Use this model to look at effects of Fam and Temp
survival.model5 = lmer(survival_norm ~ AAFam*Temp.x + (1|System:DI), data = treatment_family_model_DI2)
Anova(survival.model5, test.statistic = "F")
'''
fixed-effect model matrix is rank deficient so dropping 8 columns / coefficients
                  F  Df  Df.res  Pr(>F)    
AAFam        5.7787 102 144.504 < 2e-16 ***
Temp.x       1.1038   1   7.904 0.32449    
AAFam:Temp.x 1.5999  94 147.319 0.00526 **
'''

