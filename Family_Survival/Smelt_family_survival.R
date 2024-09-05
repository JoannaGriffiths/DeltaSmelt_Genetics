
library("dplyr")
library("ggplot2")
library(tidyr)

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign")

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

#family_count$Family <- factor(family_count$Family,levels = c(
write.table(family_count_eggs, file="Family_survival.txt",sep = "\t", row.names = F)

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
  
##NEXT STEPS
  #Check egg count to make sure its a percent survival calculted if they had low egg count numbers to begin with
  
  
#################### Models
##Redo with percent survival next
library(lme4)
library(car)
library(emmeans)
library(ggplot2)
  
ctm <- ctm %>% separate_wider_delim(Treatment, "-", names = c("DI", "Replicate", "Temp"))

##choose one or the other below
treatment_family_counts <- ctm %>%
  group_by(Temp, Replicate) %>%
  count(AAFam)

treatment_family_counts <- ctm %>%
  group_by(Temp) %>%
  count(AAFam)

treatment_family_counts$Temp_by_Fam <- paste(treatment_family_counts$Temp, treatment_family_counts$AAFam, sep="_")

##manipulate egg counts file so you can merge it with treatment_family_counts file
eggs$mean_15 <- ((eggs$R1.15 + eggs$R2.15)/2)
eggs$mean_18 <- ((eggs$R1.18 + eggs$R2.18)/2)
eggs$total_15 <- (eggs$R1.15 + eggs$R2.15)
eggs$total_18 <- (eggs$R1.18 + eggs$R2.18)

eggs_15 <- eggs[,c(1,7,9)]
colnames(eggs_15) <- c("Family", "Mean", "Total")

eggs_18 <- eggs[,c(1,8,10)]
colnames(eggs_18) <- c("Family", "Mean", "Total")

eggs_15$Temp_by_Fam <- paste("15_", eggs_15$Family, sep="")
eggs_18$Temp_by_Fam <- paste("18_", eggs_18$Family, sep="")

eggs_mean <- rbind(eggs_15, eggs_18)

treatment_family_percent <- merge(treatment_family_counts, eggs_mean, by="Temp_by_Fam")
treatment_family_percent$survival <- treatment_family_percent$n / treatment_family_percent$Total



##Graph family survival
windows()
ggplot(data=treatment_family_percent, aes(x=Family, y=survival, fill=Temp)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("lightsteelblue3", "lightsalmon")) +
  labs(y=expression("Recovered offspring per family (%)"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=10), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14), axis.text.x = element_text(angle = 90))





##model not correct 
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
