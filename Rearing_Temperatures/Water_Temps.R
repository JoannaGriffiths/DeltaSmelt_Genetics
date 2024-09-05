
setwd("~/UCDavis/FCCL/Spawning_2021/Temperatures")
library("ggplot2")

egg_temp <- read.delim("egg_temps.txt", header=T)
tray_temp <- read.delim("tray_temps.txt", header=T)

all_temps <- read.delim("life_history_temps.txt", header=T)


'''
windows()
ggplot(data=egg_temp, aes(x=Temp, fill=System)) +
  geom_histogram(binwidth=0.25, position=position_dodge(0.2)) +
  scale_fill_manual(values=c("gray", "black", "gray", "black")) +
  labs(y=expression("# of Occurrences"), x="", fill="") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=14), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))
'''
all_temps <- subset(all_temps, Location=="Tanks")
windows()
ggplot(data=all_temps, aes(x=factor(Day_Time), y=Temp, group=System, color=System)) +
  geom_line() +
  scale_color_manual(values=c("black", "red", "blue", "orange")) +
  labs(y=expression("temp (°C)"), x="", fill="") +
  theme(axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=9), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))

tray_temp$Date_Time <- paste(tray_temp$Date, tray_temp$Temp, sep = ".")

windows()
ggplot(data=tray_temp, aes(x=factor(Date_Time), y=Temp, group=System, color=System)) +
  geom_line() +
  scale_color_manual(values=c("red", "blue")) +
  labs(y=expression("temp (°C)"), x="", fill="") +
  theme(axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=9), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))



################################ tanks
tank_temp <- read.delim("Tank_temp_for_R.txt", header=T)
tank_temp$Rear_Temp <- as.factor(tank_temp$Rear_temp)
tank_temp$Replicate <- as.factor(tank_temp$Replicate)
tank_temp$System <- as.factor(tank_temp$System)
tank_temp$Temp <- as.numeric(tank_temp$Temp)

Data_1_15 <- subset(tank_temp, tank_temp$System=="1-15")
Data_2_15 <- subset(tank_temp, tank_temp$System=="2-15")
Data_1_18 <- subset(tank_temp, tank_temp$System=="1-18")
Data_2_18 <- subset(tank_temp, tank_temp$System=="2-18")
mean(Data_1_15$Temp, na.rm = TRUE) #17.62
mean(Data_2_15$Temp, na.rm = TRUE) #17.29
mean(Data_1_18$Temp, na.rm = TRUE) #18.72
mean(Data_2_18$Temp, na.rm = TRUE) #18.44

Temp_1_15 <- aggregate(Temp~Date, data=Data_1_15, FUN=function(x) c(mean=mean(x)))
Temp_1_15$System <- rep("1-15",nrow(Temp_1_15))

Temp_2_15 <- aggregate(Temp~Date, data=Data_2_15, FUN=function(x) c(mean=mean(x)))
Temp_2_15$System <- rep("2-15",nrow(Temp_2_15))

Temp_1_18 <- aggregate(Temp~Date, data=Data_1_18, FUN=function(x) c(mean=mean(x)))
Temp_1_18$System <- rep("1-18",nrow(Temp_1_18))

Temp_2_18 <- aggregate(Temp~Date, data=Data_2_18, FUN=function(x) c(mean=mean(x)))
Temp_2_18$System <- rep("2-18",nrow(Temp_2_18))

Temp_data <- rbind(Temp_1_15, Temp_2_15, Temp_1_18, Temp_2_18)

windows()
ggplot(data=Temp_data, aes(x=factor(Date), y=Temp, group=System, color=System)) +
  geom_line() +
  scale_color_manual(values=c("blue", "red", "blue", "red")) +
  labs(y=expression("temp (°C)"), x="", fill="") +
  theme(axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=6), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))



######all tank data pre and post CTMs

tank_temp <- read.delim("Tank_temp_for_R.txt", header=T)
tank_temp$Rear_Temp <- as.factor(tank_temp$Rear_temp)
tank_temp$Replicate <- as.factor(tank_temp$Replicate)
tank_temp$System <- as.factor(tank_temp$System)
tank_temp$Temp <- as.numeric(tank_temp$Temp)

Data_1_15 <- subset(tank_temp, tank_temp$System=="1-15")
Data_2_15 <- subset(tank_temp, tank_temp$System=="2-15")
Data_1_18 <- subset(tank_temp, tank_temp$System=="1-18")
Data_2_18 <- subset(tank_temp, tank_temp$System=="2-18")

Temp_1_15 <- aggregate(Temp~Date, data=Data_1_15, FUN=function(x) c(mean=mean(x)))
Temp_1_15$System <- rep("1-15",nrow(Temp_1_15))

Temp_2_15 <- aggregate(Temp~Date, data=Data_2_15, FUN=function(x) c(mean=mean(x)))
Temp_2_15$System <- rep("2-15",nrow(Temp_2_15))

Temp_1_18 <- aggregate(Temp~Date, data=Data_1_18, FUN=function(x) c(mean=mean(x)))
Temp_1_18$System <- rep("1-18",nrow(Temp_1_18))

Temp_2_18 <- aggregate(Temp~Date, data=Data_2_18, FUN=function(x) c(mean=mean(x)))
Temp_2_18$System <- rep("2-18",nrow(Temp_2_18))

Temp_data <- rbind(Temp_1_15, Temp_2_15, Temp_1_18, Temp_2_18)

#selecting for most recent days
Temp_data_recent <- subset(Temp_data, Temp_data$Date >= 58)

windows()
ggplot(data=Temp_data_recent, aes(x=factor(Date), y=Temp, group=System)) +
  geom_line() +
  scale_fill_manual(values=c("gray", "black", "blue", "red", "green")) +
  labs(y=expression("temp (°C)"), x="", fill="") +
  theme(axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=9), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14))

######################################################################
#FULL LIFE HISTORY TEMPS

full_temps <- read.delim("full_life_history_temps.txt", header=T)

egg_1_6 <- subset(full_temps, Location=="egg_chamber" & System=="1_6")
egg_7_12 <- subset(full_temps, Location=="egg_chamber" & System=="7_12")
egg_13_18 <- subset(full_temps, Location=="egg_chamber" & System=="13_18")
egg_Refuge <- subset(full_temps, Location=="egg_chamber" & System=="Refuge")

egg_Temp_1_6 <- aggregate(Temp~Dpf, data=egg_1_6, FUN=function(x) c(mean=mean(x)))
egg_Temp_1_6$System <- rep("1_6",nrow(egg_Temp_1_6))

egg_Temp_7_12 <- aggregate(Temp~Dpf, data=egg_7_12, FUN=function(x) c(mean=mean(x)))
egg_Temp_7_12$System <- rep("7_12",nrow(egg_Temp_7_12))

egg_Temp_13_18 <- aggregate(Temp~Dpf, data=egg_13_18, FUN=function(x) c(mean=mean(x)))
egg_Temp_13_18$System <- rep("13_18",nrow(egg_Temp_13_18))

egg_Temp_Refuge <- aggregate(Temp~Dpf, data=egg_Refuge, FUN=function(x) c(mean=mean(x)))
egg_Temp_Refuge$System <- rep("Refuge",nrow(egg_Temp_Refuge))

egg_Temp_data <- rbind(egg_Temp_1_6, egg_Temp_7_12, egg_Temp_13_18, egg_Temp_Refuge)
egg_Temp_data$Location <- rep("egg_chamber",nrow(egg_Temp_data))




tray <- subset(full_temps, Location=="tray")

tray_A <- subset(full_temps, System=="A")
tray_B <- subset(full_temps, System=="B")

tray_Temp_A <- aggregate(Temp~Dpf, data=tray_A, FUN=function(x) c(mean=mean(x)))
tray_Temp_A$System <- rep("A",nrow(tray_Temp_A))

tray_Temp_B <- aggregate(Temp~Dpf, data=tray_B, FUN=function(x) c(mean=mean(x)))
tray_Temp_B$System <- rep("B",nrow(tray_Temp_B))

tray_Temp_data <- rbind(tray_Temp_A, tray_Temp_B)
tray_Temp_data$Location <- rep("tray",nrow(tray_Temp_data))




tank <- subset(full_temps, Location=="tank")

Data_1_15 <- subset(tank, System=="1-15")
Data_2_15 <- subset(tank, System=="2-15")
Data_1_18 <- subset(tank, System=="1-18")
Data_2_18 <- subset(tank, System=="2-18")

Temp_1_15 <- aggregate(Temp~Dpf, data=Data_1_15, FUN=function(x) c(mean=mean(x)))
Temp_1_15$System <- rep("1-15",nrow(Temp_1_15))

Temp_2_15 <- aggregate(Temp~Dpf, data=Data_2_15, FUN=function(x) c(mean=mean(x)))
Temp_2_15$System <- rep("2-15",nrow(Temp_2_15))

Temp_1_18 <- aggregate(Temp~Dpf, data=Data_1_18, FUN=function(x) c(mean=mean(x)))
Temp_1_18$System <- rep("1-18",nrow(Temp_1_18))

Temp_2_18 <- aggregate(Temp~Dpf, data=Data_2_18, FUN=function(x) c(mean=mean(x)))
Temp_2_18$System <- rep("2-18",nrow(Temp_2_18))

tank_Temp_data <- rbind(Temp_1_15, Temp_2_15, Temp_1_18, Temp_2_18)
tank_Temp_data$Location <- rep("tank",nrow(tank_Temp_data))


all_life_temps <- rbind(egg_Temp_data, tray_Temp_data, tank_Temp_data)


windows()
ggplot(data=all_life_temps, aes(x=factor(Dpf), y=Temp, group=System, color=System)) +
  geom_line() +
  scale_color_manual(values=c("blue", "red", "red", "blue", "blue", "red", "blue", "red", "blue", "red")) +
  labs(y=expression("temp (°C)"), x="", fill="") +
  theme(axis.line = element_line(colour = "black"),
        plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        axis.text = element_text(size=6), axis.title = element_text(size=14), legend.text = element_text(size=14), 
        strip.background = element_blank(), strip.placement = "outside", strip.text = element_text(size=14)) +
  scale_y_continuous(breaks =seq(13,22,0.5))
