
library("dplyr")
library(ggrepel)
library(tidyr)
library("cowplot")
library("data.table")

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign")

## read in parent and family files
parents <- read.delim2("Parent_names_match.txt", header = T)
family <- read.delim2("Family_IDs.txt", header = T)
family <- family[,c(1:7)] #had to remove DI because one male was used in multiple crosses and was causing repeated rows

##Merge parents with family key
##rename for easier merging even though dams and sires are in same file
colnames(parents) <- c("Sire_Seq_ID", "Sire.1")
family_ID <- merge(family, parents, by="Sire.1", all.x=T)
valid_family <- family_ID[,2]

colnames(parents) <- c("Dam_Seq_ID", "Dam.1")
family_ID2 <- merge(family_ID, parents, by="Dam.1", all.x=T)

colnames(parents) <- c("Seq_ID", "Spawn_ID")

Sire <- family_ID2[,c("Sire.1", "Sire.2", "Sire_Seq_ID")]
Dam <- family_ID2[,c("Dam.1", "Dam.2", "Dam_Seq_ID")]


## Read in pedigree and offspring seq info
offspring <- read.delim2("Offspring_names_match.txt", header = T)

## Choose one pedigree below, read in pedigree.top first, remove false families, then read in .pedigree to replace false families
pedigree <- read.table("Smelt_AlphaAssigned_short_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_short_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xaa_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xaa_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xab_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xab_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xac_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xac_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xad_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xad_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xae_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xae_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xaf_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xaf_imputed.pedigree.top", header=F)

pedigree <- read.table("Smelt_AlphaAssigned_xag_imputed.pedigree", header=F)
pedigree <- read.table("Smelt_AlphaAssigned_xag_imputed.pedigree.top", header=F)


colnames(pedigree) <- c("Seq_ID", "Sire_Seq_ID", "Dam_Seq_ID")

## for pedigrees xaa to xaf, remove parents from pedigree
pedigree <- pedigree[87:586,]
## for pedigrees xag, remove parents from pedigree
pedigree <- pedigree[87:171,]

pedigree_meta <- merge(pedigree, offspring, by="Seq_ID", all.x = T)
pedigree_meta <- pedigree_meta[,c(1,2,3,4,7)]

pedigree_meta1 <- merge(pedigree_meta, Sire, by="Sire_Seq_ID", all.x=T)
pedigree_meta2 <- merge(pedigree_meta1, Dam, by="Dam_Seq_ID", all.x=T)
pedigree3 <- unique(pedigree_meta2)

##If read in pedigree with missing data:
pedigree4 <- na.omit(pedigree3)

pedigree4$AAFam <- paste(pedigree4$Dam.2, pedigree4$Sire.2, sep="")
pedigree4_AAFam <- pedigree4[,10]

fam <- factor(pedigree4_AAFam %in% valid_family)

table(fam) 
fam.table <- data.table(fam)

pedigree4_AAFam2 <- cbind(pedigree4, fam.table)

## keep correct families
true_fam <- subset(pedigree4_AAFam2, fam =="TRUE")

##extract false families from top pedigree and match with parent assignment from normal pedigree
false_fam <- subset(pedigree4_AAFam2, fam =="FALSE")

pedigree_norm <- read.table("Smelt_AlphaAssigned_xaa_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xab_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xac_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xad_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xae_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xaf_imputed.pedigree", header=F)
pedigree_norm <- read.table("Smelt_AlphaAssigned_xag_imputed.pedigree", header=F)

colnames(pedigree_norm) <- c("Seq_ID", "Sire_Seq_ID", "Dam_Seq_ID")
false_pedigree_norm <- merge(false_fam, pedigree_norm, by="Seq_ID")
false_pedigree_norm <- false_pedigree_norm[,c(1,4,5,12,13)]

colnames(false_pedigree_norm) <- c("Seq_ID", "SampleName.Fish_ID.", "Treatment", "Sire_Seq_ID","Dam_Seq_ID")

false_pedigree_meta1 <- merge(false_pedigree_norm, Sire, by="Sire_Seq_ID", all.x=T)
false_pedigree_meta2 <- merge(false_pedigree_meta1, Dam, by="Dam_Seq_ID", all.x=T)
false_pedigree3 <- unique(false_pedigree_meta2)
false_pedigree3$AAFam <- paste(false_pedigree3$Dam.2, false_pedigree3$Sire.2, sep="")

##merge true and false back together and save
true_fam <- true_fam[,c(1:10)]
final_meta_pedigree <- rbind(true_fam, false_pedigree3)

write.table(final_meta_pedigree, file = "xaa_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xab_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xac_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xad_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xae_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xaf_pedigree_meta.txt", quote = F, sep = "\t")
write.table(final_meta_pedigree, file = "xag_pedigree_meta.txt", quote = F, sep = "\t")


## read in all final pedigree files and combine
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/AlphaAssign")

xaa <- read.table("xaa_pedigree_meta.txt", header=T)
xab <- read.table("xab_pedigree_meta.txt", header=T)
xac <- read.table("xac_pedigree_meta.txt", header=T)
xad <- read.table("xad_pedigree_meta.txt", header=T)
xae <- read.table("xae_pedigree_meta.txt", header=T)
xaf <- read.table("xaf_pedigree_meta.txt", header=T)
xag <- read.table("xag_pedigree_meta.txt", header=T)

all_pedigrees <- rbind(xaa, xab, xac, xad, xae, xaf, xag)
DI <- read.delim2("Family_DIs.txt", header=T)

all_pedigree_DI <- merge(all_pedigrees, DI, by="AAFam", all.x = T)

write.table(all_pedigree_DI, file = "all_pedigree_meta_DI.txt", quote = F, sep = "\t")
