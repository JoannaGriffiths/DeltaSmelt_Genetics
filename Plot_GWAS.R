
library("qqman")
library("dplyr")
library("tidyr")

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/GWAS")

## CTM
gwas <- read.delim2("all_chr_offspring_lmm4.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_15only.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_18only.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_fixed.assoc.txt", header=T) #THIS ONE
gwas <- read.delim2("all_chr_offspring_LDfiltered_lmm4.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_fixed_system.assoc.txt")


gwas <- read.delim2("all_chr_offspring_lmm4_fixed_gxe.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_fixed_gxe_tank.assoc.txt", header=T)
gwas <- read.delim2("all_chr_offspring_lmm4_fixed_gxe_system.assoc.txt", header=T) #THIS ONE

## DI
gwas <- read.delim2("DI_all_chr_offspring_lmm4_fixed_tank_reartemp.assoc.txt", header=T)

## FL
gwas <- read.delim2("FL_all_chr_offspring_lmm4_fixed_tank_reartemp.assoc.txt", header=T)
gwas <- read.delim2("gxe_FL_all_chr_offspring_lmm4.assoc.txt", header=T)

## CTM-DI PLACO results
gwas <- read.delim2("placo_results_correlated.txt", header=T)
snps <- read.delim2("gwas_ctm_DI_shared_snps", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(9,6,3)]


## CTM-gxe PLACO results
gwas <- read.delim2("placo_results_gxe_correlated.txt", header=T)
snps <- read.delim2("gwas_ctmax_gxe_ldsc", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(14,6,3)]

## check there are no z-scores above 80 in both datasets before running placo+
gwas <- read.delim2("gwas_summary_gxe.txt", header=F)
gwas$V2 <- as.numeric(gwas$V2)
gwas$V4 <- as.numeric(gwas$V4)


## FL-FLgxe PLACO results
gwas <- read.delim2("placo_results_FL_correlated.txt", header=T)
snps <- read.delim2("gwas_FL_FLgxe_shared_snps", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(9,6,3)]


## ctm-FL PLACO results
gwas <- read.delim2("placo_results_ctm_FL_correlated.txt", header=T)
snps <- read.delim2("gwas_ctm_FL_shared_snps", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(9,6,3)]


## ctmgxe-FLgxe PLACO results
gwas <- read.delim2("placo_results_ctmgxe_FLgxe_correlated.txt", header=T)
snps <- read.delim2("gwas_ctmgxe_FLgxe_shared_snps", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(9,6,3)]


## ctmgxe-DI PLACO results
gwas <- read.delim2("placo_results_ctmgxe_DI_correlated.txt", header=T)
snps <- read.delim2("gwas_ctmgxe_DI_shared_snps", header=T)
snps <- separate(data = snps, col = SNP, into = c("chr1", "chr2", "SNP", "Allele1", "Allele2"), sep = "_")
snps$chr <- paste(snps$chr1, snps$chr2, sep="_")
gwas_snps <- cbind(gwas, snps)
gwas <- gwas_snps[c(9,6,3)]

##########################################################################

gwas <- gwas[,c(1,2,3,13)]

##for making new chr name file:
#chr <- data.frame(unique(gwas$chr))

########get top SNPs for GO
gwas$p_wald <- as.numeric(gwas$p_wald)
gwas$log_p <- -log10(gwas$p_wald)

##find top 5% and prepare for LDannot
top5 <- quantile(gwas$p_wald, 0.05) #0.0493
top1 <- quantile(gwas$p_wald, 0.01) #0.00965

top5list <- subset(gwas, p_wald < top5)
top1list <- subset(gwas, p_wald < top1)
tippytop <- subset(gwas, p_wald < 0.0001) #247 snps top 0.02%
top_log4 <- subset(gwas, log_p > 4) #247 snps, 0.02%
top_log3.5 <- subset(gwas, log_p > 3.5) #570 snps, 0.04%
bf_sig <- subset(gwas, p_wald < (0.05/1223476))

top_300 <- subset(gwas, log_p > 6.04845) #300 SNPs for CTMax GxE GWAS
top_300 <- subset(gwas, log_p > 3.878) #300 SNPs for CTMax GWAS
top_300 <- subset(gwas, log_p > 6.122) #300 SNPs for DI GWAS
top_300 <- subset(gwas, log_p > 3.953) #299 SNPs for FL GWAS
top_300 <- subset(gwas, log_p > 8.5) #300 SNPs for FL GxE GWAS

top5list <- top5list[,c(1,3)]
top1list <- top1list[,c(1,3)]
top_log3.5 <- top_log3.5[,c(1,3)]
tippytop <- tippytop[,c(1,3)]
bf_sig <- bf_sig[,c(1,3)]
top_300 <- top_300[,c(1,3)]

colnames(top5list) = c("chr", "SNP")
colnames(top1list) = c("chr", "SNP")
colnames(top_log3.5) = c("chr", "SNP")
colnames(tippytop) = c("chr", "SNP")
colnames(bf_sig) = c("chr", "SNP")
colnames(top_300) = c("chr", "SNP")

top5list$chr_SNP <- paste(top5list$chr, top5list$SNP, sep = "_")
top1list$chr_SNP <- paste(top1list$chr, top1list$SNP, sep = "_")
top_log3.5$chr_SNP <- paste(top_log3.5$chr, top_log3.5$SNP, sep = "_")
tippytop$chr_SNP <- paste(tippytop$chr, tippytop$SNP, sep = "_")
bf_sig$chr_SNP <- paste(bf_sig$chr, bf_sig$SNP, sep = "_")
top_300$chr_SNP <- paste(top_300$chr, top_300$SNP, sep = "_")


write.table(top5list,file = "gwas_top5list_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top1list,file = "gwas_top1list_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_log3.5,file = "gwas_top_log3.5_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_log3.5,file = "gwas_top_log3.5_FL_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(tippytop,file = "gwas_tippytop_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(bf_sig,file = "gwas_bf_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(bf_sig,file = "gwas_genomewide_FL_gxe_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(top5list,file = "gwas_top5list_for_LDannot_fixed_gxe", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top1list,file = "gwas_top1list_for_LDannot_fixed_gxe", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(tippytop,file = "gwas_tippytop_for_LDannot_fixed_gxe", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(bf_sig,file = "gwas_bf_for_LDannot_fixed_gxe", quote = F,row.names = F,col.names = T, sep = "\t")


write.table(top_300,file = "gwas_top300_for_LDannot_ctmax_gxe", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_300,file = "gwas_top300_for_LDannot_ctmax", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_300,file = "gwas_top300_for_LDannot_DI", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_300,file = "gwas_top300_for_LDannot_FL", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(top_300,file = "gwas_top300_for_LDannot_FL_gxe", quote = F,row.names = F,col.names = T, sep = "\t")

###############################
#### Changes made to LD-annot code
############################
L328 old: dics[z] = [(int(j[1]) - mean) , (int(j[1]) + mean)]
L335 and L336:
  up = int(max(dics[i]))
down = int(min(dics[i]))
L340 old: if int(l[1]) < down < int(l[2]) or int(l[1]) < up < int(l[2]) or down < int(l[1]) < up or down < int(l[2]) < up :

## run on the command line for input files before running LD-annot
for all files:
`:%s/NC_/NC/g`
`:%s/NW_/NW/g`


##################################
#######  Fix Chromosome names
###############################

names <- read.delim2("DS_chr_num_names.txt", header=F)
colnames(names) <- c("chr_num", "chr")

gwas2 <- merge(gwas, names, by="chr", all.x = T)
rm(gwas)

gwas2$chr <- as.numeric(gwas2$chr)
gwas2$ps <- as.numeric(gwas2$ps)
gwas2$p_wald <- as.numeric(gwas2$p_wald)
gwas2$P_wald_log <- (-log10(gwas2$p_wald))
gwas2$p.placo.plus <- as.numeric(gwas2$p.placo.plus)
gwas2$P_wald_log <- (-log10(gwas2$p.placo.plus)) #for placo results
gwas2$SNP <- as.numeric(gwas2$SNP)
gwas2$snp <- paste(gwas2$chr_num, gwas2$ps, sep = "_")

##############################################
### bonferroni corrections for each chromosome
###################################################

snps_chr <- gwas2 %>%
  group_by(chr_num) %>%
  count(chr_num)

t_snps_chr <- data.frame(t(snps_chr))
#header_names <- as.character(t_snps_chr[1, ])
#colnames(t_snps_chr) <- header_names
t_snps_chr <- t_snps_chr[-1, ]

#test code with 1 chr
sig_chr1 <- subset(gwas2, chr_num ==1 & P_wald_log > (-log10(0.05/t_snps_chr$X1)))

##for all chr, with all deets
for (i in 1:316) {
  assign(
    paste0("sig_chr", i),
    subset(
      gwas2,
      chr_num == i & P_wald_log > (-log10(0.05 / t_snps_chr[[paste0("X", i)]]))
    )
  )
}

##for all chr, with only chr name
for (i in 1:316) {
  assign(
    paste0("sig_chr", i),
    subset(
      gwas2,
      chr_num == i & P_wald_log > (-log10(0.05 / t_snps_chr[[paste0("X", i)]]))
    #)[, 7] #for gemma
    )[, 5] #for placo
  )
}


sig_vec_0.05 <- c()
for (i in 1:316) {
  sig_vec_0.05 <- c(sig_vec_0.05, get(paste0("sig_chr", i)))
}




for (i in 1:316) {
  assign(
    paste0("sig_chr", i),
    subset(
      gwas2,
      chr_num == i & P_wald_log > (-log10(0.1 / t_snps_chr[[paste0("X", i)]]))
    )[, 7]
  )
}
sig_vec_0.1 <- c()
for (i in 1:316) {
  sig_vec_0.1 <- c(sig_vec_0.1, get(paste0("sig_chr", i)))
}

## CTM
save(sig_vec_0.05, sig_vec_0.1, file = "chr_sig.RData")
save(sig_vec_0.05, sig_vec_0.1, file = "chr_sig_gxe.RData")

## DI
save(sig_vec_0.05, sig_vec_0.1, file = "DI_chr_sig.RData")

## FL
save(sig_vec_0.05, sig_vec_0.1, file = "FL_chr_sig.RData")
save(sig_vec_0.05, sig_vec_0.1, file = "FL_gxe_chr_sig.RData")

## CTM-FL
save(sig_vec_0.05, sig_vec_0.1, file = "CTM-FL_chr_sig.RData")

## CTM-DI
save(sig_vec_0.05, sig_vec_0.1, file = "CTM-DI_chr_sig.RData")


##########################
##### Create LD-annot input, don't run before manhattan plots
##########################
load("chr_sig.RData")
load("chr_sig_gxe.RData")
load("DI_chr_sig.RData")
load("FL_chr_sig.RData")
load("CTM-FL_chr_sig.RData")
load("CTM-DI_chr_sig.RData")
load("FL_gxe_chr_sig.RData")

library("tidyr")

sig_vec_0.05 <- data.frame(sig_vec_0.05)
sig_vec_0.05 <- separate(data = sig_vec_0.05, col = sig_vec_0.05, into = c("chr_num", "SNP"), sep = "_")
sig_vec_0.05_names <- merge(sig_vec_0.05, names, by="chr_num")
sig_vec_0.05_names$chr_SNP <- paste(sig_vec_0.05_names$chr, sig_vec_0.05_names$SNP, sep = "_")

sig_vec_0.05_names$chr_num <- as.numeric(sig_vec_0.05_names$chr_num)
#sig_vec_0.05_names <- subset(sig_vec_0.05_names, chr_num < 27)
sig_vec_0.05_names <- sig_vec_0.05_names[c(3,2,4)]


write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_fixed_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_gxe_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_gxe_for_LDannot_fixed_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_DI", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_DI_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_FL", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_FL_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_CTM-DI", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.05_names,file = "gwas_chr_sig_0.05_for_LDannot_CTM-FL", quote = F,row.names = F,col.names = T, sep = "\t")


sig_vec_0.1 <- data.frame(sig_vec_0.1)
sig_vec_0.1 <- separate(data = sig_vec_0.1, col = sig_vec_0.1, into = c("chr_num", "SNP"), sep = "_")
sig_vec_0.1_names <- merge(sig_vec_0.1, names, by="chr_num")
sig_vec_0.1_names$chr_SNP <- paste(sig_vec_0.1_names$chr, sig_vec_0.1_names$SNP, sep = "_")

sig_vec_0.1_names$chr_num <- as.numeric(sig_vec_0.1_names$chr_num)
#sig_vec_0.1_names <- subset(sig_vec_0.1_names, chr_num < 27)
sig_vec_0.1_names <- sig_vec_0.1_names[c(3,2,4)]

write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_for_LDannot_fixed_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_gxe_for_LDannot_fixed", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_gxe_for_LDannot_fixed_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_for_LDannot_DI", quote = F,row.names = F,col.names = T, sep = "\t")
write.table(sig_vec_0.1_names,file = "gwas_chr_sig_0.1_for_LDannot_DI_largeChr", quote = F,row.names = F,col.names = T, sep = "\t")

gwas$p_wald <- as.numeric(gwas$p_wald)
gwas_genomewide_sig <- subset(gwas, p_wald < (0.05/1224047)) ##for DI only
gwas_genomewide_sig$snp <- paste(gwas_genomewide_sig$chr, gwas_genomewide_sig$ps, sep="_")
gwas_genomewide_sig <- gwas_genomewide_sig[c("chr", "ps", "snp")]

write.table(gwas_genomewide_sig,file = "gwas_genomewide_for_LDannot_DI", quote = F,row.names = F,col.names = T, sep = "\t")

## on the command line before running LD-annot
for all files:
  `:%s/NC_/NC/g`
`:%s/NW_/NW/g`

#######################
####### sig overlap
######################
ctm_sig_vec_0.05 <- sig_vec_0.05
ctm_sig_vec_0.1 <- sig_vec_0.1

gxe_sig_vec_0.05 <- sig_vec_0.05
gxe_sig_vec_0.1 <- sig_vec_0.1

DI_sig_vec_0.05 <- sig_vec_0.05
DI_sig_vec_0.1 <- sig_vec_0.1


# compare the chr lists of cands & all genes and turn into factor as to whether they overlap or not
overlap <- factor(as.integer(ctm_sig_vec_0.05 %in% gxe_sig_vec_0.05))
table(overlap) #only 2 snps overlap

overlap <- factor(as.integer(ctm_sig_vec_0.1 %in% gxe_sig_vec_0.05))
table(overlap) #only 2 snps overlap

overlap2 <- factor(as.integer(ctm_sig_vec_0.05 %in% DI_sig_vec_0.05))
table(overlap2) #0 overlap

overlap3 <- factor(as.integer(ctm_sig_vec_0.1 %in% DI_sig_vec_0.05))
table(overlap3) #0 overlap

overlap4 <- factor(as.integer(gxe_sig_vec_0.05 %in% DI_sig_vec_0.05))
table(overlap4) #0 overlap



################################
### Manhattan Plots
################################
## CTM
windows()
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", suggestiveline = -log10(0.05/1223476), genomewideline = 4, highlight = sig_vec_0.05, logp = T) #7.39 #OG
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", suggestiveline = -log10(0.05/1223476), genomewideline = 3.5, highlight = sig_vec_0.05, logp = T) #7.39 #OG

manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1223476), genomewideline = F, logp = T)
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.1, suggestiveline = -log10(0.05/1223476), genomewideline = -log10(0.1/1223476), logp = T)

## DI
windows()
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1224047), genomewideline = F, logp = T)
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.1, suggestiveline = -log10(0.05/1224047), genomewideline = -log10(0.1/1224047), logp = T)

## FL
windows()
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1223467), genomewideline = 3.5, logp = T)
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1223467), genomewideline = F, logp = T)

## CTM-FL
windows()
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1223467), genomewideline=F, logp = T)
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.1, suggestiveline = -log10(0.05/1223467), logp = T)
manhattan(subset(gwas2, chr_num==8), chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1223467), logp = T)

## CTM-DI
windows()
manhattan(gwas2, chr="chr_num", snp="snp", bp="ps", p="p_wald", highlight = sig_vec_0.05, suggestiveline = -log10(0.05/1224047), genomewideline = F, logp = T)

## PLACO CTM-DI
gwas2$SNP <- as.numeric(gwas2$SNP)
gwas2$p.placo.plus <- as.numeric(gwas2$p.placo.plus)
gwas2$p.placo.plus <- abs(gwas2$p.placo.plus)
windows()
manhattan(gwas2, chr="chr_num", snp="SNP", bp="SNP", p="p.placo.plus", suggestiveline = -log10(0.05/1219600), logp = T, )


## PLACO CTM-gxe, FL-FLgxe, ctm-FL, ctmgxe-FLgxe
gwas2$SNP <- as.numeric(gwas2$SNP)
gwas2$p.placo.plus <- as.numeric(gwas2$p.placo.plus)
gwas2$p.placo.plus <- abs(gwas2$p.placo.plus)
gwas3 <- subset(gwas2, p.placo.plus > 0)
windows()
manhattan(gwas3, chr="chr_num", snp="SNP", bp="SNP", p="p.placo.plus", suggestiveline = -log10(0.05/1223467), logp = T, )


## PLACO CTMgxe-DI
gwas2$SNP <- as.numeric(gwas2$SNP)
gwas2$p.placo.plus <- as.numeric(gwas2$p.placo.plus)
gwas2$p.placo.plus <- abs(gwas2$p.placo.plus)
gwas3 <- subset(gwas2, p.placo.plus > 0)
windows()
manhattan(gwas3, chr="chr_num", snp="SNP", bp="SNP", p="p.placo.plus", suggestiveline = -log10(0.05/1219600), logp = T, )




gwas2$ps <- as.numeric(gwas2$ps)
windows()
manhattan(subset(gwas2, chr== 254038 & ps < 1000000), chr="chr", snp="ps", bp="ps", p="p_wald", suggestiveline = -log10(0.05/1223476), genomewideline = -log10(0.1/1223476), logp = T) #7.39

##LD filtered
windows()
manhattan(gwas2, chr="chr", snp="ps", bp="ps", p="p_wald", suggestiveline = -log10(0.05/615336), genomewideline = -log10(0.1/615336), logp = T) #7.09





###########################
### R script to do GO functional enrichment
###########################

setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/GWAS")

#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("biomaRt")
#BiocManager::install("topGO")

library("biomaRt")
library(topGO)

listMarts()

# look at top 10 databases      
head(biomaRt::listMarts(host = "http://www.ensembl.org/"), 10)  

#collect gene names from biomart
mart=biomaRt::useMart('ENSEMBL_MART_ENSEMBL', host = "http://www.ensembl.org/")
listDatasets(mart)

#using zebrafish genes for now
#Make sure the host includes 'www' otherwise fetching the attributes below doesn't work!
mart2 <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                          dataset = "drerio_gene_ensembl",
                          host = 'http://www.ensembl.org/',
                          verbose=TRUE)


#head(biomaRt::listAttributes(mart=mart2), 500)
## Get ensembl gene ids and GO terms
GTOGO <- biomaRt::getBM(attributes = c("go_id","external_gene_name"), 
                        mart = mart2,
                        verbose=TRUE)

#examine result
head(GTOGO)

#Remove blank entries
GTOGO <- GTOGO[GTOGO$external_gene_name != '',]

# convert from table format to list format
geneID2GO <- by(GTOGO$go_id,
                GTOGO$external_gene_name,
                function(x) as.character(x))

#examine result
head(geneID2GO)



#select genes of interest from the gene list
library(dplyr)
library(tidyr)

## select one output file
#LDannot <- read.delim2("gwas_tippytop_LDannot_output")
#LDannot <- read.delim2("gwas_tippytop_LDannot_output_fixed") #140, 83 unique
#LDannot <- read.delim2("gwas_top1list_LDannot_output_fixed") #7677, 3212 unique
LDannot <- read.delim2("gwas_log3.5_LDannot_output") #363, 209 unique

LDannot <- read.delim2("gwas_bf_LDannot_output_fixed_gxe") #3
LDannot <- read.delim2("gwas_top1list_LDannot_output_fixed_gxe") #7689

LDannot <- read.delim2("gwas_chr_sig_0.05_LDannot_output") #38, but 26 unique, no enrichment
LDannot <- read.delim2("gwas_chr_sig_0.05_largeChr_LDannot_output") #2, 1 unique, no overlap with zebra fish gene names
LDannot <- read.delim2("gwas_chr_sig_0.1_LDannot_output") #60, but 38 unique, no enrichment
LDannot <- read.delim2("gwas_chr_sig_0.1_largeChr_LDannot_output") #2, 1 unique, no overlap with zebra fish gene names
LDannot <- read.delim2("gwas_top300_ctmax_LDannot_output") #300 snps, 105 unique genes

LDannot <- read.delim2("gwas_chr_sig_gxe_0.05_LDannot_output") #507, but 247 unique
LDannot <- read.delim2("gwas_chr_sig_gxe_0.05_largeChr_LDannot_output") #150, 84 unique
LDannot <- read.delim2("gwas_chr_sig_gxe_0.1_LDannot_output")  #803, 361 unique
LDannot <- read.delim2("gwas_chr_sig_gxe_0.1_largeChr_LDannot_output")  #244, 128 unique
LDannot <- read.delim2("gwas_top300_ctmax_gxe_LDannot_output") #300 snps, 90 unique

## DI
LDannot <- read.delim2("gwas_genomewide_LDannot_output_DI") #110, 69 unique, 5 functions enriched
LDannot <- read.delim2("gwas_chr_sig_0.05_DI_LDannot_output") #281, 192 unique, 18 functions enriched
LDannot <- read.delim2("gwas_chr_sig_0.1_DI_LDannot_output") #384, 238 unique, 25 functions enriched
LDannot <- read.delim2("gwas_chr_sig_0.05_DI_largeChr_LDannot_output") #224, 165 unique, 12 function enriched
LDannot <- read.delim2("gwas_chr_sig_0.1_DI_largeChr_LDannot_output") #237, 199 unique, 17 functions enriched
LDannot <- read.delim2("gwas_top300_LDannot_output_DI") #300 snps, 144 unique genes


## FL
LDannot <- read.delim2("gwas_chr_sig_0.05_FL_LDannot_output") #36, 15 unique genes
LDannot <- read.delim2("gwas_chr_sig_0.05_FL_LargeChr_LDannot_output") #21, 6 unique genes
LDannot <- read.delim2("gwas_log3.5_FL_LDannot_output") #393, 192 unique genes
LDannot <- read.delim2("gwas_genomewide_FL_gxe_LDannot_output") #514, 268 unique
LDannot <- read.delim2("gwas_top300_for_FL_LDannot_output") #299snps, 98 unique
LDannot <- read.delim2("gwas_top300_for_FL_gxe_LDannot_output") #300snps, 91 unique

## CTM-DI (not actually gxe)
LDannot <- read.delim2("gwas_chr_sig_gxe_0.05_CTM-DI_LDannot_output") #219, 148 unique genes

## CTM-FL
LDannot <- read.delim2("gwas_chr_sig_0.05_CTM-FL_LDannot_output") #50, 31 unique genes

## FST
LDannot <- read.delim2("../Fst/Fst_outflank_LDannot_output")  #5, 5 unique, no enrichment, but also 1 gene overlap with zebra fish gene names
LDannot <- read.delim2("../Fst/Fst_top0.3_LDannot_output") #3926, but 210 unique

## Background genes
LDannot <- read.delim2("vcf_snps_gene_overlap") #23835 unique genes for possible overlap with snps in vcf


############# get snp location relative to gene
library(tidyr)
LDannot2 <- LDannot %>%
  separate(col = SNP, into = c("CHR", "pos"), sep = "_")
LDannot2$gene_length <- LDannot2$gene_end - LDannot2$gene_start
LDannot2$pos <- as.numeric(LDannot2$pos)
LDannot2$snp_loc <- LDannot2$gene_start - LDannot2$pos
LDannot2$snp_loc2 <- LDannot2$region_start - LDannot2$pos

LDannot3 <- subset(LDannot2, snp_loc > -10000)
LDannot <- LDannot3
############

gene_ids <- separate(data = LDannot, col=annotation, into= c("gene_LOC", "GeneID"), sep = ";", remove=T)
gene_ids$gene_LOC <- gsub("ID=gene-","",as.character(gene_ids$gene_LOC))


int.genes.df <- data.frame(gene_ids[,("gene_LOC")])
colnames(int.genes.df) <- "gene_id"
int.genes.df <- unique(int.genes.df) ##new addition


############################ only run below for overlap of gene ids between traits
ctm0.05_gene_ids <- int.genes.df
ctmlog3.5_gene_ids <- int.genes.df
gxe_gene_ids <- int.genes.df
DI_gene_ids <- int.genes.df
Fst_gene_ids <- int.genes.df
FLlog3.5_gene_ids <- int.genes.df
FL_gxe_gene_ids <- int.genes.df
ctm_DI_gene_ids <- int.genes.df
ctm_FL_gene_ids <- int.genes.df

ctm_top300_gene_ids <- int.genes.df
ctm_gxe_top300_gene_ids <- int.genes.df
DI_top300_gene_ids <- int.genes.df
FL_top300_gene_ids <- int.genes.df
FL_gxe_top300_gene_ids <- int.genes.df

# compare the chr lists of cands & all genes and turn into factor as to whether they overlap or not
setwd("C:/Users/joann/OneDrive/Documents/UCDavis/Whitehead_lab/Smelt_sequencing/2021_spawning/Analyses/GWAS")
load("gwas_gene_ids.RData")
overlap <- merge(ctmlog3.5_gene_ids, gxe_gene_ids, by="gene_id")
table(overlap) #17 overlap

overlap <- merge(ctm0.05_gene_ids, gxe_gene_ids, by="gene_id")
table(overlap) #3 overlap

overlap <- merge(DI_gene_ids, ctm0.05_gene_ids, by="gene_id")
table(overlap) #2 overlap

overlap <- merge(DI_gene_ids, ctmlog3.5_gene_ids, by="gene_id")
table(overlap) #6 overlap

overlap <- merge(DI_gene_ids, gxe_gene_ids, by="gene_id")
table(overlap) #3 overlap

overlap <- merge(Fst_gene_ids, ctm0.05_gene_ids, by="gene_id")
table(overlap) #1 overlap

overlap <- merge(Fst_gene_ids, ctmlog3.5_gene_ids, by="gene_id")
table(overlap) #7 overlap

overlap <- merge(Fst_gene_ids, gxe_gene_ids, by="gene_id")
table(overlap) #2 overlap

overlap <- merge(Fst_gene_ids, DI_gene_ids, by="gene_id")
table(overlap) #6 overlap

overlap <- merge(ctm0.05_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #0

overlap <- merge(ctmlog3.5_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #8

overlap <- merge(ctm0.05_gene_ids, FL_gxe_gene_ids, by="gene_id")
table(overlap) #1

overlap <- merge(ctmlog3.5_gene_ids, FL_gxe_gene_ids, by="gene_id")
table(overlap) #9

overlap <- merge(gxe_gene_ids, FL_gxe_gene_ids, by="gene_id")
table(overlap) #3

overlap <- merge(DI_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #7

overlap <- merge(Fst_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #5

overlap <- merge(Fst_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #5

overlap <- merge(ctm_FL_gene_ids, FLlog3.5_gene_ids, by="gene_id")
table(overlap) #8

overlap <- merge(FLlog3.5_gene_ids, FL_gxe_gene_ids, by="gene_id")
table(overlap) #10

###
overlap <- merge(ctm_top300_gene_ids, ctm_gxe_top300_gene_ids, by="gene_id")
table(overlap) #6

overlap <- merge(ctm_top300_gene_ids, DI_top300_gene_ids, by="gene_id")
table(overlap) #2

overlap <- merge(ctm_gxe_top300_gene_ids, DI_top300_gene_ids, by="gene_id")
table(overlap) #1

overlap <- merge(FL_gxe_top300_gene_ids, FL_top300_gene_ids, by="gene_id")
table(overlap) #1

overlap <- merge(ctm_top300_gene_ids, FL_top300_gene_ids, by="gene_id")
table(overlap) #3

overlap <- merge(ctm_gxe_top300_gene_ids, FL_top300_gene_ids, by="gene_id")
table(overlap) #2

overlap <- merge(ctm_top300_gene_ids, FL_gxe_top300_gene_ids, by="gene_id")
table(overlap) #3

overlap <- merge(ctm_gxe_top300_gene_ids, FL_gxe_top300_gene_ids, by="gene_id")
table(overlap) #2

##
save(Fst_gene_ids, DI_gene_ids, ctm0.05_gene_ids, ctmlog3.5_gene_ids, gxe_gene_ids, FLlog3.5_gene_ids, FL_gxe_gene_ids, ctm_DI_gene_ids, ctm_FL_gene_ids,ctm_top300_gene_ids,ctm_gxe_top300_gene_ids,DI_top300_gene_ids,FL_top300_gene_ids,FL_gxe_top300_gene_ids, file="gwas_gene_ids.RData")



####################################################
## GO functional enrichment continued
#################################################

#keep only unique gene IDs from biomaRt--but which GO term gets chosen, I think it get's chosen in topGO later
all.genes <- sort(unique(as.character(GTOGO$external_gene_name)))

cands.chr <-as.character(int.genes.df$gene_id)

# compare the chr lists of cands & all genes and turn into factor as to whether they overlap or not
int.genes <- factor(as.integer(all.genes %in% cands.chr))
names(int.genes) = all.genes
table(int.genes)

library(data.table)
int.genes.df <- data.table(int.genes)

# #create GO object (for running in topGO!)
library(topGO)
test.stat <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")
test.stat2 <- new("weightCount", testStatistic = GOFisherTest, name = "Fisher test", sigRatio = "ratio")

# GO for Biological Process
go.obj.bp <- new("topGOdata", ontology='BP'
                 , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
                 , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
                 , gene2GO = geneID2GO # this is our object we just created with biomaRt
)

resultFisher <- getSigGroups(go.obj.bp, test.stat)
resultWeight <- getSigGroups(go.obj.bp, test.stat2)

allRes_bp <- GenTable(go.obj.bp, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic",topNodes=length(resultFisher@score),numChar=100)
filtRes.BP <- allRes_bp[allRes_bp$classic<0.05 & allRes_bp$Significant>2,]
filtRes.BP$Ontology <- "BP"

# GO for Molecular Function
go.obj.mf <- new("topGOdata", ontology='MF'
                 , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
                 , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
                 , gene2GO = geneID2GO # this is our object we just created with biomaRt
)

resultFisher <- getSigGroups(go.obj.mf, test.stat)
resultWeight <- getSigGroups(go.obj.mf, test.stat2)
allRes_mf <- GenTable(go.obj.mf, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic",topNodes=length(resultFisher@score),numChar=100)
filtRes.MF <- allRes_mf[allRes_mf$classic<0.05 & allRes_mf$Significant>2,]
filtRes.MF$Ontology <- "MF"

# GO for Cellular Component
go.obj.cc <- new("topGOdata", ontology='CC'
                 , allGenes = int.genes # this is factor containing interesting genes (gene IDs for outliers)
                 , annot = annFUN.gene2GO # we use this function as we supply our own annotations 
                 , gene2GO = geneID2GO # this is our object we just created with biomaRt
)
resultFisher <- getSigGroups(go.obj.cc, test.stat)
resultWeight <- getSigGroups(go.obj.cc, test.stat2)
allRes_cc <- GenTable(go.obj.cc, classic = resultFisher, 
                   weight = resultWeight, orderBy = "weight", 
                   ranksOf = "classic", topNodes=length(resultFisher@score))
filtRes.CC <- allRes_cc[allRes_cc$classic<0.05 & allRes_cc$Significant>2,]
filtRes.CC$Ontology <- "CC"

##Output all three
filt.all <- rbind(filtRes.BP,filtRes.MF,filtRes.CC)

#Write to correspinging input file
write.csv(filt.all,"GO.gwas.tippytop.csv") #10
write.csv(filt.all,"GO.gwas.log3.5.csv") #14, 16 if remove snps more than 10,000 upstream of gene
write.csv(filt.all,"GO.gwas.top1list.csv") #270
write.csv(filt.all,"GO.gwas.top1list.gxe.csv") #74

write.csv(filt.all,"GO.gwas.chr_sig_0.05.gxe.csv") #34
write.csv(filt.all,"GO.gwas.chr_sig_0.05.gxe.largeChr.csv") #27
write.csv(filt.all,"GO.gwas.chr_sig_0.1.gxe.csv") #53
write.csv(filt.all,"GO.gwas.chr_sig_0.1.gxe.largeChr.csv") #41

write.csv(filt.all,"GO.gwas.genomewide.CTM-DI.csv") #9

write.csv(filt.all,"GO.gwas.genomewide.DI.csv") #5
write.csv(filt.all,"GO.gwas.chr_sig_0.05.DI.csv") #18
write.csv(filt.all,"GO.gwas.chr_sig_0.1.DI.csv") #25
write.csv(filt.all,"GO.gwas.chr_sig_0.05.DI.largeChr.csv") #12
write.csv(filt.all,"GO.gwas.chr_sig_0.1.DI.largeChr.csv") #17

write.csv(filt.all,"GO.gwas.log3.5.FL.csv") #13
write.csv(filt.all,"GO.gwas.gxe.FL.csv") #14

write.csv(filt.all,"GO.gwas.Fst_top0.3.csv") #10

