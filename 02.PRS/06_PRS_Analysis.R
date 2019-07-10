####################
### PRS Analysis ###
####################
# In this R script, we will analyze and plot all generated polygenic risk scores

# 1. Load packages (for own R installation - version 3)
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
library(viridis) # for colorblind color palette
library(data.table)
library(ggplot2)
library(corrplot)
library(psych)
library(pheatmap)
library(limma)
library(rcompanion)
library(oligo)
library(RColorBrewer)

# 2. Load data (all different generated PRS)
# load, filter and save phenotype data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_core_gene_summaries_annotated_combat_batch_corrected.RData")
index <- which(is.na(pData(geneSummaries)$C1))  # remove all samples with C1 = NA
geneSummaries <- geneSummaries[, -c(index)]
index <- which(is.na(pData(geneSummaries)$Type_of_Death))  # remove all samples with Type_of_Death = NA
geneSummaries <- geneSummaries[, -c(index)]
dim(geneSummaries)
pData <- pData(geneSummaries)

# load PRS LDpred-inf scores. 
# Remember to add a comma "," at the end of the header line (or remove all commas from the rest of the lines) otherwise read.table or fread can not read in the data correctly. 
# load SCZ2 PRS
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load MDD2 PRS
prs.mdd2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/MDD2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_mdd2_filtered_PRS_LDpred_p1.0000e+00.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load BIP-2018 PRS
prs.bip <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIP-2018/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_filtered_PRS_LDpred_p1.0000e-02.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load BIP + SCZ vs. CONT
prs.bip.scz.vs.ctl <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIP+SCZvsCONT/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_scz_vs_cont_filtered_PRS_LDpred_p1.0000e+00.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load BIP vs SCZ
prs.bip.vs.scz <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIPvsSCZ/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_vs_scz_filtered_PRS_LDpred_p1.0000e-03.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load ADHD-2017
prs.adhd <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/ADHD-2017/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_adhd_filtered_PRS_LDpred_p1.0000e-02.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load ASD
prs.asd <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/ASD/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_asd_filtered_PRS_LDpred_p1.0000e-03.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load educational attainment
prs.ea <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/EA/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_ea_filtered_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load BMI
#prs.bmi <- fread("", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load cross disorder
prs.cross.disorder <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/Cross_Disorder/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_cd_filtered_PRS_LDpred_p1.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load GR Score and filter for our datasets
gr.score <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/03.TRPS/03_TRPS.txt", sep = "\t", header = T, stringsAsFactors = F, select = c("IID", "gr_score"))
gr.score <- gr.score[!startsWith(gr.score$IID, "MD") & !endsWith(gr.score$IID, "ZNP")]
gr.score <- gr.score[gr.score$IID %in% pData$IID]
gr.score <- gr.score[, IID:=as.integer(IID)]

# Add diagnosis column to the PRS
pData$Diagnosis <- as.character(pData$Dx)
pData[pData$Diagnosis != "CTL", "Diagnosis"] <- "CASE"
prs.scz2 <- merge(prs.scz2, pData, by = "IID")
prs.mdd2 <- merge(prs.mdd2, pData, by = "IID")
prs.bip <- merge(prs.bip, pData, by = "IID")
prs.adhd <- merge(prs.adhd, pData, by = "IID")
prs.asd <- merge(prs.asd, pData, by = "IID")
prs.bip.vs.scz <- merge(prs.bip.vs.scz, pData, by = "IID")
prs.bip.scz.vs.ctl <- merge(prs.bip.scz.vs.ctl, pData, by = "IID")
prs.ea <- merge(prs.ea, pData, by = "IID")
prs.cd <- merge(prs.cross.disorder, pData, by = "IID")
prs.gr.score <- merge(gr.score, pData, by = "IID")
#prs.bmi <- merge(prs.bmi, pData, by = "IID", all.x = T)


prs.test <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test4_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.test <- merge(prs.test, pData, by = "IID")
lm.test <- lm(prs.test$PRS ~ prs.test$Diagnosis)
nagelkerke(lm.test)
summary(lm.test)
# generate data for boxplot
prs.test.plot <- prs.test[, c("PRS", "Dx")]
cases <- prs.test.plot[prs.test.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.test.plot <- rbind(prs.test.plot, cases)
prs.test.plot$Dx <- factor(prs.test.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
# plot boxplots
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_SCZ2_eQTL_based.pdf", width = 12, height = 10)
ggplot(prs.test.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on schizophrenia GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()


# after smoothing
prs.smoothed <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed4_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.smoothed <- merge(prs.smoothed, pData, by = "IID")
lm.smoothed <- lm(prs.smoothed$PRS ~ prs.smoothed$Diagnosis)
nagelkerke(lm.smoothed)
summary(lm.smoothed)
# generate data for boxplot
prs.smoothed.plot <- prs.smoothed[, c("PRS", "Dx")]
cases <- prs.smoothed.plot[prs.smoothed.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.smoothed.plot <- rbind(prs.smoothed.plot, cases)
prs.smoothed.plot$Dx <- factor(prs.smoothed.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
# plot boxplots
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_SCZ2_smoothed.pdf", width = 12, height = 10)
ggplot(prs.smoothed.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on schizophrenia GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

#### Calculate statistics ###
# SCZ
lm.scz2 <- lm(prs.scz2$PRS ~ prs.scz2$Diagnosis)
nagelkerke(lm.scz2)
summary(lm.scz2)
# BPD
lm.bip <- lm(prs.bip$PRS ~ prs.bip$Diagnosis)
summary(lm.bip)
# MDD
lm.mdd2 <- lm(prs.mdd2$PRS ~ prs.mdd2$Diagnosis)
summary(lm.mdd2)
# CD
lm.cd <- lm(prs.cd$PRS ~ prs.cd$Diagnosis)
summary(lm.cd)
# ADHD
lm.adhd <- lm(prs.adhd$PRS ~ prs.adhd$Diagnosis)
summary(lm.adhd)
# ASD
lm.asd <- lm(prs.asd$PRS ~ prs.asd$Diagnosis)
summary(lm.asd)
# EA
lm.ea <- lm(prs.ea$PRS ~ prs.ea$Diagnosis)
summary(lm.ea)
# BIP vs SCZ
lm.bip.vs.scz <- lm(prs.bip.vs.scz$PRS ~ prs.bip.vs.scz$Diagnosis)
summary(lm.bip.vs.scz)
# BIP + SCZ vs CTL
lm.bip.scz.vs.ctl <- lm(prs.bip.scz.vs.ctl$PRS ~ prs.bip.scz.vs.ctl$Diagnosis)
summary(lm.bip.scz.vs.ctl)


### 3. Generate plots ###
# get six colors (for plots)
col <- viridis(20)

# change and relevel Dx column of SCZ2
prs.scz2.plot <- prs.scz2[, c("PRS", "Dx")]
cases <- prs.scz2.plot[prs.scz2.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.scz2.plot <- rbind(prs.scz2.plot, cases)
prs.scz2.plot$Dx <- factor(prs.scz2.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_SCZ2.pdf", width = 12, height = 10)
ggplot(prs.scz2.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on schizophrenia GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel Dx column of MDD2
prs.mdd2.plot <- prs.mdd2[, c("PRS", "Dx")]
cases <- prs.mdd2.plot[prs.mdd2.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.mdd2.plot <- rbind(prs.mdd2.plot, cases)
prs.mdd2.plot$Dx <- factor(prs.mdd2.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of MDD2 PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_MDD2.pdf", width = 12, height = 10)
ggplot(prs.mdd2.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on MDD GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of BIP 2018
prs.bip.plot <- prs.bip[, c("PRS", "Dx")]
cases <- prs.bip.plot[prs.bip.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.bip.plot <- rbind(prs.bip.plot, cases)
prs.bip.plot$Dx <- factor(prs.bip.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of BIP-2018 PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_BIP2018.pdf", width = 12, height = 10)
ggplot(prs.bip.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on bipolar disorder GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of ADHD
prs.adhd.plot <- prs.adhd[, c("PRS", "Dx")]
cases <- prs.adhd.plot[prs.adhd.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.adhd.plot <- rbind(prs.adhd.plot, cases)
prs.adhd.plot$Dx <- factor(prs.adhd.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of ADHD PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_ADHD.pdf", width = 12, height = 10)
ggplot(prs.adhd.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on ADHD GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)")
dev.off()

# Change and relevel the Dx column of ASD
prs.asd.plot <- prs.asd[, c("PRS", "Dx")]
cases <- prs.asd.plot[prs.asd.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.asd.plot <- rbind(prs.asd.plot, cases)
prs.asd.plot$Dx <- factor(prs.asd.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of ASD PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_ASD.pdf", width = 12, height = 10)
ggplot(prs.asd.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on ASD GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of BIP vs SCZ
prs.bip.vs.scz.plot <- prs.bip.vs.scz[, c("PRS", "Dx")]
cases <- prs.bip.vs.scz.plot[prs.bip.vs.scz.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.bip.vs.scz.plot <- rbind(prs.bip.vs.scz.plot, cases)
prs.bip.vs.scz.plot$Dx <- factor(prs.bip.vs.scz.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of BIP vs SCZ PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_BIPvsSCZ.pdf", width = 12, height = 10)
ggplot(prs.bip.vs.scz.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on BIP vs. SCZ GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of BIP + SCZ vs. CTL
prs.bip.scz.vs.ctl.plot <- prs.bip.scz.vs.ctl[, c("PRS", "Dx")]
cases <- prs.bip.scz.vs.ctl.plot[prs.bip.scz.vs.ctl.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.bip.scz.vs.ctl.plot <- rbind(prs.bip.scz.vs.ctl.plot, cases)
prs.bip.scz.vs.ctl.plot$Dx <- factor(prs.bip.scz.vs.ctl.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of BIP + SCZ vs. CTL PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_BIP+SCZvsCTL.pdf", width = 12, height = 10)
ggplot(prs.bip.scz.vs.ctl.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on BIP+SCZ vs. CTL GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of GR Score
prs.gr.score.plot <- prs.gr.score[, c("gr_score", "Dx")]
cases <- prs.gr.score.plot[prs.gr.score.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.gr.score.plot <- rbind(prs.gr.score.plot, cases)
prs.gr.score.plot$Dx <- factor(prs.gr.score.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of GR Score
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_GRScore.pdf", width = 12, height = 10)
ggplot(prs.gr.score.plot, aes(x = Dx, y = gr_score, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on GRScore summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of Educational Attainment (EA)
prs.ea.plot <- prs.ea[, c("PRS", "Dx")]
cases <- prs.ea.plot[prs.ea.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.ea.plot <- rbind(prs.ea.plot, cases)
prs.ea.plot$Dx <- factor(prs.ea.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of Educational Attainment (EA)
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_EA.pdf", width = 12, height = 10)
ggplot(prs.ea.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on educational attainment GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Dx column of Cross Disorder (CD) PRS
prs.cd.plot <- prs.cd[, c("PRS", "Dx")]
cases <- prs.cd.plot[prs.cd.plot$Dx != "CTL", ]
cases$Dx <- as.character(cases$Dx)
cases$Dx <- "CASE"
prs.cd.plot <- rbind(prs.cd.plot, cases)
prs.cd.plot$Dx <- factor(prs.cd.plot$Dx, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of Cross Disorder (CD) PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/boxplot_PRS_CD.pdf", width = 12, height = 10)
ggplot(prs.cd.plot, aes(x = Dx, y = PRS, fill = Dx)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on cross-disorder GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()


# Correlation Plot
# Save only important columns of pData and add PRS
pData <- pData[, c("IID", "Age", "Sex", "Dx", "pH", "PMI", "Suicide", "RIN", "Type_of_Death", "C1", "C2", "C3", "C4")]
pData$RIN2 <- pData$RIN^2
# Add PRS for different summary statistic files to pData
prs.scz2 <- data.frame(IID = prs.scz2$IID, PRS_SCZ2 = prs.scz2$PRS)
pData <- merge(pData, prs.scz2, by = "IID", all.x = T)
prs.mdd2 <- data.frame(IID = prs.mdd2$IID, PRS_MDD2 = prs.mdd2$PRS)
pData <- merge(pData, prs.mdd2, by = "IID", all.x = T)
prs.bip <- data.frame(IID = prs.bip$IID, PRS_BIP = prs.bip$PRS)
pData <- merge(pData, prs.bip, by = "IID", all.x = T)
prs.adhd <- data.frame(IID = prs.adhd$IID, PRS_ADHD = prs.adhd$PRS)
pData <- merge(pData, prs.adhd, by = "IID", all.x = T)
prs.asd <- data.frame(IID = prs.asd$IID, PRS_ASD = prs.asd$PRS)
pData <- merge(pData, prs.asd, by = "IID", all.x = T)
prs.bip.vs.scz <- data.frame(IID = prs.bip.vs.scz$IID, PRS_BIPvsSCZ = prs.bip.vs.scz$PRS)
pData <- merge(pData, prs.bip.vs.scz, by = "IID", all.x = T)
prs.bip.scz.vs.ctl <- data.frame(IID = prs.bip.scz.vs.ctl$IID, PRS_BIPSCZvsCTL = prs.bip.scz.vs.ctl$PRS)
pData <- merge(pData, prs.bip.scz.vs.ctl, by = "IID", all.x = T)
prs.gr.score <- data.frame(IID = prs.gr.score$IID, PRS_GRScore = prs.gr.score$gr_score)
pData <- merge(pData, prs.gr.score, by = "IID", all.x = T)
prs.ea <- data.frame(IID = prs.ea$IID, PRS_EA = prs.ea$PRS)
pData <- merge(pData, prs.ea, by = "IID", all.x = T)
prs.cd <- data.frame(IID = prs.cd$IID, PRS_CD = prs.cd$PRS)
pData <- merge(pData, prs.cd, by = "IID", all.x = T)
# 4. Remove all rows including an NA in any column and remove the IID column (NOT needed for the corrplot!)
pData.new <- pData[complete.cases(pData),]
pData.new$IID <- NULL
# 5. change BPD and MDD diagnosis in Dx to SCZ
pData.new$Dx <- as.character(pData.new$Dx)  # change factor column to character
pData.new[pData.new$Dx == "BPD" | pData.new$Dx == "MDD", "Dx"] <- "SCZ"  # change BPD and MDD to SCZ
pData.new[, c("Sex", "Dx", "Suicide")] <- lapply(pData.new[, c("Sex", "Dx", "Suicide")], factor)
#sapply(pData.new, class) # check type of columns
# relevel disease (Dx) to control (CTL)
pData.new$Dx <- factor(pData.new$Dx, levels = c("CTL", "SCZ"))
# save categorical and continuous variable names
classes <- sapply(pData.new, class) # get type of columns
categorical.cov = names(classes[classes == "factor"])
continuous.cov = names(classes[classes == "numeric"])
# generate contrast for design matrix
contrast <- lapply(categorical.cov, function(x) {fac = pData.new[, x]; return(contrasts(fac))})
names(contrast) <- categorical.cov
# generate design matrix
design <- model.matrix(~ PRS_SCZ2 + PRS_MDD2 + PRS_BIP + PRS_ADHD + PRS_ASD + PRS_BIPSCZvsCTL + PRS_CD + Age + Sex + Dx + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData.new, contrasts.arg = contrast)
design.cor <- corr.test(design[, -1], use = 'pairwise.complete.obs', method = "pearson", adjust = "none")
# change colnames
colnames(design.cor$r) <- c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Attention deficite hyperactivity disorder", "Autism", "Bipolar disorder + schizophrenia", "Cross-disorder", "Age", "Sex", "Diagnosis", "Suicide", "Non-violent type of death", "Violent type of death", "C1", "C2", "C3", "C4")
rownames(design.cor$r) <- c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Attention deficite hyperactivity disorder", "Autism", "Bipolar disorder + schizophrenia", "Cross-disorder", "Age", "Sex", "Diagnosis", "Suicide", "Non-violent type of death", "Violent type of death", "C1", "C2", "C3", "C4")
png("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/corrplot_covariates_PRS.png", type = "cairo", res = 300)
corrplot(as.matrix(design.cor$r), method = "circle", col = viridis(200, direction = -1), tl.col = "black")
dev.off()


# heatmap for all PRS GWAS
# save only IID and PRS columns
pData.new <- pData[, c("IID", "PRS_BIP", "PRS_SCZ2", "PRS_MDD2")]
# Remove all rows including an NA in any column and remove the IID column 
pData.new <- pData.new[complete.cases(pData.new),]
# set IIDs as rownames and remove IID column afterwards
row.names(pData.new) <- pData.new$IID
pData.new$IID <- NULL
# generate heatmap
samples <- as.data.frame(pData[, c("IID", "Dx", "Sex", "Age")])
samples <- samples[complete.cases(samples),]
samples$Dx <- as.character(samples$Dx)
samples[samples$Dx != "CTL", "Dx"] <- "Cases"
samples[samples$Dx == "CTL", "Dx"] <- "Controls"
row.names(samples) <- samples$IID
samples$IID <- NULL
colnames(samples) <- c("Diagnosis", "Sex", "Age")
ann_colors = list(Diagnosis = c(Cases = viridis(2)[2], Controls = viridis(2)[1]), Age = brewer.pal(n = 9, name = "Blues"), Sex = c(F = "darkorchid1", M = "darkblue"))
#pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/heatmap_PRS_all_GWAS.pdf", width = 12, height = 9)
pheatmap(t(pData.new), scale = "row", annotation_col = samples, annotation_colors = ann_colors, font_size = 0.6, show_colnames = F, 
         cex = 1.5, cutree_cols = 5, labels_row = c("PGC bipolar disorder GWAS (2018)", "PGC schizophrenia GWAS (2014)", "PGC major depression disorder GWAS (2018)"))#, color = rev(brewer.pal(n = 11, name = "RdBu")))
#dev.off()


# heatmap only SCZ, BIP, MDD and CrossDisorder PRS
pData.new <- pData[, c("IID", "PRS_SCZ2", "PRS_MDD2", "PRS_BIP", "PRS_CD")]
pData.new <- pData.new[complete.cases(pData.new),]
# set IIDs as rownames and remove IID column afterwards
row.names(pData.new) <- pData.new$IID
pData.new$IID <- NULL
# generate heatmap
samples <- as.data.frame(pData[, c("IID", "Dx", "Sex", "Age")])
samples <- samples[complete.cases(samples),]
samples$Dx <- as.character(samples$Dx)
samples[samples$Dx != "CTL", "Dx"] <- "Cases"
samples[samples$Dx == "CTL", "Dx"] <- "Controls"
row.names(samples) <- samples$IID
samples$IID <- NULL
colnames(samples) <- c("Diagnosis", "Sex", "Age")
ann_colors = list(Diagnosis = c(Cases = viridis(2)[2], Controls = viridis(2)[1]), Age = brewer.pal(n = 9, name = "Blues"), Sex = c(F = "darkmagenta", M = "darkblue"))
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/figures/heatmap_PRS_only_SCZ_BPD_MDD_CD_GWAS.pdf", width = 14, height = 8)
pheatmap(t(pData.new), scale = "row", annotation_col = samples, annotation_colors = ann_colors, font_size = 0.6, show_colnames = F,
         cex = 1.5, labels_row = c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Cross-disorder"), color = viridis(200))
dev.off()


# session info
sessionInfo()

