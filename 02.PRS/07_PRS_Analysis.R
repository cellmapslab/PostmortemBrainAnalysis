####################
### PRS Analysis ###
####################
# In this R script, we will analyze and plot all generated polygenic risk scores

# Load packages (for own R installation - version 3.6)
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")
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

# function to calculate Nagelkerkes R2 and p-value 
# and generate a boxplot figure for a given PRS file (best PRS of PRSice-2)
EvaluatePRS <- function(input, gwas, seperator, output){
  # 1. read file
  prs <- fread(file = input, sep = seperator, header = T, stringsAsFactors = F, select = c("IID", "PRS"))
  # 2. add pData information to PRS
  prs <- merge(prs, pData, by = "IID")
  # change Dx to numeric
  prs[prs$Dx == "CTL", "Dx"] <- "0"
  prs[prs$Dx == "CASE", "Dx"] <- "1"
  prs$Dx <- as.numeric(as.character(prs$Dx))
  # 3. Get Nagelkerkes R2 and p-value
  lm.prs <- lm(prs$Dx ~ prs$PRS)
  r2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
  r2.2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
  sum.prs <- summary(lm.prs)
  p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
  print(paste0(gwas, ":\nNagelkerke R2 = ", r2, "\nCox and Snell R2 = ", r2.2, "\nP-value = ", p.value))
  # 4. Generate boxplot
  col <- viridis(20) # get six colors (for plots)
  # change and relevel Diagnosis column of SCZ2
  prs.plot <- prs[, c("PRS", "Diagnosis")]
  cases <- prs.plot[prs.plot$Diagnosis != "CTL", ]
  cases$Diagnosis <- as.character(cases$Diagnosis)
  cases$Diagnosis <- "CASE"
  prs.plot <- rbind(prs.plot, cases)
  prs.plot$Diagnosis <- factor(prs.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
  # boxplot of PRS
  ggplot(prs.plot, aes(x = Diagnosis, y = PRS, fill = Diagnosis)) +
    scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
    theme_classic() +
    theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
    geom_boxplot(alpha = 0.7) + 
    scale_fill_manual(values = viridis(5), guide = F) +
    ggtitle(paste0("PRS based on ", gwas, " GWAS summary statistic")) + 
    ylab("Polygenic risk score (PRS)") 
  ggsave(output, width = 12, height = 10, units = "in")
  return(prs)
}

# get PRS output path and generate the figures folder
path = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS"
# generate output (05_figures) folder
system(paste0("mkdir ", path, "/05_figures"))

# load, filter and save phenotype data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
index <- which(is.na(pData(gene_summaries)$C1))  # remove all samples with C1 = NA
dim(gene_summaries)
pData <- pData(gene_summaries)
rm(gene_summaries)

### PRSice-2 ###
adhd_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/ADHD-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ADHD-2017_PRSice-2_PRS.best"),
                           gwas = "ADHD-2017", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_ADHD-2017.pdf"))

asd_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/ASD-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ASD-2019_PRSice-2_PRS.best"),
                          gwas = "ASD-2017", seperator = " ",
                          output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_ASD-2017.pdf"))

bpd_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/BPD-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BPD-2018_PRSice-2_PRS.best"),
                          gwas = "BPD-2018", seperator = " ",
                          output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_BPD-2018.pdf"))

bpd_vs_scz_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/BPDvsSCZ-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BPDvsSCZ-2018_PRSice-2_PRS.best"),
                                 gwas = "BPDvsSCZ-2018", seperator = " ",
                                 output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_BPDvsSCZ-2018.pdf"))

bpd_scz_vs_ctl_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/BPDSCZvsCONT-2018/"),
                                     gwas = "BPDSCZvsCONT-2018", seperator = " ",
                                     output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_BPDSCZvsCONT-2018.pdf"))

bmi_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/BMI-2015/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BMI-2015_PRSice-2_PRS.best"),
                          gwas = "BMI-2015", seperator = " ",
                          output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_BMI-2015.pdf"))

cdg1_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/CDG1-2013/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_CDG1-2013_PRSice-2_PRS.best"),
                           gwas = "CDG1-2013", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_CDG1-2013.pdf"))

cdg2_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/CDG2-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_CDG2-2018_PRSice-2_PRS.best"),
                           gwas = "CDG2-2018", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_CDG2-2018.pdf"))

ea_prsice <- EvaluatePRS(input = paste0(path, "/PRSice-2/EA-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_EA-2018_PRSice-2_PRS.best"),
                         gwas = "EA-2018", seperator = " ",
                         output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_EA-2018.pdf"))

intelligence_prsice <- EvaluatePRS(input = paste0(path, "/PRSice-2/Intelligence-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_Intelligence-2018_PRSice-2_PRS.best"),
                                   gwas = "Intelligence-2018", seperator = " ",
                                   output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_Intelligence-2018.pdf"))

mdd1_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/MDD1-2012/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD1-2012_PRSice-2_PRS.best"), 
                           gwas = "MDD1-2012", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_MDD1-2012.pdf"))

mdd2_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/MDD2-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD2-2018_PRSice-2_PRS.best"), 
                           gwas = "MDD2-2018", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_MDD2-2018.pdf"))

mdd3_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/MDD3-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD3-2019_PRSice-2_PRS.best"), 
                           gwas = "MDD3-2019", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_MDD3-2019.pdf"))

scz2_prsice <- EvaluatePRS(input = paste0(path, "/04_PRSice-2/SCZ2-2014/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_SCZ2-2014_PRSice-2_PRS.best"),
                           gwas = "SCZ2-2014", seperator = " ",
                           output = paste0(path, "/05_figures/boxplot_PRSice-2_PRS_SCZ2-2014.pdf"))

### Load ALL Scores and get R2 + p-value for threshold of 1, 0.01 and 0.00000005 (5e-08)
scz2.scores <- fread(input = paste0(path, "/04_PRSice-2/SCZ2-2014/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_SCZ2-2014_PRSice-2_PRS.all.score"))
names <- colnames(scz2.scores)[-1]
scz2.scores$`1` <- NULL
colnames(scz2.scores) <- names
scz2.scores <- merge(scz2.scores, pData, by = "IID")
scz2.scores[scz2.scores$Dx == "CTL", "Dx"] <- "0"
scz2.scores[scz2.scores$Dx == "CASE", "Dx"] <- "1"
scz2.scores$Dx <- as.numeric(as.character(scz2.scores$Dx))
# 3. Get Nagelkerkes R2 and p-value
prs.scz2 <- lm(scz2.scores$Dx ~ scz2.scores$`5e-08`)
r2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
r2.2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
sum.prs <- summary(prs.scz2)
p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
print(paste0("5e-08: Nagelkerke R2 = ", r2, "\nCos and Snell R2 = ", r2.2, "\nP-value = ", p.value))
prs.scz2 <- lm(scz2.scores$Dx ~ scz2.scores$`0.01` + scz2.scores$Age + scz2.scores$Sex + scz2.scores$PMI + scz2.scores$Suicide + scz2.scores$Type_of_Death)
prs.scz2 <- lm(scz2.scores$Dx ~ scz2.scores$`0.0222001`)
r2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
r2.2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
sum.prs <- summary(prs.scz2)
p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
print(paste0("0.01: Nagelkerke R2 = ", r2, "\nCos and Snell R2 = ", r2.2, "\nP-value = ", p.value))
prs.scz2 <- lm(scz2.scores$Dx ~ scz2.scores$`1`)
               r2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
r2.2 <- nagelkerke(prs.scz2)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
sum.prs <- summary(prs.scz2)
p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
print(paste0("1: Nagelkerke R2 = ", r2, "\nCos and Snell R2 = ", r2.2, "\nP-value = ", p.value))

## Load PRSice-2 SNP files and get number of SNPs
scz_snps <- fread(paste0(path, "/PRSice-2/SCZ2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_scz2_filtered_PRSice-2_PRS.snp"))
length(scz_snps[as.numeric(as.character(scz_snps$P)) < 0.0222001, ]$P)
length(scz_snps[as.numeric(as.character(scz_snps$P)) < 5e-08, ]$P)
length(scz_snps[as.numeric(as.character(scz_snps$P)) < 0.05, ]$P)
length(scz_snps[as.numeric(as.character(scz_snps$P)) < 1, ]$P)

mdd_snps <- fread(paste0(path, "/PRSice-2/MDD2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_mdd2_filtered_PRSice-2_PRS.snp"))
length(mdd_snps[as.numeric(as.character(mdd_snps$P)) < 0.0354001, ]$P)
length(mdd_snps[as.numeric(as.character(mdd_snps$P)) < 5e-08, ]$P)
length(mdd_snps[as.numeric(as.character(mdd_snps$P)) < 0.05, ]$P)
length(mdd_snps[as.numeric(as.character(mdd_snps$P)) < 1, ]$P)

bpd_snps <- fread(paste0(path, "/PRSice-2/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_filtered_PRSice-2_PRS.snp"))
length(bpd_snps[as.numeric(as.character(bpd_snps$P)) < 0.0154501, ]$P)
length(bpd_snps[as.numeric(as.character(bpd_snps$P)) < 5e-08, ]$P)
length(bpd_snps[as.numeric(as.character(bpd_snps$P)) < 0.05, ]$P)
length(bpd_snps[as.numeric(as.character(bpd_snps$P)) < 1, ]$P)

cd_snps <- fread(paste0(path, "/PRSice-2/Cross_Disorder/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_cd_filtered_PRSice-2_PRS.snp"))
length(cd_snps[as.numeric(as.character(cd_snps$P)) < 0.00005005, ]$P)
length(cd_snps[as.numeric(as.character(cd_snps$P)) < 5e-08, ]$P)
length(cd_snps[as.numeric(as.character(cd_snps$P)) < 0.05, ]$P)
length(cd_snps[as.numeric(as.character(cd_snps$P)) < 1, ]$P)


# load PRS LDpred-inf scores
# Remember to add a comma "," at the end of the header line (or remove all commas from the rest of the lines) otherwise read.table or fread can not read in the data correctly. 
# load SCZ2 PRS
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load MDD2 PRS
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
prs.cross.disorder <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/Cross_Disorder/03_validate/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_cd_filtered_PRS_LDpred_p1.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load GR Score and filter for our datasets
gr.score <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/03.TRPS/03_TRPS.txt", sep = "\t", header = T, stringsAsFactors = F, select = c("IID", "gr_score"))
gr.score <- gr.score[!startsWith(gr.score$IID, "MD") & !endsWith(gr.score$IID, "ZNP")]
gr.score <- gr.score[gr.score$IID %in% pData$IID]
gr.score <- gr.score[, IID:=as.integer(IID)]

# Add diagnosis column to the PRS
prs.scz2 <- merge(prs.scz2, pData, by = "IID")
prs.adhd <- merge(prs.adhd, pData, by = "IID")
prs.asd <- merge(prs.asd, pData, by = "IID")
prs.bip.vs.scz <- merge(prs.bip.vs.scz, pData, by = "IID")
prs.bip.scz.vs.ctl <- merge(prs.bip.scz.vs.ctl, pData, by = "IID")
prs.ea <- merge(prs.ea, pData, by = "IID")
prs.gr.score <- merge(gr.score, pData, by = "IID")
#prs.bmi <- merge(prs.bmi, pData, by = "IID", all.x = T)

# after smoothing
prs.smoothed <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed4_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.smoothed <- merge(prs.smoothed, pData, by = "IID")
lm.smoothed <- lm(prs.smoothed$PRS ~ prs.smoothed$Diagnosis)
nagelkerke(lm.smoothed)
summary(lm.smoothed)
# generate data for boxplot
prs.smoothed.plot <- prs.smoothed[, c("PRS", "Diagnosis")]
cases <- prs.smoothed.plot[prs.smoothed.plot$Diagnosis != "CTL", ]
cases$Diagnosis <- as.character(cases$Diagnosis)
cases$Diagnosis <- "CASE"
prs.smoothed.plot <- rbind(prs.smoothed.plot, cases)
prs.smoothed.plot$Diagnosis <- factor(prs.smoothed.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
# plot boxplots
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/05_figures/boxplot_PRS_SCZ2_smoothed.pdf", width = 12, height = 10)
ggplot(prs.smoothed.plot, aes(x = Diagnosis, y = PRS, fill = Diagnosis)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on schizophrenia GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()


### 3. Generate plots ###
# get six colors (for plots)
col <- viridis(20)

# change and relevel Diagnosis column of SCZ2
prs.scz2.plot <- prs.scz2[, c("PRS", "Diagnosis")]
cases <- prs.scz2.plot[prs.scz2.plot$Diagnosis != "CTL", ]
cases$Diagnosis <- as.character(cases$Diagnosis)
cases$Diagnosis <- "CASE"
prs.scz2.plot <- rbind(prs.scz2.plot, cases)
prs.scz2.plot$Diagnosis <- factor(prs.scz2.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_figures/boxplot_LDpred_PRS_SCZ2.pdf", width = 12, height = 10)
ggplot(prs.scz2.plot, aes(x = Diagnosis, y = PRS, fill = Diagnosis)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on schizophrenia GWAS summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()

# Change and relevel the Diagnosis column of GR Score
prs.gr.score.plot <- prs.gr.score[, c("gr_score", "Diagnosis")]
cases <- prs.gr.score.plot[prs.gr.score.plot$Diagnosis != "CTL", ]
cases$Diagnosis <- as.character(cases$Diagnosis)
cases$Diagnosis <- "CASE"
prs.gr.score.plot <- rbind(prs.gr.score.plot, cases)
prs.gr.score.plot$Diagnosis <- factor(prs.gr.score.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))

# Boxplots of GR Score
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_figures/boxplot_PRS_GRScore.pdf", width = 12, height = 10)
ggplot(prs.gr.score.plot, aes(x = Diagnosis, y = gr_score, fill = Diagnosis)) +
  scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
  theme_classic() +
  theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
  geom_boxplot(alpha = 0.7) + 
  scale_fill_manual(values = viridis(5), guide = F) +
  ggtitle("PRS based on GRScore summary statistic") + 
  ylab("Polygenic risk score (PRS)") 
dev.off()


# Correlation Plot
# Add PRS for different summary statistic files to pData
prs_scz2 <- data.frame(IID = scz2_prsice$IID, PRS_SCZ = scz2_prsice$PRS)
pData <- merge(pData, prs_scz2, by = "IID", all.x = T)
prs_mdd2 <- data.frame(IID = mdd2_prsice$IID, PRS_MDD = mdd2_prsice$PRS)
pData <- merge(pData, prs_mdd2, by = "IID", all.x = T)
prs_bip <- data.frame(IID = bip_prsice$IID, PRS_BIP = bip_prsice$PRS)
pData <- merge(pData, prs_bip, by = "IID", all.x = T)
prs_adhd <- data.frame(IID = adhd_prsice$IID, PRS_ADHD = adhd_prsice$PRS)
pData <- merge(pData, prs_adhd, by = "IID", all.x = T)
prs_asd <- data.frame(IID = asd_prsice$IID, PRS_ASD = asd_prsice$PRS)
pData <- merge(pData, prs_asd, by = "IID", all.x = T)
prs_bip_vs_scz <- data.frame(IID = bip_vs_scz_prsice$IID, PRS_BIPvsSCZ = bip_vs_scz_prsice$PRS)
pData <- merge(pData, prs_bip_vs_scz, by = "IID", all.x = T)
prs_bip_scz_vs_ctl <- data.frame(IID = bip_scz_vs_ctl_prsice$IID, PRS_BIPSCZvsCTL = bip_scz_vs_ctl_prsice$PRS)
pData <- merge(pData, prs_bip_scz_vs_ctl, by = "IID", all.x = T)
prs_gr_score <- data.frame(IID = gr.score$IID, PRS_GRScore = gr.score$gr_score)
pData <- merge(pData, prs_gr_score, by = "IID", all.x = T)
prs_ea <- data.frame(IID = ea_prsice$IID, PRS_EA = ea_prsice$PRS)
pData <- merge(pData, prs_ea, by = "IID", all.x = T)
prs_cd <- data.frame(IID = cd_prsice$IID, PRS_CD = cd_prsice$PRS)
pData <- merge(pData, prs_cd, by = "IID", all.x = T)
prs_bmi <- data.frame(IID = bmi_prsice$IID, PRS_BMI = bmi_prsice$PRS)
pData <- merge(pData, prs_bmi, by = "IID", all.x = T)
# 4. Remove all rows including an NA in any column and remove the IID column (NOT needed for the corrplot!)
pData.new <- pData[, c("Age", "Sex", "Dx", "Suicide", "Type_of_Death", "C1", "C2", "C3", "C4", "PRS_SCZ", "PRS_MDD", "PRS_BIP", "PRS_CD", "PRS_BIPSCZvsCTL", "PRS_ADHD", "PRS_ASD", "PRS_EA", "PRS_BMI")]
# 5. change BPD and MDD diagnosis in Diagnosis to SCZ
# save categorical and continuous variable names
classes <- sapply(pData.new, class) # get type of columns
categorical.cov = names(classes[classes == "factor"])
continuous.cov = names(classes[classes == "numeric"])
# generate contrast for design matrix
contrast <- lapply(categorical.cov, function(x) {fac = pData.new[, x]; return(contrasts(fac))})
names(contrast) <- categorical.cov
# generate design matrix
design <- model.matrix(~ PRS_SCZ + PRS_MDD + PRS_BIP + PRS_CD + PRS_ADHD + PRS_ASD + PRS_EA + PRS_BMI + Age + Sex + Dx + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData.new, contrasts.arg = contrast)
design.cor <- corr.test(design[, -1], use = 'pairwise.complete.obs', method = "pearson", adjust = "none")
# change colnames
colnames(design.cor$r) <- c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Cross-disorder", "Attention deficite hyperactivity disorder", "Autism", "Educational attainment", "Body mass index", "Age", "Sex", "Diagnosis", "Suicide", "Non-violent type of death", "Violent type of death", "C1", "C2", "C3", "C4")
rownames(design.cor$r) <- c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Cross-disorder", "Attention deficite hyperactivity disorder", "Autism", "Educational attainment", "Body mass index", "Age", "Sex", "Diagnosis", "Suicide", "Non-violent type of death", "Violent type of death", "C1", "C2", "C3", "C4")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_figures/corrplot_covariates_PRS.pdf")
corrplot(as.matrix(design.cor$r), method = "circle", col = viridis(200, direction = -1), tl.col = "black")
dev.off()

# heatmap for all PRS GWAS
# save only IID and PRS columns
#pData.new <- pData[, c("IID", "PRS_BIP", "PRS_SCZ", "PRS_MDD", "PRS_CD", "PRS_ADHD", "PRS_ASD", "PRS_EA", "PRS_BMI")]
pData.new <- pData[, c("IID", "PRS_MDD", "PRS_BIP", "PRS_SCZ", "PRS_CD")]
# set IIDs as rownames and remove IID column afterwards
row.names(pData.new) <- pData.new$IID
# generate heatmap
samples <- as.data.frame(pData[, c("IID", "Diagnosis", "Sex", "Age")])
samples <- samples[complete.cases(samples),]
samples$Diagnosis <- as.character(samples$Diagnosis)
samples[samples$Diagnosis != "CTL", "Diagnosis"] <- "Cases"
samples[samples$Diagnosis == "CTL", "Diagnosis"] <- "Controls"
row.names(samples) <- samples$IID
samples$IID <- NULL
colnames(samples) <- c("Diagnosis", "Sex", "Age")
ann_colors = list(Diagnosis = c(Cases = viridis(8)[2], Controls = "gray"), Age = brewer.pal(n = 9, name = "Blues"), Sex = c(F = "violetred", M = "darkblue"))
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_figures/heatmap_PRS_SCZ_BIP_MDD_CD_GWAS.pdf", width = 12, height = 9)
pheatmap(t(pData.new[, 2:length(pData.new)]), scale = "row", cluster_rows = F, annotation_col = samples, annotation_colors = ann_colors, font_size = 0.6, show_colnames = F, 
         cex = 1.5, labels_row = c("Major depression disorder", "Bipolar disorder", "Schizophrenia", "Cross-Disorder"))#, color = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
p <- pheatmap(t(pData.new[, 2:length(pData.new)]), scale = "row", cluster_rows = F, annotation_col = samples, annotation_colors = ann_colors, font_size = 0.6, show_colnames = F, 
              cex = 1.5, labels_row = c("Major depression disorder", "Bipolar disorder", "Schizophrenia", "Cross-Disorder")) #, color = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
iid.order <- p$tree_col$order
row_names <- as.character(pData.new[iid.order, "IID"])
pData.new2 <- pData[, c("IID", "PRS_ADHD", "PRS_ASD", "PRS_EA", "PRS_BMI")]
row.names(pData.new2) <- as.character(pData.new2$IID)
pData.new2$IID <- NULL
pData.plot <- t(pData.new2[row_names, ])
# generate heatmap
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_figures/heatmap_PRS_ADHD_ASD_EA_BMI_GWAS.pdf", width = 12, height = 9)
pheatmap(pData.plot, scale = "row", cluster_rows = F, cluster_cols = F, font_size = 0.6, show_colnames = F, 
              cex = 1.5, labels_row = c("Attention deficite hyperactivity disorder", "Autism", "Educational attainment", "Body mass index"))#, color = rev(brewer.pal(n = 11, name = "RdBu")))
dev.off()
pheatmap(t(pData.new2)[, row_names], scale = "row", cluster_rows = F, cluster_cols = F, font_size = 0.6, show_colnames = F, 
         cex = 1.5, labels_row = c("Attention deficite hyperactivity disorder", "Autism", "Educational attainment", "Body mass index"))#, color = rev(brewer.pal(n = 11, name = "RdBu")))

# heatmap only SCZ, BIP, MDD and CrossDisorder PRS
pData.new <- pData[, c("IID", "PRS_SCZ2", "PRS_MDD2", "PRS_BIP", "PRS_CD")]
pData.new <- pData.new[complete.cases(pData.new),]
# set IIDs as rownames and remove IID column afterwards
row.names(pData.new) <- pData.new$IID
pData.new$IID <- NULL
# generate heatmap
samples <- as.data.frame(pData[, c("IID", "Diagnosis", "Sex", "Age")])
samples <- samples[complete.cases(samples),]
samples$Diagnosis <- as.character(samples$Diagnosis)
samples[samples$Diagnosis != "CTL", "Diagnosis"] <- "Cases"
samples[samples$Diagnosis == "CTL", "Diagnosis"] <- "Controls"
row.names(samples) <- samples$IID
samples$IID <- NULL
colnames(samples) <- c("Diagnosis", "Sex", "Age")
ann_colors = list(Diagnosis = c(Cases = viridis(2)[2], Controls = viridis(2)[1]), Age = brewer.pal(n = 9, name = "Blues"), Sex = c(F = "darkmagenta", M = "darkblue"))
pdf("/storage/scratch/users/karolina.worf/MPIP/02.PRS/05_figures/heatmap_PRS_only_SCZ_BPD_MDD_CD_GWAS.pdf", width = 14, height = 8)
pheatmap(t(pData.new), scale = "row", annotation_col = samples, annotation_colors = ann_colors, font_size = 0.6, show_colnames = F,
         cex = 1.5, labels_row = c("Schizophrenia", "Major depression disorder", "Bipolar disorder", "Cross-disorder"), color = viridis(200))
dev.off()


# session info
sessionInfo()

