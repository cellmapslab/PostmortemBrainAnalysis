####################
### PRS Analysis ###
####################
# In this R script, we will analyze and plot all generated polygenic risk scores

# 1. Load packages (for own R installation - version 3.6)
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

# 2. Load data (all different generated PRS)
# load, filter and save phenotype data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
dim(geneSummaries)
pData <- pData(geneSummaries)

# load smoothed PRS scores
# Remember to add a comma "," at the end of the header line (or remove all commas from the rest of the lines) otherwise read.table or fread can not read in the data correctly. 
# load smoothed SCZ2 PRS (normal and smoothed weights)
prs.scz2.1 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed1_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load smoothed SCZ2 PRS (zeros and smoothed weights)
prs.scz2.2 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed2_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load LD-smoothed SCZ2 PRS (normal and smoothed weights)
prs.scz2.3 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed3_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# load LD-smoothed SCZ2 PRS (zeros and smoothed weights)
prs.scz2.4 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed4_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))

prs.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.gene.cases <-  fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/gene_cis-eqtl_cases_only_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/transcript_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/transcript_cis-eqtl_cases_only_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/exon_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/exon_cis-eqtl_cases_only_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.gene.negative <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/negative_gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))

prs.bios.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/bios_blood_gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.bios.exon <-  fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/bios_blood_exon_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.cmc.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/cmc_gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.cmc.iso <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/cmc_isoform_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.gtex.wb <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/gtex_wb_gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.gtex.brain <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/results/gtex_brains_gene_cis-eqtl_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))


# Add diagnosis column to the PRS
prs.scz2.1 <- merge(prs.scz2.1, pData, by = "IID")
prs.scz2.2 <- merge(prs.scz2.2, pData, by = "IID")
prs.scz2.3 <- merge(prs.scz2.3, pData, by = "IID")
prs.scz2.4 <- merge(prs.scz2.4, pData, by = "IID")

prs.gene <- merge(prs.gene, pData, by = "IID")
prs.gene.cases <- merge(prs.gene.cases, pData, by = "IID")
prs.transcript <- merge(prs.transcript, pData, by = "IID")
prs.transcript.cases <- merge(prs.transcript.cases, pData, by = "IID")
prs.exon <- merge(prs.exon, pData, by = "IID")
prs.exon.cases <- merge(prs.exon.cases, pData, by = "IID")
prs.gene.negative <- merge(prs.gene.negative, pData, by = "IID")

prs.bios.gene <- merge(prs.bios.gene, pData, by = "IID")
prs.bios.exon <- merge(prs.bios.exon, pData, by = "IID")
prs.cmc.gene <- merge(prs.cmc.gene, pData, by = "IID")
prs.cmc.iso <- merge(prs.cmc.iso, pData, by = "IID")
prs.gtex.wb <- merge(prs.gtex.wb, pData, by = "IID")
prs.gtex.brain <- merge(prs.gtex.brain, pData, by = "IID")

# get p-value & R²
lm.smoothed <- lm(prs.exon$PRS ~ prs.exon$Diagnosis)
nagelkerke(lm.smoothed) # for R²
summary(lm.smoothed) # for p-value

### 3. Boxplots of different PRS (with diagnosis, sex, age and suicide) ###
# function: transform data for the plot using the prs + pData dataframe as input
TransformData <- function(prs){ 
  # 1. only take neccessary columns
  prs <- prs[, c("PRS", "Dx", "Sex", "Age", "Suicide", "Type_of_Death")]
  # 2. get and transform the diagnosis data
  dx <- data.frame(Label = "Diagnosis", Variable = prs$Dx, PRS = prs$PRS)
  dx$Variable <- as.character(dx$Variable)
  dx[dx$Variable != "CTL", "Variable"] <- "Cases"
  dx[dx$Variable == "CTL", "Variable"] <- "Controls"
  # 3. get and transform the different diagnosis (MDD, BPD, SCZ) for the cases dataset
  cases <- data.frame(Label = "Cases", Variable = prs[prs$Dx != "CTL", "Dx"], PRS = prs[prs$Dx != "CTL", "PRS"])
  cases$Dx <- as.character(cases$Dx)
  cases[cases$Dx == "SCZ", "Dx"] <- "Schizophrenia"
  cases[cases$Dx == "MDD", "Dx"] <- "Major depression disorder"
  cases[cases$Dx == "BPD", "Dx"] <- "Bipolar disorder"
  colnames(cases) <- c("Label", "Variable", "PRS")
  # 4. get and transform the sex dataset
  sex <- data.frame(Label = "Sex", Variable = prs$Sex, PRS = prs$PRS)
  sex$Variable <- as.character(sex$Variable)
  sex[sex$Variable == "M", "Variable"] <- "Males"
  sex[sex$Variable == "F", "Variable"] <- "Females"
  # 5. get and transform the age dataset (transform into two groups: once <= 51.9 and once bigger or equal 51.9)
  age <- data.frame(Label = "Age", Variable = prs$Age, PRS = prs$PRS)
  age[age$Variable >= 51.9, "Variable"] <- ">= 51.9"
  age[age$Variable != ">= 51.9", "Variable"] <- "< 51.9"
  # 6. get and transform the suicide dataset
  suicide <- data.frame(Label = "Suicide", Variable = prs$Suicide, PRS = prs$PRS)
  suicide$Variable <- as.character(suicide$Variable)
  suicide[suicide$Variable  == "Y", "Variable"] <- "Yes"
  suicide[suicide$Variable  == "N", "Variable"] <- "No"
  # 7. get and transform the type of death dataset
  tod <- data.frame(Label = "Type of Death", Variable = prs$Type_of_Death, PRS = prs$PRS)
  # 8. merge all datasets for generating the plot
  prs.plot <- rbind(dx, cases, sex, age, suicide, tod)
  prs.plot$Variable <- factor(prs.plot$Variable, levels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia", "Males", "Females", "< 51.9", ">= 51.9", "No", "Yes", "Natural", "Non-violent", "Violent"))
  plot <- ggplot(prs.plot, aes(x = Variable, y = PRS, fill = Label)) +
    theme_classic() +
    theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
    scale_fill_manual(values = plasma(6)) +
    geom_boxplot(alpha = 0.7) +
    ylab("Polygenic risk score (PRS)") 
  return(plot)  # return the plot
}

# boxplot of unsmoothed/normal SCZ PRS
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2 <- merge(prs.scz2, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/normal_SCZ_PRS.pdf", width = 16, height = 9)
test <- TransformData(prs.scz2)
test + ggtitle("Normal PRS")
dev.off()

# boxplot of SCZ PRS with filtered SNPs nearby genes
prs.scz2.t1 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test3_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t1 <- merge(prs.scz2.t1, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/filtering_SNPs_nearby_genes_SCZ_PRS.pdf", width = 16, height = 9)
plot.t1 <- TransformData(prs.scz2.t1)
plot.t1 + ggtitle("Filtered PRS (SNPs nearby genes)")
dev.off()

# boxplot of SCZ PRS with filtered significant gene-level eQTL SNPs
prs.scz2.t2 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test1_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t2 <- merge(prs.scz2.t2, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/filtering_sig_gene_eQTL_SNPs_SCZ_PRS.pdf", width = 16, height = 9)
plot.t2 <- TransformData(prs.scz2.t2)
plot.t2 + ggtitle("Filtered PRS (significant gene eSNPs)")
dev.off()

# boxplot of SCZ PRS with reweighted SNPs connected through significant gene-level eQTLs
prs.scz2.t3 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed1_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t3 <- merge(prs.scz2.t3, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/reweighting_sig_gene_eQTL_SNPs_SCZ_PRS.pdf", width = 16, height = 9)
plot.t3 <- TransformData(prs.scz2.t3)
plot.t3 + ggtitle("Reweighted PRS (smoothing: significant gene eSNPs)")
dev.off()

# boxplot of SCZ PRS with reweighted SNPs connected through significant gene-level eQTLs + LD SNPs
prs.scz2.t4 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed3_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t4 <- merge(prs.scz2.t4, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/reweighting_sig_gene_eQTL_SNPs_with_LD_SCZ_PRS.pdf", width = 16, height = 9)
plot.t4 <- TransformData(prs.scz2.t4)
plot.t4 + ggtitle("Reweighted PRS (smoothing: sig. gene eSNPs + LD SNPs)")
dev.off()

# boxplot of SCZ PRS with filtering + reweighted SNPs connected through significant gene-level eQTLs
prs.scz2.t5 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed2_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t5 <- merge(prs.scz2.t5, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/filtering+reweighting_sig_gene_eQTL_SNPs_SCZ_PRS.pdf", width = 16, height = 9)
plot.t5 <- TransformData(prs.scz2.t5)
plot.t5 + ggtitle("Filtering + reweighted PRS (smoothing: significant gene eSNPs)")
dev.off()

# boxplot of SCZ PRS with filtering + reweighted SNPs connected through significant gene-level eQTLs + LD SNPs
prs.scz2.t6 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed4_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t6 <- merge(prs.scz2.t6, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/filtering+reweighting_sig_gene_eQTL_SNPs_with_LD_SCZ_PRS.pdf", width = 16, height = 9)
plot.t6 <- TransformData(prs.scz2.t6)
plot.t6 + ggtitle("Filtering + reweighted PRS (smoothing: sig. gene eSNPs + LD SNPs)")
dev.off()

# boxplot of SCZ PRS with randomized SNP values
prs.scz2.t7 <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test3_randomized_smoothed_PRS_LDpred_p3.0000e-01.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
prs.scz2.t7 <- merge(prs.scz2.t7, pData, by = "IID")  # Add diagnosis column to the PRS
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/randomized_SNP_values_SCZ_PRS.pdf", width = 16, height = 9)
plot.t7 <- TransformData(prs.scz2.t7)
plot.t7 + ggtitle("Randomized SNP values PRS")
dev.off()

#cases <- prs.scz2.1[prs.scz2.1$Dx != "CTL", "PRS"]
#ctrls <- prs.scz2.1[prs.scz2.1$Dx == "CTL", "PRS"]
#t.test(x = cases, y = ctrls, alternative = "greater")

#cases <- prs.scz2.2[prs.scz2.2$Dx != "CTL", "PRS"]
#ctrls <- prs.scz2.2[prs.scz2.2$Dx == "CTL", "PRS"]
#t.test(x = cases, y = ctrls, alternative = "greater")

#cases <- prs.scz2.3[prs.scz2.3$Dx != "CTL", "PRS"]
#ctrls <- prs.scz2.3[prs.scz2.3$Dx == "CTL", "PRS"]
#t.test(x = cases, y = ctrls, alternative = "greater")

#cases <- prs.scz2.4[prs.scz2.4$Dx != "CTL", "PRS"]
#ctrls <- prs.scz2.4[prs.scz2.4$Dx == "CTL", "PRS"]
#t.test(x = cases, y = ctrls, alternative = "greater")


# 3. Generate plots
# get six colors (for plots)
col <- viridis(20)

# change prs.scz2 from data.table to data.frame for releveling the Dx column
prs.scz2.1 <- as.data.frame(prs.scz2.1)
prs.scz2.1$Dx <- as.factor(prs.scz2.1$Dx)
prs.scz2.1$Dx <- factor(prs.scz2.1$Dx, levels = c("CTL", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
png("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/figures/boxplot_PRS_SCZ2_1.png", type = "cairo", width = 10, height = 6, units = 'in', res = 300)
ggplot(prs.scz2.1, aes(x = Dx, y = PRS, fill = Dx)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(col[3], col[6], col[11], col[17]), name = "Diagnosis") +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("p-value derived effect prs") +
  theme_grey(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.75), axis.title.x = element_blank(), axis.title.y = element_text(size=20))
dev.off()

# change prs.scz2 from data.table to data.frame for releveling the Dx column
prs.scz2.2 <- as.data.frame(prs.scz2.2)
prs.scz2.2$Dx <- as.factor(prs.scz2.2$Dx)
prs.scz2.2$Dx <- factor(prs.scz2.2$Dx, levels = c("CTL", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
png("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/figures/boxplot_PRS_SCZ2_2.png", type = "cairo", width = 10, height = 6, units = 'in', res = 300)
ggplot(prs.scz2.2, aes(x = Dx, y = PRS, fill = Dx)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(col[3], col[6], col[11], col[17]), name = "Diagnosis") +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("p-value derived effect prs") +
  theme_grey(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.75), axis.title.x = element_blank(), axis.title.y = element_text(size=20))
dev.off()

# change prs.scz2 from data.table to data.frame for releveling the Dx column
prs.scz2.3 <- as.data.frame(prs.scz2.3)
prs.scz2.3$Dx <- as.factor(prs.scz2.3$Dx)
prs.scz2.3$Dx <- factor(prs.scz2.3$Dx, levels = c("CTL", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
png("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/figures/boxplot_PRS_SCZ2_3.png", type = "cairo", width = 10, height = 6, units = 'in', res = 300)
ggplot(prs.scz2.3, aes(x = Dx, y = PRS, fill = Dx)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(col[3], col[6], col[11], col[17]), name = "Diagnosis") +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("p-value derived effect prs") +
  theme_grey(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.75), axis.title.x = element_blank(), axis.title.y = element_text(size=20))
dev.off()

# change prs.scz2 from data.table to data.frame for releveling the Dx column
prs.scz2.4 <- as.data.frame(prs.scz2.4)
prs.scz2.4$Dx <- as.factor(prs.scz2.4$Dx)
prs.scz2.4$Dx <- factor(prs.scz2.4$Dx, levels = c("CTL", "MDD", "BPD", "SCZ"))

# Boxplots of SCZ2 PRS
png("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/figures/boxplot_PRS_SCZ2_4.png", type = "cairo", width = 10, height = 6, units = 'in', res = 300)
ggplot(prs.scz2.4, aes(x = Dx, y = PRS, fill = Dx)) + 
  geom_boxplot() + 
  scale_fill_manual(values = c(col[3], col[6], col[11], col[17]), name = "Diagnosis") +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("p-value derived effect prs") +
  theme_grey(base_size = 18) +
  theme(plot.title = element_text(hjust = 0.75), axis.title.x = element_blank(), axis.title.y = element_text(size=20))
dev.off()

# SCATTERPLOTS
# Normal PRS SNP weights vs. eQTL smoothed PRS SNP weights
prs.scz2.weights <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, select = c("sid", "raw_beta", "ldpred_beta"))
prs.scz2.weights.smoothed <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed3_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, select = c("sid", "raw_beta", "ldpred_beta"))

# BOXPLOT: SNP weights without eQTL vs. SNP with eQTL
smoothed <- data.frame()
for (i in seq_along(prs.scz2.weights.smoothed$sid)) {
  if(prs.scz2.weights.smoothed$sid[i] == prs.scz2.weights$sid[i] & prs.scz2.weights.smoothed$ldpred_beta[i] != prs.scz2.weights$ldpred_beta[i]){
    smoothed <- rbind(smoothed, prs.scz2.weights[i, ])
  }
  if((i %% 100000) == 0){
    print(i)
  }
}
head(smoothed)
not.smoothed <- anti_join(prs.scz2.weights, smoothed) #setdiff(prs.scz2.weights, smoothed)

# generate new data.frame with extra column indicating "smoothed" and "not smoothed"
smoothed$eQTL <- "Yes"
not.smoothed$eQTL <-"No"
data <- rbind(smoothed, not.smoothed)

ggplot(data, aes(x = eQTL, y = raw_beta, fill = eQTL)) + 
  geom_boxplot() +
  scale_fill_manual(values = c(col[3], col[6])) +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("raw beta") +
  theme_grey(base_size = 18)

# SCATTERPLOT smoothed vs. not smoothed
colnames(prs.scz2.weights.smoothed) <- c("sid", "raw_beta_smoothed", "ldpred_beta_smoothed")
prs.scz2 <- merge(prs.scz2.weights, prs.scz2.weights.smoothed, by = "sid") # merge the two data.sets

# scatterplot unsmoothed (raw beta vs LDpred beta values)
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/Scatterplot_LDpred_beta_vs_raw_beta_unsmoothed.pdf")
ggplot(prs.scz2.weights, aes(x = raw_beta, y = ldpred_beta)) + geom_point(alpha = 0.3) + xlab("raw beta value") + ylab("LDpred beta value")
dev.off()

# scatterplot smoothed (raw beta vs LDpred beta)
pdf("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/Scatterplot_LDpred_beta_vs_raw_beta_smoothed.pdf")
ggplot(prs.scz2.weights.smoothed, aes(x = raw_beta_smoothed, y = ldpred_beta_smoothed)) + geom_point(alpha = 0.3) + xlab("raw beta value") + ylab("reweighted LDpred beta value")
dev.off()

# session info
sessionInfo()

