# Load locally installed packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
setwd("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing")

# Load packages
library(data.table)
library(bigmemory)
library(plyr)
library(dplyr)
library(magrittr)
library(bigmemory)

# 1. Load data
# 1.1 Load PRS weights
prs_scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, data.table = F)
#prs_bip <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIP-2018/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_filtered_LDpred_p1.0000e-02.txt", header = T, stringsAsFactors = F, data.table = F)

# 1.2 Load cis eQTLs gene-level
# cis-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_cis_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
cis.genes.all <- data[,  c("SNP", "gene")]  # no filter ... (use all cis-eQTLs)
cis.genes.sig <- data[data$FDR < 0.05, c("SNP", "gene")]  # filter for FDR < 0.05 (use all cis-eQTLs)
#cis.genes.all[duplicated(cis.genes.all$SNP), ] # get SNPs lying on multiple genes
rm(data)
# trans-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_trans_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
trans.genes.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with FDR < 0.05
rm(data)

# combine datasets
genes.all <- rbind(cis.genes.all, trans.genes.all)
genes.sig <- rbind(cis.genes.sig, trans.genes.all)
rm(cis.genes.all, cis.genes.sig, trans.genes.all)

# 1.2 Load cis eQTLs exon-level
# cis-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_cis_results_exons.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
cis.exons.all <- data[,  c("SNP", "gene")]  # no filter ... (use all cis-eQTLs)
cis.exons.sig <- data[data$FDR < 0.05, c("SNP", "gene")]  # filter for FDR < 0.05 (use all cis-eQTLs)
#cis.exons.all[duplicated(cis.exons.all$SNP), ] # get SNPs lying on multiple exons
rm(data)
# trans-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_trans_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
trans.exons.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with FDR < 0.05
rm(data)

# combine datasets
exons.all <- rbind(cis.exons.all, trans.exons.all)
exons.sig <- rbind(cis.exons.sig, trans.exons.all)
rm(cis.exons.all, cis.exons.sig, trans.exons.all)

# Preprocess eQTLs
# load SNP position info
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/SNP_locations.txt", header = T, stringsAsFactors = FALSE, sep = "\t") # snp chr pos
# recode all SNP IDs to chr:pos 
snpspos$id <- paste0("chr", snpspos$chr, ":", snpspos$pos) 
snpspos <- snpspos[, c("snp", "id")]
head(snpspos)
# map eqtl SNP IDs to changed SNP IDs (rs ID to chr:pos ID)
genes.all <- merge(genes.all, snpspos, by.x = "SNP", by.y = "snp")
genes.sig <- merge(genes.sig, snpspos, by.x = "SNP", by.y = "snp")
head(genes.all)
rm(snpspos)


# Number of all cis-eQTLs and filtered trans-eQTLs
length(genes.all$SNP)

# 1. Take all significant cis-/trans-eQTLs (gene-level) and take only their weights for PRS
weights1 <- prs_scz2[prs_scz2$sid %in% genes.sig$id, ]
write.table(weights1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test1_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 2. Take all significant cis-/trans-eQTLs (gene-level) and multiply their weights with 2
weights2 <- prs_scz2
weights2[weights2$sid %in% genes.sig$id, "ldpred_beta"] <- weights2[weights2$sid %in% genes.sig$id, "ldpred_beta"] * 2
write.table(weights2, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test2_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 3. Take ALL cis-eQTLs and significant trans-eQTLs (gene-level) and take only their weights for PRS
weights3 <- prs_scz2
weights3[!(weights3$sid %in% genes.all$id), "ldpred_beta"] <- 0
weights3$ldpred_beta <- sample(weights3$ldpred_beta)
#weights3 <- prs_scz2[prs_scz2$sid %in% genes.all$id, ]
write.table(weights3, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test3_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 4. Take ALL cis-eQTLs and significant trans-eQTLs (gene-level) and multiply their weights with 2 
weights4 <- prs_scz2
weights2[weights2$sid %in% genes.all$id, "ldpred_beta"] <- weights2[weights2$sid %in% genes.all$id, "ldpred_beta"] * 2
write.table(weights4, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test4_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 5. Take all significant cis-eQTLs (exon-level) and take only their weights
exons.sig <- merge(cis.exons.sig, snpspos, by.x = "SNP", by.y = "snp")
head(exons.sig)
weights1 <- prs_scz2[prs_scz2$sid %in% exons.sig$id, ]
write.table(weights1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test5_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

