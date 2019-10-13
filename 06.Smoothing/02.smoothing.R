# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")
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

# 1.2 Load eQTLs
# cis-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_cis_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
cis.genes.all <- data[,  c("SNP", "gene")]  # no filter ... (use all cis-eQTLs)
#cis.genes.all[duplicated(cis.genes.all$SNP), ] # get SNPs lying on multiple genes
rm(data)
# trans-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_trans_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
trans.genes.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with FDR < 0.05
rm(data)

genes.all <- rbind(cis.genes.all, trans.genes.all)
# Number of all cis-eQTLs and filtered trans-eQTLs
length(genes.all$SNP)

# 2. Preprocessings
# 2.1 Preprocess eQTLs
# load SNP position info
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/SNP_locations.txt", header = T, stringsAsFactors = FALSE, sep = "\t") # snp chr pos
# recode all SNP IDs to chr:pos 
snpspos$id <- paste0("chr", snpspos$chr, ":", snpspos$pos) 
snpspos <- snpspos[, c("snp", "id")]
head(snpspos)
# map eqtl SNP IDs to changed SNP IDs (rs ID to chr:pos ID)
genes.all <- merge(genes.all, snpspos, by.x = "SNP", by.y = "snp", all.x = T)
head(genes.all)
rm(snpspos)

# 2.2 Build eQTL matrix
unique.snps <- as.character(unique(genes.all$id))
# number of unique SNPs
length(unique.snps)
unique.genes <- as.character(unique(genes.all$gene))
# number of unique genes/exons
length(unique.genes)

# generate matrix with zeros (rows = 42049 SNPs, columns = 1006 genes)
zero_matrix <- big.matrix(length(unique.snps), length(unique.genes), init = 0, dimnames = list(unique.genes, unique.snps), backingfile = "big1.bin", descriptorfile = "big1.desc")
#options(bigmemory.allow.dimnames = TRUE) 
#colnames(zero_matrix) <- unique.genes  # name columns using SNP IDs
#row.names(zero_matrix) <- unique.snps  # name rows using gene names
zero_matrix[1:6, 1:6]

# generate dictionary like data, where each SNP as key gets a list with all found genes as value
gene_snp <- tapply(genes.all$id, genes.all$gene, list)
genes.all.list <- as.list(genes.all %>% group_by(gene) %>% summarise(id_list = list(id)))
names(genes.all.list$id_list) <- genes.all.list$gene
genes.all.list <- genes.all.list$id_list

# generate SNP x gene (42049 x 1006) matrix
for(i in seq_along(genes.all.list)){
  gene_name = names(genes.all.list)[i]
  ids = genes.all.list[gene_name][[1]]
  zero_matrix[ids, gene_name] = 1
}
#eqtl.matrix <- zero_matrix
#eqtl.matrix.transp <- t(zero_matrix)
#rm(zero_matrix)

# generate SNP x SNP eQTL based matrix
snps.smoothed <- zero_matrix %*% t(zero_matrix)
#rm(eqtl.matrix, eqtl.matrix.transp)

weights <- as.data.frame(prs_scz2[, c("sid", "ldpred_beta")]) # get SNP ID and LDpred weight
eqtl.snps <- as.data.frame(row.names(snps.smoothed))  # get eQTL SNP IDs
colnames(eqtl.snps) <- "sid"

# get only weights of SNPs with an eQTL
weights.eqtl.snps <- join(eqtl.snps, weights, by = "sid")
weights.eqtl.snps[is.na(weights.eqtl.snps$ldpred_beta), "ldpred_beta"] <- 0  # change all NA's into 0
identical(as.character(weights.eqtl.snps$sid), row.names(snps.smoothed))  # check if order is identical (needs to be TRUE)
weights.eqtl.snps.matrix <- as.matrix(weights.eqtl.snps$ldpred_beta)  # transform into matrix
row.names(weights.eqtl.snps.matrix) <- weights.eqtl.snps$sid

# calculate eQTL based weights
weights.smoothed <- t(weights.eqtl.snps.matrix) %*% snps.smoothed

# 1. Version: change the values in the weights file and leave the rest as it is
prs_scz2_1 <- prs_scz2
for(i in seq_along(weights.smoothed)){
  id <- colnames(weights.smoothed)[i]
  if(nrow(prs_scz2_1[prs_scz2_1$sid == id, ]) != 0){
    prs_scz2_1[prs_scz2_1$sid == id, "ldpred_beta"] = as.numeric(weights.smoothed[i])
  }
}
head(prs_scz2_1)
# save new LDpred weights file
write.table(prs_scz2_1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed1_ALL_cis_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 2. Version: change the values in the weights file and set the rest of the weights to 0
prs_scz2_2 <- prs_scz2
prs_scz2_2$ldpred_beta <- 0.0
for(i in seq_along(weights.smoothed)){
  id <- colnames(weights.smoothed)[i]
  if(nrow(prs_scz2_2[prs_scz2_2$sid == id, ]) != 0){
    prs_scz2_2[prs_scz2_2$sid == id, "ldpred_beta"] <- as.numeric(weights.smoothed[i])
  }
}
head(prs_scz2_2)
# save new LDpred weights file
write.table(prs_scz2_2, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed2_ALL_cis_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)




#############################
### TEST: ABSOLUTE VALUES ###
#############################
# load SCZ2 PRS weights
prs_scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, data.table = F)
# get absolute value of LDpred weights
prs_scz2$ldpred_beta_new <- abs(prs_scz2$ldpred_beta)
# remove old LDpred values
prs_scz2$ldpred_beta <- NULL
# change colnames
colnames(prs_scz2) <- c("chrom", "pos", "sid", "nt1", "nt2", "raw_beta", "ldpred_beta")
# save new LDpred weights file
write.table(prs_scz2, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_absolute_weights_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 3. Run Validation (get score for each patient)
#/storage/groups/ccm01/tools/ldpred-master/LDpred.py score --gf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_scz2_filtered --rf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_absolute_weights --f 0.3 --pf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/phenotypes.txt --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_absolute_weights_PRS &> /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/validate.log



# load SCZ2 PRS weights
prs_scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, data.table = F)
# change directions (+ to - and vice versa) of LDpred weights
prs_scz2$ldpred_beta_new <- prs_scz2$ldpred_beta * (-1)
# remove old LDpred values
prs_scz2$ldpred_beta <- NULL
# change colnames
colnames(prs_scz2) <- c("chrom", "pos", "sid", "nt1", "nt2", "raw_beta", "ldpred_beta")
# save new LDpred weights file
write.table(prs_scz2, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_reverse_weights_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 3. Run Validation (get score for each patient)
#/storage/groups/ccm01/tools/ldpred-master/LDpred.py score --gf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_scz2_filtered --rf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_reverse_weights --f 0.3 --pf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/phenotypes.txt --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_reverse_weights_PRS &> /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/validate_reverse.log

