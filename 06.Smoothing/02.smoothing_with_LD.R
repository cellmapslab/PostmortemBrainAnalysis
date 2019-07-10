# Load locally installed packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
setwd("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing")

# Load packages
library(data.table)
library(bigmemory)
library(plyr)
library(dplyr)
library(magrittr)

# 1. Load data
# 1.1 Load PRS weights
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, data.table = F)

# 1.2 Load eQTLs
# cis-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_SCZ_CTL_matrixEQTL_cis_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
cis.genes.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with p-value < 0.01
rm(data)
# trans-eQTLs
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_SCZ_CTL_matrixEQTL_trans_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
trans.genes.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with p-value < 0.01
rm(data)

# 2. Preprocessing
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

prs.scz2.weights <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, select = c("sid", "ldpred_beta"))
snps.with.eqtls <- prs.scz2.weights[prs.scz2.weights$sid %in% genes.all$id, ]
snps.without.eqtls <- prs.scz2.weights[!(prs.scz2.weights$sid %in% genes.all$id), ]

# generate new data.frame with extra column indicating "smoothed" and "not smoothed"
snps.with.eqtls$eQTL <- "Yes"
snps.without.eqtls$eQTL <-"No"
data <- rbind(snps.with.eqtls, snps.without.eqtls)

ggplot(data, aes(x = eQTL, y = ldpred_beta, fill = eQTL)) + 
  geom_boxplot() +
  scale_fill_manual(values = c(col[3], col[6])) +
  ggtitle("PRS based on SCZ2 summary statistic") + 
  ylab("LDpred weights") +
  theme_grey(base_size = 18)


# 2.2 Check, for which eQTL SNP we do NOT have a LDPred weight given (are present in the LDpred weights file)
overlap <- genes.all[!(genes.all$id %in% prs.scz2$sid), "id"]
non.overlap <- unique(genes.all[!(genes.all$id %in% prs.scz2$sid), "id"])

# 2.3 For all non overlapping SNPs get the corresponding LD SNPs (which are present in/overlapping with the LDpred weight SNPs)
# 2.3.1 Run PLINK to get LD SNPs: 
# /storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --r2 --ld-window-r2 0.8 --out /storage/scratch/users/karolina.worf/MPIP/06.Smoothing/input/LD_SNPs
# 2.3.2 load LD SNP file (only chromosome and position) and transform for the two SNPs the IDs as chr:pos
ld.snps <- fread("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/input/LD_SNPs.ld", header = T, stringsAsFactors = F, data.table = F, select = c("CHR_A", "BP_A", "CHR_B", "BP_B"))
ld.snps$SNP_A <- paste("chr", ld.snps$CHR_A, ":", ld.snps$BP_A, sep = "")
ld.snps$SNP_B <- paste("chr", ld.snps$CHR_B, ":", ld.snps$BP_B, sep = "")
# 2.3.3 reduce the data to only the SNP IDs which do NOT overlap between eQTLs and LDpred weights
ld.snps <- ld.snps[ld.snps$SNP_A %in% non.overlap$id | ld.snps$SNP_B %in% non.overlap$id, c("SNP_A", "SNP_B")]
# 2.3.4 generate list of lists, where each list name is one of the not overlapping SNP IDs (chr:pos) and the list includes all LD SNPs of this SNP ID
ld.snps.changed <- data.frame(SNP_A = ld.snps$SNP_B, SNP_B = ld.snps$SNP_A) # change SNP_A to SNP_B and wise vera
ld.snps <- rbind(ld.snps, ld.snps.changed) # concatenate the two data.frames
rm(ld.snps.changed)
# generate dictionary like data, where each non overlapping SNP as key gets a list with all found LD SNPs as value
ld.snps.list <- as.list(ld.snps %>% group_by(SNP_A) %>% summarise(SNP_B = list(SNP_B)))
names(ld.snps.list$SNP_B) <- ld.snps.list$SNP_A
ld.snps.list <- ld.snps.list$SNP_B
head(ld.snps.list)
# 2.3.5 check LD SNPs for overlap to SNPs of the LDpred weight file and if multiple hits save the SNP with the best LDpred weight
ld.snps.hits <- data.frame()
for (i in seq_along(ld.snps.list)) {
  count = 0
  hit = data.frame()
  for(j in ld.snps.list[i][[1]]) {
    if(j %in% prs.scz2$sid){
      count = count + 1
      if(count > 1 && hit$ldpred_beta < prs.scz2[prs.scz2$sid == j, "ldpred_beta"]){
        hit = prs.scz2[prs.scz2$sid == j, c("sid", "ldpred_beta")]
      } else {
        hit = prs.scz2[prs.scz2$sid == j, c("sid", "ldpred_beta")]
      }
    }
  }
  if(count > 0){
    snps <- data.frame(eQTL_SNP = names(ld.snps.list[i]), LD_SNP = hit$sid)
    ld.snps.hits <- rbind(ld.snps.hits, snps)
  }
  if((i %% 1000) == 0){
    print(i)
  }
}
head(ld.snps.hits)
# 2.3.6 Change eQTL SNPs not overlapping with LDpred weight SNPs with the LD SNPs (if available)
ld.snps.hits[, ] <- lapply(ld.snps.hits[, ], as.character)
head(genes.all)
for(i in seq_along(ld.snps.hits$eQTL_SNP)){
  genes.all[genes.all$id == ld.snps.hits[i, "eQTL_SNP"], "id"] <- ld.snps.hits[i, "LD_SNP"]
}
head(genes.all)
# 2.3.7 Reduce genes.all to only SNPs also present in the LDpred weights
genes.all <- genes.all[genes.all$id %in% prs.scz2$sid, ]

# 2.2 Build eQTL matrix
unique.snps <- as.character(unique(genes.all$id))
unique.genes <- as.character(unique(genes.all$gene))
# generate matrix with zeros (rows = 16565 SNPs, columns = 785 genes)
zero_matrix <- matrix(0, length(unique.snps), length(unique.genes)) 
colnames(zero_matrix) <- unique.genes  # name columns using SNP IDs
row.names(zero_matrix) <- unique.snps  # name rows using gene names
zero_matrix[1:6, 1:6]

# generate dictionary like data, where each SNP as key gets a list with all found genes as value
#gene_snp <- tapply(genes.all$id, genes.all$gene, list)
genes.all.list <- as.list(genes.all %>% group_by(gene) %>% summarise(id_list = list(id)))
names(genes.all.list$id_list) <- genes.all.list$gene
genes.all.list <- genes.all.list$id_list

# generate SNP x gene (42049 x 1006) matrix
for(i in seq_along(genes.all.list)){
  gene_name = names(genes.all.list)[i]
  ids = genes.all.list[gene_name][[1]]
  zero_matrix[ids, gene_name] = 1
}
eqtl.matrix <- zero_matrix
eqtl.matrix.transp <- t(zero_matrix)
rm(zero_matrix)

# generate SNP x SNP eQTL based matrix
snps.smoothed <- eqtl.matrix %*% eqtl.matrix.transp
rm(eqtl.matrix, eqtl.matrix.transp)

weights <- as.data.frame(prs.scz2[, c("sid", "ldpred_beta")]) # get SNP ID and LDpred weight
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
prs.scz2.1 <- prs.scz2
for(i in seq_along(weights.smoothed)){
  id <- colnames(weights.smoothed)[i]
  if(nrow(prs.scz2.1[prs.scz2.1$sid == id, ]) != 0){
    prs.scz2.1[prs.scz2.1$sid == id, "ldpred_beta"] = as.numeric(weights.smoothed[i])
  }
}
head(prs.scz2.1)
# save new LDpred weights file
write.table(prs.scz2.1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed3_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 2. Version: change the values in the weights file and set the rest of the weights to 0
prs.scz2.2 <- prs.scz2
prs.scz2.2$ldpred_beta <- 0.0
for(i in seq_along(weights.smoothed)){
  id <- colnames(weights.smoothed)[i]
  if(nrow(prs.scz2.2[prs.scz2.2$sid == id, ]) != 0){
    prs.scz2.2[prs.scz2.2$sid == id, "ldpred_beta"] <- as.numeric(weights.smoothed[i])
  }
}
head(prs.scz2.2)
# save new LDpred weights file
write.table(prs.scz2.2, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_eqtl_smoothed4_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
