# Use packages from own library (R version 3.6)
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# First, install all necessary packages as well as lassosum (only done once in the beginning)
#install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies=TRUE)
# Then, download lassosum in the terminal using following command:
# wget https://github.com/tshmak/lassosum/releases/download/v0.4.4/lassosum_0.4.4.tar.gz
#install.packages("/home/icb/karolina.worf/R/lassosum_0.4.4.tar.gz", repos = NULL, type = "source")

# 1. Load packages (for own R installation - version 3)
library(RcppArmadillo) # for colorblind color palette
library(data.table)
library(Matrix)
library(lassosum)
library(rcompanion)
library(parallel)
library(ggplot2)
library(oligo)
library(viridis)

# function, to run lassosum for a given summary statistics file
# @param ss_p    p-values given in the summary statistics file
# @param ss_or   odds ratios given in the summary statistics file
# @param ss_chr  chromosomes given in the summary statistics file
# @param ss_pos  base pair positions given in the summary statistics file
# @param ss_a1   first allele given in the summary statistics file
# @param ss_a2   second allele given in the summary statistics file
RunLassosum <- function(ss_p, ss_or, ss_chr, ss_pos, ss_a1, ss_a2){
  ### Specify the PLINK file stub of the reference panel ###
  ref_bfile <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_cohort/1.4_dups_removed/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final"
  ### Specify the PLINK file stub of the target data ###
  target_bfile <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final"
  ### Read LD region file ###
  # Other alternatives available: type ?lassosum.pipeline for more details
  LDblocks <- "EUR.hg19" # This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome
  ### Calculate SNP-wise correlations from p-values ###
  # n is sample size
  cor <- p2cor(p = ss_p, n = 150063, sign = log(ss_or))
  ### Run lassosum in parallel ###
  cl <- makeCluster(2, type = "FORK")
  out <- lassosum.pipeline(cor = cor, chr = ss_chr, pos = ss_pos, 
                           A1 = ss_a1, A2 = ss_a2, # A2 is not required but advised
                           ref.bfile = ref_bfile, test.bfile = target_bfile, 
                           LDblocks = LDblocks, cluster = cl)
  return(out)
}

GeneratePRSWithoutCovar <- function(input, gwas, output, pheno){
  # generate PRS
  v <- validate(input, pheno = pheno) # Use defined phenotype dataframe (from above) with FID, IID and Phenotype column
  # get Diagnosis and PRS
  res <- as.data.frame(cbind(v$pheno, v$best.pgs))
  colnames(res) <- c("Dx", "PRS")
  # get Pseudo R squared and p-value
  lm.prs <- lm(res$Dx ~ res$PRS)
  r2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
  r2.2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
  sum.prs <- summary(lm.prs)
  p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
  print(paste0(gwas, ":\nNagelkerke R2 = ", r2, "\nCos and Snell R2 = ", r2.2, "\nP-value = ", p.value))
  # generate boxplot
  col <- viridis(20) # get six colors (for plots)
  prs.plot <- merge(v$results.table[, c("IID", "best.pgs")], pData[, c("IID", "Diagnosis")], by = "IID")
  cases <- prs.plot[prs.plot$Diagnosis != "CTL", ]
  cases$Diagnosis <- as.character(cases$Diagnosis)
  cases$Diagnosis <- "CASE"
  prs.plot <- rbind(prs.plot, cases)
  prs.plot$Diagnosis <- factor(prs.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
  # boxplot of PRS
  ggplot(prs.plot, aes(x = Diagnosis, y = best.pgs, fill = Diagnosis)) +
    scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
    theme_classic() +
    theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
    geom_boxplot(alpha = 0.7) + 
    scale_fill_manual(values = viridis(5), guide = F) +
    ggtitle(paste0("PRS based on ", gwas, " GWAS summary statistic")) + 
    ylab("Polygenic risk score (PRS)") 
  ggsave(output, width = 12, height = 10, units = "in")
  return(v$results.table)
}

GeneratePRSWithCovar <- function(input, gwas, output, pheno, cov){
  # generate PRS
  v <- validate(input, pheno = pheno, covar = covar) # Use defined phenotype dataframe (from above) with FID, IID and Phenotype column
  # get Diagnosis and PRS
  res <- as.data.frame(cbind(v$pheno, v$best.pgs))
  colnames(res) <- c("Dx", "PRS")
  # get Pseudo R squared and p-value
  lm.prs <- lm(res$Dx ~ res$PRS)
  r2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"]
  r2.2 <- nagelkerke(lm.prs)$Pseudo.R.squared.for.model.vs.null["Cox and Snell (ML)", "Pseudo.R.squared"]
  sum.prs <- summary(lm.prs)
  p.value <- pf(sum.prs$fstatistic[1], sum.prs$fstatistic[2], sum.prs$fstatistic[3], lower.tail = FALSE)
  print(paste0(gwas, ":\nNagelkerke R2 = ", r2, "\nCos and Snell R2 = ", r2.2, "\nP-value = ", p.value))
  # generate boxplot
  col <- viridis(20) # get six colors (for plots)
  prs.plot <- merge(v$results.table[, c("IID", "best.pgs")], pData[, c("IID", "Diagnosis")], by = "IID")
  cases <- prs.plot[prs.plot$Diagnosis != "CTL", ]
  cases$Diagnosis <- as.character(cases$Diagnosis)
  cases$Diagnosis <- "CASE"
  prs.plot <- rbind(prs.plot, cases)
  prs.plot$Diagnosis <- factor(prs.plot$Diagnosis, levels = c("CTL", "CASE", "MDD", "BPD", "SCZ"))
  # boxplot of PRS
  ggplot(prs.plot, aes(x = Diagnosis, y = best.pgs, fill = Diagnosis)) +
    scale_x_discrete(labels = c("Controls", "Cases", "Major depression disorder", "Bipolar disorder", "Schizophrenia")) +
    theme_classic() +
    theme(text = element_text(size = 26), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 26), plot.title = element_text(hjust = 0.5), axis.title.x = element_blank()) +
    geom_boxplot(alpha = 0.7) + 
    scale_fill_manual(values = viridis(5), guide = F) +
    ggtitle(paste0("PRS based on ", gwas, " GWAS summary statistic")) + 
    ylab("Polygenic risk score (PRS)") 
  ggsave(output, width = 12, height = 10, units = "in")
  return(v$results.table)
}

### 1. Load phenotype file as data.frame
pheno <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/phenotypes.txt", data.table = F)
bim <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final.fam", select = c(1, 2), data.table = F)
pheno <- merge(bim, pheno, by = "V1")
colnames(pheno) <- c("FID", "IID", "Phenotype")

### 2. Reduce phenotypes to samples, for which expression data is available too (169 instead of 215 samples)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
index <- which(is.na(pData(geneSummaries)$C1))  # remove all samples with C1 = NA
dim(geneSummaries)
pData <- pData(geneSummaries)
rm(geneSummaries)
pheno <- pheno[pheno$IID %in% pData$IID, ]

### 3. Get covariates
covar <- pData[, c("FID", "IID", "Age", "Sex", "Suicide", "Type_of_Death", "C1", "C2", "C3", "C4")]

### 4. Run Lassosum ###
# splitvalidate makes the results worse...
path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/figures/"

### SCZ GWAS
scz <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs/SCZ2.sumstats")
head(scz)
scz_out <- RunLassosum(ss_p = scz$p, ss_or = scz$or, ss_chr = scz$hg19chrc, ss_pos = scz$bp, ss_a1 = scz$a1, ss_a2 = scz$a2)
scz_v <- GeneratePRSWithoutCovar(input = scz_out, gwas = "schizophrenia", output = paste0(path, "boxplot_lassosum_without_cov_PRS_SCZ2.pdf"), pheno = pheno)
write.table(scz_v, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/SCZ2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_scz2_lassosum_PRS_without_cov.best", sep = "\t", quote = FALSE)
scz_v2 <- GeneratePRSWithCovar(input = scz_out, gwas = "schizophrenia", output = paste0(path, "boxplot_lassosum_with_cov_PRS_SCZ2.pdf"), pheno = pheno, cov = covar)
write.table(scz_v2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/SCZ2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_scz2_lassosum_PRS_with_cov.best", sep = "\t", quote = FALSE)

### BPD GWAS
bpd <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs/BIP_2018.sumstats")
head(bpd)
bpd_out <- RunLassosum(ss_p = bpd$P, ss_or = bpd$OR, ss_chr = bpd$CHR, ss_pos = bpd$BP, ss_a1 = bpd$A1, ss_a2 = bpd$A2)
bpd_v <- GeneratePRSWithoutCovar(input = bpd_out, gwas = "bipolar disorder", output = paste0(path, "boxplot_lassosum_without_cov_PRS_BPD.pdf"), pheno = pheno)
write.table(bpd_v, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_lassosum_PRS_without_cov.best", sep = "\t", quote = FALSE)
bpd_v2 <- GeneratePRSWithCovar(input = bpd_out, gwas = "bipolar disorder", output = paste0(path, "boxplot_lassosum_with_cov_PRS_BPD.pdf"), pheno = pheno, cov = covar)
write.table(bpd_v2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_lassosum_PRS_with_cov.best", sep = "\t", quote = FALSE)

### MDD GWAS
mdd <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs/MDD2_2018.sumstats")
head(mdd)
mdd_out <- RunLassosum(ss_p = mdd$P, ss_or = mdd$OR, ss_chr = mdd$CHR, ss_pos = mdd$BP, ss_a1 = mdd$A1, ss_a2 = mdd$A2)
mdd_v <- GeneratePRSWithoutCovar(input = mdd_out, gwas = "major depression disorder", output = paste0(path, "boxplot_lassosum_without_cov_PRS_MDD2.pdf"), pheno = pheno)
write.table(mdd_v, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/MDD2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_mdd2_lassosum_PRS_without_cov.best", sep = "\t", quote = FALSE)
mdd_v2 <- GeneratePRSWithCovar(input = mdd_out, gwas = "major depression disorder", output = paste0(path, "boxplot_lassosum_with_cov_PRS_MDD2.pdf"), pheno = pheno, cov = covar)
write.table(mdd_v2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/MDD2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_mdd2_lassosum_PRS_with_cov.best", sep = "\t", quote = FALSE)

### CD GWAS
cd <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs/Cross_disorder_hg19.sumstats")
cd$hg19chr <- unlist(lapply(cd$hg19chr, function(x) strsplit(x, "chr")[[1]][2])) # remove chr from chromosome name
head(cd)
cd_out <- RunLassosum(ss_p = cd$p, ss_or = cd$or, ss_chr = cd$hg19chr, ss_pos = cd$bp, ss_a1 = cd$a1, ss_a2 = cd$a2)
cd_v <- GeneratePRSWithoutCovar(input = cd_out, gwas = "schizophrenia", output = paste0(path, "boxplot_lassosum_without_cov_PRS_CD.pdf"), pheno = pheno)
write.table(cd_v, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/Cross_Disorder/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_cd_lassosum_PRS_without_cov.best", sep = "\t", quote = FALSE)
cd_v2 <- GeneratePRSWithCovar(input = cd_out, gwas = "schizophrenia", output = paste0(path, "boxplot_lassosum_with_cov_PRS_CD.pdf"), pheno = pheno, cov = covar)
write.table(cd_v2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/lassosum/Cross_Disorder/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_cd_lassosum_PRS_with_cov.best", sep = "\t", quote = FALSE)
