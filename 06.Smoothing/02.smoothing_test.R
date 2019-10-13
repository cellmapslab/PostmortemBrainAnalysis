# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")
setwd("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing")

# Load packages
library(data.table)
library(bigmemory)
library(plyr)
library(dplyr)
library(magrittr)
library(bigmemory)

# 1. Load data
# 1.1 Load PRS weights
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_LDpred_p3.0000e-01.txt", header = T, stringsAsFactors = F, data.table = F)
#prs_bip <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIP-2018/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_filtered_LDpred_p1.0000e-02.txt", header = T, stringsAsFactors = F, data.table = F)

# 1.2 Transform PRS data into SNP-based GRanges objects
prs.scz2.gr <- GRanges(seqnames = unlist(lapply(prs.scz2$chrom, function(x) strsplit(x, "chrom_")[[1]][2])), 
                       ranges = IRanges(start = as.numeric(as.character(prs.scz2$pos)), 
                                        end = as.numeric(as.character(prs.scz2$pos))),
                       snp_id = prs.scz2$sid)
prs.scz2.gr <- prs.scz2.gr[!duplicated(prs.scz2.gr)]  # 3,273,750 SNPs


# 1.3 Load cis eQTLs of gene, transcript and exon-level
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

# 1.4 Transform eQTL data into SNP-based GRanges objects
# Gene eQTLs
eqtls.gene.gr <- GRanges(seqnames = eqtls.gene$SNP_location, snp_id = eqtls.gene$SNP, gene_name = eqtls.gene$Gene_symbol, p_value = eqtls.gene$`p-value`, fdr = eqtls.gene$FDR)
eqtls.gene.gr <- eqtls.gene.gr[!duplicated(eqtls.gene.gr)]  # 36,637 SNPs
# Gene cases eQTLs
eqtls.gene.cases.gr <- GRanges(seqnames = eqtls.gene.cases$SNP_location, snp_id = eqtls.gene.cases$SNP, gene_name = eqtls.gene.cases$Gene_symbol, p_value = eqtls.gene.cases$`p-value`, fdr = eqtls.gene.cases$FDR)
eqtls.gene.cases.gr <- eqtls.gene.cases.gr[!duplicated(eqtls.gene.cases.gr)]  # 12,605 SNPs
# Transcript eQTLs
eqtls.transcript.gr <- GRanges(seqnames = eqtls.transcript$SNP_location, snp_id = eqtls.transcript$SNP, transcript_id = eqtls.transcript$Transcript, p_value = eqtls.transcript$pvalue, fdr = eqtls.transcript$FDR)
eqtls.transcript.gr <- eqtls.transcript.gr[!duplicated(eqtls.transcript.gr)]  # 43,039 SNPs
# Transcript cases eQTLs
eqtls.transcript.cases.gr <- GRanges(seqnames = eqtls.transcript.cases$SNP_location, snp_id = eqtls.transcript.cases$SNP, transcript_id = eqtls.transcript.cases$Transcript, p_value = eqtls.transcript.cases$pvalue, fdr = eqtls.transcript.cases$FDR)
eqtls.transcript.cases.gr <- eqtls.transcript.cases.gr[!duplicated(eqtls.transcript.cases.gr)]  # 16,237 SNPs
# Exon eQTLs
eqtls.exon.gr <- GRanges(seqnames = eqtls.exon$SNP_location, snp_id = eqtls.exon$SNP, exon_id = eqtls.exon$Exon, p_value = eqtls.exon$pvalue, fdr = eqtls.exon$FDR)
eqtls.exon.gr <- eqtls.exon.gr[!duplicated(eqtls.exon.gr)]  # 190,570 SNPs
# Exon cases eQTLs
eqtls.exon.cases.gr <- GRanges(seqnames = eqtls.exon.cases$SNP_location, snp_id = eqtls.exon.cases$SNP, exon_id = eqtls.exon.cases$Exon, p_value = eqtls.exon.cases$pvalue, fdr = eqtls.exon.cases$FDR)
eqtls.exon.cases.gr <- eqtls.exon.cases.gr[!duplicated(eqtls.exon.cases.gr)]  # 79,108 SNPs
# remove unnecessary data
rm(eqtls.gene, eqtls.gene.cases, eqtls.transcript, eqtls.transcript.cases, eqtls.exon, eqtls.exon.cases)


### Filter for eQTLs ###
# 1. Take only SNPs with gene cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.gene.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 2. Take only SNPs with gene cis-eQTLs of cases only
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.gene.cases.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/gene_cis-eqtl_cases_only_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 3. Take only SNPs with transcript cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.transcript.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/transcript_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 4. Take only SNPs with transcript cis-eQTLs of cases only
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.transcript.cases.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/transcript_cis-eqtl_cases_only_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 5. Take only SNPs with exon cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.exon.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/exon_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 6. Take only SNPs with exon cis-eQTLs of cases only
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.exon.cases.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/exon_cis-eqtl_cases_only_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

# 7. Take only SNPs with gene cis-eQTLs and reverse the weight (multiply with -1)
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, eqtls.gene.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
weights$ldpred_beta <-  weights$ldpred_beta * (-1)
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/negative_gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)


### Load public data
### 2.2.1 BIOS whole blood gene eQTLs
bios.gene <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/gene_level_eQTLs.txt", sep = "\t")
bios.gene[, length(bios.gene)] <- NULL # remove last column
bios.gene <- bios.gene[bios.gene$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.gene.gr <- GRanges(seqnames = bios.gene$SNPChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.gene$SNPChrPos)), 
                                         end = as.numeric(as.character(bios.gene$SNPChrPos))),
                        snp_id = bios.gene$SNPName)  # generate GRanges object
bios.gene.gr <- bios.gene.gr[!duplicated(bios.gene.gr)]  # remove duplicated rows: 5,526,636 => 1,957,252 SNPs
seqlevelsStyle(bios.gene.gr) <- "NCBI"
rm(bios.gene)

# 8. Take only SNPs with BIOS blood gene cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, bios.gene.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/bios_blood_gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

### 2.2.2 BIOS whole blood exon eQTLs
bios.exon <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/exon_level_eQTLs_independent_effects.txt", sep = "\t")
bios.exon <- bios.exon[bios.exon$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.exon.gr <- GRanges(seqnames = bios.exon$SNPChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.exon$SNPChrPos)), 
                                         end = as.numeric(as.character(bios.exon$SNPChrPos))),
                        snp_id = bios.exon$SNPName)  # generate GRanges object
bios.exon.gr <- bios.exon.gr[!duplicated(bios.exon.gr)]  # remove duplicated rows: 171,889 => 102,146 SNPs
seqlevelsStyle(bios.exon.gr) <- "NCBI"
rm(bios.exon)

# 9. Take only SNPs with BIOS blood exon cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, bios.exon.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/bios_blood_exon_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

### 2.2.3 CMC (Fromer et al.) gene eQTLs
cmc.gene <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc.gene <- cmc.gene[cmc.gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.gene.gr <- GRanges(seqnames = cmc.gene$SNP_location,
                       snp_id = cmc.gene$SNP)  # generate GRanges object 
cmc.gene.gr <- cmc.gene.gr[!duplicated(cmc.gene.gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
seqlevelsStyle(cmc.gene.gr) <- "NCBI"
rm(cmc.gene)

# 10. Take only SNPs with CMC genecis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, cmc.gene.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/cmc_gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

### 2.2.4 CMC (Fromer et al.) isoform eQTLs
cmc.iso <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt")
cmc.iso <- cmc.iso[cmc.iso$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.iso.gr <- GRanges(seqnames = cmc.iso$SNP_location,
                      snp_id = cmc.iso$SNP)  # generate GRanges object 
cmc.iso.gr <- cmc.iso.gr[!duplicated(cmc.iso.gr)]  # remove duplicated rows: 3,394,514 => 1,422,121 SNPs
seqlevelsStyle(cmc.iso.gr) <- "NCBI"
rm(cmc.iso)

# 11. Take only SNPs with CMC isoform/transcript cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, cmc.iso.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/cmc_isoform_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

### 2.2.5 GTEx Whole Blood
gtex.wb <- fread("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/Whole_Blood.v7.egenes.txt")
gtex.wb <- gtex.wb[gtex.wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
gtex.wb.gr <- GRanges(seqnames = gtex.wb$chr, 
                      ranges = IRanges(start = as.numeric(as.character(gtex.wb$pos)), 
                                       end = as.numeric(as.character(gtex.wb$pos))),
                      snp_id = gtex.wb$rs_id_dbSNP147_GRCh37p13)  # generate GRanges object only with gene ID, name and location
gtex.wb.gr <- gtex.wb.gr[!duplicated(gtex.wb.gr)]  # remove duplicated rows: 8,663 => 8,367 SNPs
rm(gtex.wb)

# 12. Take only SNPs with GTEx whole blood gene cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, gtex.wb.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/gtex_wb_gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)

### 2.2.6 GTEx Brain (all brain tissues)
brains <- c("Brain_Amygdala.v7.egenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt", "Brain_Caudate_basal_ganglia.v7.egenes.txt", "Brain_Cerebellar_Hemisphere.v7.egenes.txt", "Brain_Cerebellum.v7.egenes.txt", "Brain_Cortex.v7.egenes.txt", "Brain_Frontal_Cortex_BA9.v7.egenes.txt", "Brain_Hippocampus.v7.egenes.txt", "Brain_Hypothalamus.v7.egenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v7.egenes.txt", "Brain_Putamen_basal_ganglia.v7.egenes.txt", "Brain_Spinal_cord_cervical_c-1.v7.egenes.txt", "Brain_Substantia_nigra.v7.egenes.txt")
gtex.brains <- data.frame()
for (i in brains) {
  filename <- paste0("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/", i)
  data <- fread(filename)
  data <- data[data$qval < 0.05, ]  # filter for significant hits
  gtex.brains <- rbind(gtex.brains, data) # merge different brain files together
}
gtex.brains.gr <- GRanges(seqnames = gtex.brains$chr, 
                          ranges = IRanges(start = as.numeric(as.character(gtex.brains$pos)), 
                                           end = as.numeric(as.character(gtex.brains$pos))), 
                          snp_id = gtex.brains$rs_id_dbSNP147_GRCh37p13)  # generate GRanges object only with gene ID, name and location
gtex.brains.gr <- gtex.brains.gr[!duplicated(gtex.brains.gr)]  # remove duplicated rows: 58,589 => 43,986 SNPs
rm(brains, filename, data, i, gtex.brains)

# 13. Take only SNPs with GTEx brain gene cis-eQTLs
tmp <- prs.scz2.gr[overlapsAny(prs.scz2.gr, gtex.brains.gr)]
weights <- prs.scz2[prs.scz2$sid %in% tmp$snp_id, ]
write.table(weights, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/output/gtex_brains_gene_cis-eqtl_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)
rm(tmp, weights)


### OLD: with trans eQTLs ###
# trans-eQTLs
#data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_SCZ_CTL_matrixEQTL_trans_results.txt", header = T, stringsAsFactors = F, select = c("SNP", "gene", "p-value", "FDR"))
#trans.genes.all <- data[data$FDR < 0.05,  c("SNP", "gene")]  # filter for only eQTLs with FDR < 0.05
#rm(data)

# combine datasets
#genes.all <- rbind(cis.genes.all, trans.genes.all)
#genes.sig <- rbind(cis.genes.sig, trans.genes.all)
#rm(cis.genes.all, cis.genes.sig, trans.genes.all)

# 1. Take all significant cis-/trans-eQTLs (gene-level) and take only their weights for PRS
#weights1 <- prs_scz2[prs_scz2$sid %in% genes.sig$id, ]
#write.table(weights1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test1_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 2. Take all significant cis-/trans-eQTLs (gene-level) and multiply their weights with 2
#weights2 <- prs_scz2
#weights2[weights2$sid %in% genes.sig$id, "ldpred_beta"] <- weights2[weights2$sid %in% genes.sig$id, "ldpred_beta"] * 2
#write.table(weights2, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test2_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 3. Take ALL cis-eQTLs and significant trans-eQTLs (gene-level) and take only their weights for PRS
#weights3 <- prs_scz2
#weights3[!(weights3$sid %in% genes.all$id), "ldpred_beta"] <- 0
#weights3$ldpred_beta <- sample(weights3$ldpred_beta)
##weights3 <- prs_scz2[prs_scz2$sid %in% genes.all$id, ]
#write.table(weights3, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test3_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 4. Take ALL cis-eQTLs and significant trans-eQTLs (gene-level) and multiply their weights with 2 
#weights4 <- prs_scz2
#weights2[weights2$sid %in% genes.all$id, "ldpred_beta"] <- weights2[weights2$sid %in% genes.all$id, "ldpred_beta"] * 2
#write.table(weights4, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test4_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

# 5. Take all significant cis-eQTLs (exon-level) and take only their weights
#exons.sig <- merge(cis.exons.sig, snpspos, by.x = "SNP", by.y = "snp")
#head(exons.sig)
#weights1 <- prs_scz2[prs_scz2$sid %in% exons.sig$id, ]
#write.table(weights1, "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/output/test5_smoothed_LDpred_p3.0000e-01.txt", row.names = F, col.names = T, sep = "    ", quote = F)

