# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)
library(GenomicRanges)
library(tidyr)
library(arules)
library(S4Vectors)
library(reshape2)
library(rtracklayer)

# function, to calculate the gene annotation to a public dataset, but with respect to the MAF bins
# @param own Own cis-eQTL data set as GRanges object
# @param public A public cis-eQTL data set as GRanges object
# @param nperm Number of permutations (usually a 1000 times)
# @param background All background genes coming from the gene expression (here: 14,048 unqiue genes)
# @return vector including the number of overlapping SNPs (1 = overlap), the resampling p-value (2 = p_value) and the resampling odds ratio (3 = or)
AnnotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  # resample nperm times
  resampling <- lapply(1:nperm, function(x)
    sum(sapply(1:11, function(x)
      sum(
        IRanges::overlapsAny(sample(background_bins[[x]], own_bin_lengths[x]), public)
      ))))
  # IF 0 => actually smaller than 1e^(-5/6)
  p_value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } 
  else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling overlapping SNPs, p-value and odds ratio
}

# path, where all datasets are stored (as .rds files)
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"

### Load public data ###
cmc_eqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_eQTL_SNPs.rds"))  # CMC eQTL
cmc_isoqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_isoQTL_SNPs.rds"))  # CMC isoQTL
gtex_wb_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))  # GTEx Whole Blood eQTL
gtex_wb_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))  # GTEx Whole Blood splice QTL
gtex_brains_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))  # GTEx All Brain regions eQTL
gtex_brains_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))  ## GTEx All Brain regions splice QTL

### Load eQTL results ###
gene_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_eQTL_SNPs.rds"))  # gene based eQTLs
gene_cases_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only gene based eQTLs
transcript_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # transcript based eQTLs (similar to splice QTLs)
transcript_cases_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only transcript eQTLs (similar to splice QTLs)
exon_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # exon based eQTLs (eeQTLs)
exon_cases_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only exon based eQTLs (eeQTLs)

### Load Background data ###
background <- readRDS(paste0(ipath, "background_SNPs.rds"))

### Annotation ###
# First, set parameters
results <- data.frame()
perm <- 10000

### annotation to GTEx whole blood gene cis-eQTLs
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, gtex_wb_eqtl, perm, background)
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, gtex_wb_eqtl, perm, background)
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, gtex_wb_eqtl, perm, background)
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, gtex_wb_eqtl, perm, background)
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, gtex_wb_eqtl, perm, background)
results[c("GTEx_whole_blood_eQTL_overlap", "GTEx_whole_blood_eQTL_pVal", "GTEx_whole_blood_eQTL_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, gtex_wb_eqtl, perm, background)
print("Done GTEx wb eQTL")

### annotation to GTEx whole blood isoform cis-eQTLs (sQTL)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, gtex_wb_sqtl, perm, background)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, gtex_wb_sqtl, perm, background)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, gtex_wb_sqtl, perm, background)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, gtex_wb_sqtl, perm, background)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, gtex_wb_sqtl, perm, background)
results[c("GTEx_whole_blood_sQTL_eQTL_overlap", "GTEx_whole_blood_sQTL_eQTL_pVal", "GTEx_whole_blood_sQTL_eQTL_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, gtex_wb_sqtl, perm, background)
print("Done GTEx wb sQTL")

### annotation to GTEx brain gene cis-eQTLs
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, gtex_brains_eqtl, perm, background)
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, gtex_brains_eqtl, perm, background)
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, gtex_brains_eqtl, perm, background)
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, gtex_brains_eqtl, perm, background)
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, gtex_brains_eqtl, perm, background)
results[c("GTEx_brain_eQTL_overlap", "GTEx_brain_eQTL_pVal", "GTEx_brain_eQTL_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, gtex_brains_eqtl, perm, background)
print("Done GTEx brain eQTL")

### annotation to GTEx brain isoform cis-eQTLs (sQTL)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, gtex_brains_sqtl, perm, background)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, gtex_brains_sqtl, perm, background)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, gtex_brains_sqtl, perm, background)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, gtex_brains_sqtl, perm, background)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, gtex_brains_sqtl, perm, background)
results[c("GTEx_brain_sQTL_eQTL_overlap", "GTEx_brain_sQTL_eQTL_pVal", "GTEx_brain_sQTL_eQTL_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, gtex_brains_sqtl, perm, background)
print("Done GTEx brain sQTL")

# 2.4.6 write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

### annotation to CMC (Fromer et al.) DPC gene cis-eQTLs
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, cmc_eqtl, perm, background)
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, cmc_eqtl, perm, background)
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, cmc_eqtl, perm, background)
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, cmc_eqtl, perm, background)
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, cmc_eqtl, perm, background)
results[c("CMC_eQTL_SNP_overlap", "CMC_eQTL_SNP_pVal", "CMC_eQTL_SNP_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, cmc_eqtl, perm, background)
print("Done CMC gene")

### annotation to CMC (Fromer et al.) DPC isoform cis-eQTLs
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Gene"] <- AnnotatePublicData(gene_eqtl, cmc_isoqtl, perm, background)
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, cmc_isoqtl, perm, background)
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Transcript"] <- AnnotatePublicData(transcript_eqtl, cmc_isoqtl, perm, background)
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, cmc_isoqtl, perm, background)
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Exon"] <- AnnotatePublicData(exon_eqtl, cmc_isoqtl, perm, background)
results[c("CMC_isoQTL_eQTL_SNP_overlap", "CMC_isoQTL_eQTL_SNP_pVal", "CMC_isoQTL_eQTL_SNP_odds_ratio"), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, cmc_isoqtl, perm, background)
print("Done CMC iso")

# 2.4.6 write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_CMC_GTEx.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")


