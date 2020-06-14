# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(arules)
  library(foreach)
  library(doSNOW)
  library(S4Vectors)
})


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

### Public data ###
### VEP data ###
vepGR_hg19 <- readRDS("/storage/groups/ccm01/datasets/2019_VEP/VEP_hg19_collapsed_MAF01.rds")
seqlevels(vepGR_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(vepGR_hg19), c("X", "Y", "M"))

### LOAD eQTL results ###
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData/"
#ipath=""
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls_gene_gr <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_transcript_gr <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_exon_gr <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_gene_cases_gr <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_transcript_cases_gr <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_exon_cases_gr <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))

### Load Background data ###
background_gr <- readRDS(paste0(ipath,"background_SNPs.rds"))
length(background_gr@ranges) #[1] 6,798,048
background_gr <- subsetByOverlaps(background_gr, vepGR_hg19)
length(background_gr@ranges) # 3,382,170

### Load public data
# CMC (Fromer et al.) gene eQTLs and isoQTLs
cmc_eqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_eQTL_SNPs.rds"))
cmc_isoqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_isoQTL_SNPs.rds"))

### GTEx Whole Blood eQTLs and sQTLs
gtex_wb_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
gtex_wb_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))

### GTEx Brain (all brain tissues) eQTLs and sQTLs
gtex_brains_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
gtex_brains_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))

# Vep codes
states <- names(table(vepGR_hg19$code))

### Get overlap of VEP to CMC, GTEx wb and brain
results <- data.frame()
for (state in states) {
  vep_selected <- vepGR_hg19[vepGR_hg19$code == state]
  results[substr(state, 1, nchar(state) - 2), "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, vep_selected))
  results[substr(state, 1, nchar(state) - 2), "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, vep_selected))
  results[substr(state, 1, nchar(state) - 2), "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, vep_selected))
  results[substr(state, 1, nchar(state) - 2), "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, vep_selected))
  results[substr(state, 1, nchar(state) - 2), "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, vep_selected))
  results[substr(state, 1, nchar(state) - 2), "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, vep_selected))
}

### write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/VEP/CMC_GTEx_SNP_overlap_to_VEP_all_states.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

### Annotation ###
# First, set parameters
perm <- 10000

# 2.4.1 annotation to VEP (parallelized)
cl <- makeCluster(17)
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above

results2 <- foreach (i = seq_along(states), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                       results <- data.frame()
                       overlap <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_overlap")
                       pval <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_pVal")
                       or <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_or")
                       vep_selected <- vepGR_hg19[vepGR_hg19$code == states[i]]
                       results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr,  vep_selected, perm, background_gr)
                       results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr,  vep_selected, perm, background_gr)
                       results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr,  vep_selected, perm, background_gr)
                       results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr,  vep_selected, perm, background_gr)
                       results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, vep_selected, perm, background_gr)
                       results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr,  vep_selected, perm, background_gr)
                       write.table(results, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/VEP/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_", strsplit(states[i], "_d$")[[1]],".tsv"), row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
                       results
                     }
stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/VEP/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_all_states.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

print(sessionInfo())

