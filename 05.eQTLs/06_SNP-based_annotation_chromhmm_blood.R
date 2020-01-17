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
annotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  # resample nperm times
  resampling <- lapply(1:nperm, function(x)
    sum(sapply(1:11, function(x)
      sum(
        overlapsAny(sample(background_bins[[x]], own_bin_lengths[x]), public)
      ))))
  # IF 0 => actually smaller than 1e^(-5/6)
  p_value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } 
  else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling p-value
}


### Public data ###
### ChromHMM ###
# read data
hmm <- readRDS("/storage/groups/ccm01/datasets/2019_ChromHMM/chromHMM_all_states.Rds")
# remove chromosome X, Y and M
seqlevels(hmm, pruning.mode = "coarse") <- setdiff(seqlevels(hmm), c("chrX", "chrY", "chrM"))
# load ChromHMM data information
hmm_info <- read.csv2("/storage/groups/ccm01/datasets/2019_ChromHMM/epigenomes.tsv", header = T, stringsAsFactors = F, sep = "\t")

# get ChromHMM blood data
# Blood has multiple  cell lines --> filter for Group == "Blood & T-cell"
blood <- hmm_info[hmm_info$GROUP == "Blood & T-cell", "EID"]  # get IDs for all blood related ChromHMM data
hmm_blood <- hmm[hmm$code %in% blood, ]  # get blood data
seqlevelsStyle(hmm_blood) <- "NCBI"
rm(hmm, hmm_info)


### LOAD eQTL results ###
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls_gene_gr <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_transcript_gr <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_exon_gr <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_gene_cases_gr <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_transcript_cases_gr <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_exon_cases_gr <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))


### Load Background data ###
background_gr <- readRDS(paste0(ipath, "background_SNPs.rds"))
length(background_gr@ranges) #[1] 6,798,048
# get only SNPs overlapping with ChromHMM
background_gr <- subsetByOverlaps(background_gr, hmm_brain)
length(background_gr@ranges) # ???

# Test frequency of eQTLs on chromosomes
#eqtls_chr <- as.data.frame(table(seqnames(eqtls_transcript_cases_gr)))
#eqtls_hmm <- as.data.frame(table(seqnames(hmm_blood)))
#chromosomes <- merge(eqtls_chr, eqtls_hmm, by = "Var1")
#wilcox.test(chromosomes$Freq.x, chromosomes$Freq.y, alternative = "two.sided")

### Annotation ###
# First, set parameters
perm <- 10000

# ChromHMM states
states <- names(table(hmm_blood$name))

# 2.4.1 annotation to ChromHMM Blood (parallelized)
cl <- makeCluster(15)
registerDoSNOW(cl)

results2 <- foreach (i = seq_along(states), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
  results <- data.frame()
  overlap <- paste0("ChromHMM_blood_", states[i], "_SNP_overlap")
  pval <- paste0("ChromHMM_blood_", states[i], "_SNP_pVal")
  or <- paste0("ChromHMM_blood_", states[i], "_SNP_or")
  hmm_blood_selected <- hmm_blood[hmm_blood$name == states[i]]
  results[c(overlap, pval, or), "Gene"] <- annotatePublicData(eqtls_gene_gr, hmm_blood_selected, perm, background_gr)
  results[c(overlap, pval, or), "Gene_cases"] <- annotatePublicData(eqtls_gene_cases_gr, hmm_blood_selected, perm, background_gr)
  results[c(overlap, pval, or), "Transcript"] <- annotatePublicData(eqtls_transcript_gr, hmm_blood_selected, perm, background_gr)
  results[c(overlap, pval, or), "Transcript_cases"] <- annotatePublicData(eqtls_transcript_cases_gr, hmm_blood_selected, perm, background_gr)
  results[c(overlap, pval, or), "Exon"] <- annotatePublicData(eqtls_exon_gr, hmm_blood_selected, perm, background_gr)
  results[c(overlap, pval, or), "Exon_cases"] <- annotatePublicData(eqtls_exon_cases_gr, hmm_blood_selected, perm, background_gr)
  results
}
stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")


### Cortex meQTL ###
# Load meQTLs
cortex_meqtl <- fread("/storage/groups/ccm01/datasets/2019_cortex_meQTL/GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv")
cortex_meqtl <- cortex_meqtl[cortex_meqtl$p < 0.05, ]
cortex_meqtl_gr <- GRanges(seqnames = cortex_meqtl$snpChr, 
                       ranges = IRanges(start = as.numeric(as.character(cortex_meqtl$snpPos)), 
                                        end = as.numeric(as.character(cortex_meqtl$snpPos))),
                       snp_id = cortex_meqtl$snpRsNum)  # generate GRanges object only with gene ID, name and location
cortex_meqtl_gr <- cortex_meqtl_gr[!duplicated(cortex_meqtl_gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
seqlevelsStyle(cortex_meqtl_gr) <- "NCBI"
rm(cortex_meqtl)

results <- data.frame()
perm <- 10000

# annotate meQTLs
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene"] <- annotatePublicData(eqtls_gene_gr, cortex_meqtl_gr, perm, background_gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene_cases"] <- annotatePublicData(eqtls_gene_cases_gr, cortex_meqtl_gr, perm, background_gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript"] <- annotatePublicData(eqtls_transcript_gr, cortex_meqtl_gr, perm, background_gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript_cases"] <- annotatePublicData(eqtls_transcript_cases_gr, cortex_meqtl_gr, perm, background_gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon"] <- annotatePublicData(eqtls_exon_gr, cortex_meqtl_gr, perm, background_gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon_cases"] <- annotatePublicData(eqtls_exon_cases_gr, cortex_meqtl_gr, perm, background_gr)

# write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_meQTLs.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

