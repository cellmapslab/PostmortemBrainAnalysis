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

# function to generate an eQTL GRanges object (with MAF bins)
# @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
# @param mdata Merging data (e.g. Background SNPs)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObject <- function(input, mdata, ofile){
  # merge with MAF info
  input <- merge(input, mdata, by = "SNP")
  # generate GRanges object
  input.gr <- GRanges(seqnames = input$CHR, 
                      ranges = IRanges(start = as.numeric(as.character(input$POS)), 
                                       end = as.numeric(as.character(input$POS))),
                      snp_id = input$SNP, 
                      bin = input$bin)
  # remove duplicates
  input.gr <- input.gr[!duplicated(input.gr)] 
  # save GRanges object to an rds file
  saveRDS(input.gr, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/", ofile))
  # return finalized GRanges object
  return(input.gr)
}

# function, to calculate the gene annotation to a public dataset, but with respect to the MAF bins
# @param own Own cis-eQTL data set as GRanges object
# @param public A public cis-eQTL data set as GRanges object
# @param nperm Number of permutations (usually a 1000 times)
# @param background All background genes coming from the gene expression (here: 14,048 unqiue genes)
annotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  background.bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own.bin.lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  # resample nperm times
  resampling <- lapply(1:nperm, function(x)
    sum(sapply(1:11, function(x)
      sum(
        overlapsAny(sample(background.bins[[x]], own.bin.lengths[x]), public)
       ))))
  # IF 0 => actually smaller than 1e^(-5/6)
  p.value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } 
  else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p.value), as.character(or))) # return the resampling p-value
}

### Load Background data ###
# We first need to get MAF values for each SNP using PLINK
# /storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --freq --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/annotation/pmb_info06_all_chr_maf0.05_merged_maf_information
# Now read in the MAF values for each SNP
mafs <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/annotation/pmb_info06_all_chr_maf0.05_merged_maf_information.frq", data.table=F)

# generate bins of MAF in 0.05 steps
boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
mafs$bin <- discretize(mafs$MAF, method = "fixed", breaks = boundaries, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# load background SNPs
snps <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos

### LOAD eQTL results ###
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

### SNP-based anotation ###
# Generate background SNP data as GRanges object ###
background.gr <- GenerateGrangesObject(input = snps, 
                                       mdata = mafs, 
                                       ofile = "background_SNPs_granges_object.rds")
length(background.gr$snp_id) # 6,798,048
rm(snps)

# Generate merging file for eQTL data from background SNPs #
background <- as.data.frame(background.gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")
head(background)

### Transform eQTL data into SNP-based GRanges objects ###
# Gene eQTLs
eqtls.gene.gr <- GenerateGrangesObject(input = eqtls.gene,
                                       mdata = background,
                                       ofile = "pmb_cis-eQTL_gene_hits_cases_vs_controls_granges_object.rds")
length(eqtls.gene.gr$snp_id) # 52,496 SNPs
# Gene cases eQTLs
eqtls.gene.cases.gr <- GenerateGrangesObject(input = eqtls.gene.cases,
                                             mdata = background,
                                             ofile = "pmb_cis-eQTL_gene_hits_cases_only_granges_object.rds")
length(eqtls.gene.cases.gr$snp_id) # 18,944 SNPs
# Transcript eQTLs
eqtls.transcript.gr <- GenerateGrangesObject(input = eqtls.transcript,
                                             mdata = background,
                                             ofile = "pmb_cis-eQTL_transcript_hits_cases_vs_controls_granges_object.rds")
length(eqtls.transcript.gr$snp_id) # 112,861 SNPs
# Transcript cases eQTLs
eqtls.transcript.cases.gr <- GenerateGrangesObject(input = eqtls.transcript.cases,
                                                   mdata = background,
                                                   ofile = "pmb_cis-eQTL_transcript_hits_cases_only_granges_object.rds")
length(eqtls.transcript.cases.gr$snp_id) # 40,791 SNPs
# Exon eQTLs
eqtls.exon.gr <- GenerateGrangesObject(input = eqtls.exon,
                                       mdata = background,
                                       ofile = "pmb_cis-eQTL_exon_hits_cases_vs_controls_granges_object.rds")
length(eqtls.exon.gr$snp_id) # 190,313 SNPs
# Exon cases eQTLs
eqtls.exon.cases.gr <- GenerateGrangesObject(input = eqtls.exon.cases,
                                             mdata = background,
                                             ofile = "pmb_cis-eQTL_exon_hits_cases_only_granges_object.rds")
length(eqtls.exon.cases.gr$snp_id) # 79,524 SNPs
rm(eqtls.gene, eqtls.transcript, eqtls.exon, eqtls.gene.cases, eqtls.transcript.cases, eqtls.exon.cases, background)

### Public data ###
### ChromHMM ###
# read data
hmm <- readRDS("/storage/groups/ccm01/datasets/2019_ChromHMM/chromHMM_all_states.Rds")
# remove chromosome X, Y and M
seqlevels(hmm, pruning.mode="coarse") <- setdiff(seqlevels(hmm), c("chrX", "chrY", "chrM"))
# load ChromHMM data information
hmm_info <- read.csv2("/storage/groups/ccm01/datasets/2019_ChromHMM/epigenomes.tsv", header = T, stringsAsFactors = F, sep = "\t")

# get ChromHMM brain data --> Brain_Dorsolateral_Prefrontal_Cortex
#hmm_brain <- hmm[hmm$code == "E073"]
#seqlevelsStyle(hmm_brain) <- "NCBI"

# get ChromHMM brain data
# Blood has multiple  cell lines --> filter for Group == "Blood & T-cell"
brain <- hmm_info[hmm_info$GROUP == "Brain", "EID"]  # get IDs for all brain related ChromHMM data
hmm_brain <- hmm[hmm$code %in% brain, ]  # get brain data
seqlevelsStyle(hmm_brain) <- "NCBI"
rm(hmm, hmm_info)

# Test frequency of eQTLs on chromosomes
#eqtls.chr <- as.data.frame(table(seqnames(eqtls.transcript.cases.gr)))
#eqtls.hmm <- as.data.frame(table(seqnames(hmm_brain)))
#chromosomes <- merge(eqtls.chr, eqtls.hmm, by = "Var1")
#wilcox.test(chromosomes$Freq.x, chromosomes$Freq.y, alternative = "two.sided")

### Annotation ###
# First, set parameters
perm <- 10000

# ChromHMM states
states <- names(table(hmm_brain$name))

# 2.4.1 annotation to ChromHMM Brain (parallelized)
cl <- makeCluster(15)
registerDoSNOW(cl)

results2 <- foreach (i = seq_along(states), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
  results <- data.frame()
  overlap <- paste0("ChromHMM_brain_", states[i], "_SNP_overlap")
  pval <- paste0("ChromHMM_brain_", states[i], "_SNP_pVal")
  or <- paste0("ChromHMM_brain_", states[i], "_SNP_or")
  hmm_brain.selected <- hmm_brain[hmm_brain$name == states[i]]
  results[c(overlap, pval, or), "Gene"] <- annotatePublicData(eqtls.gene.gr, hmm_brain.selected, perm, background.gr)
  results[c(overlap, pval, or), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, hmm_brain.selected, perm, background.gr)
  results[c(overlap, pval, or), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, hmm_brain.selected, perm, background.gr)
  results[c(overlap, pval, or), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, hmm_brain.selected, perm, background.gr)
  results[c(overlap, pval, or), "Exon"] <- annotatePublicData(eqtls.exon.gr, hmm_brain.selected, perm, background.gr)
  results[c(overlap, pval, or), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, hmm_brain.selected, perm, background.gr)
  results
}

stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

### Load public data
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"

# CMC (Fromer et al.) gene eQTLs and isoQTLs
cmc_eqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_eQTL_SNPs.rds"))
cmc_isoqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_isoQTL_SNPs.rds"))

### GTEx Whole Blood eQTLs and sQTLs
gtex_wb_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
gtex_wb_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))

### GTEx Brain (all brain tissues) eQTLs and sQTLs
gtex_brains_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
gtex_brains_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))

### Get overlap to CMC, GTEx wb and brain
results <- data.frame()
for (state in states) {
  hmm_brain_selected <- hmm_brain[hmm_brain$name == state]
  results[strsplit(state, "_d")[[1]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, hmm_brain_selected))
  results[strsplit(state, "_d")[[1]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, hmm_brain_selected))
  results[strsplit(state, "_d")[[1]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, hmm_brain_selected))
  results[strsplit(state, "_d")[[1]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, hmm_brain_selected))
  results[strsplit(state, "_d")[[1]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, hmm_brain_selected))
  results[strsplit(state, "_d")[[1]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, hmm_brain_selected))
}

### write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/CMC_GTEx_SNP_overlap_to_ChromHMM.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
