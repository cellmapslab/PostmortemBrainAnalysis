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

# function, to calculate the SNP annotation to a public dataset, but with respect to the MAF bins
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

# function, to calculate the annotation for one public dataset to all eQTL datasets
# @param input  Public dataset (e.g. CMC eQTLs, GWAS AD-2019)
# @param name   Name of the public dataset
# @param res    Data frame in which the results should be stored
# @return       the results data frame including the newly calculated results
RunAnnotationForAllEQTLs <- function(input, name, results){
  res <- data.frame()
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene"] <- AnnotatePublicData(gene_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript"] <- AnnotatePublicData(transcript_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon"] <- AnnotatePublicData(exon_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, input, perm, background)
  # Write single runs into a file
  write.table(res, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_", name, ".tsv"), row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
  results <- rbind(results, res)
  print(paste0("Done ", name))
  return(results)
}

# path, where all datasets are stored (as .rds files)
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"

### Load eQTL results ###
gene_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # gene based eQTLs
gene_cases_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only gene based eQTLs
transcript_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # transcript based eQTLs (similar to splice QTLs)
transcript_cases_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only transcript eQTLs (similar to splice QTLs)
exon_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # exon based eQTLs (eeQTLs)
exon_cases_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only exon based eQTLs (eeQTLs)

### Load Background data ###
background <- readRDS(paste0(ipath, "background_SNPs.rds"))


####################
### Cortex meQTL ###
####################
### Load public data: CMC ###
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


###########
### CMC ###
###########
### Load public data: CMC ###
cmc_eqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_eQTL_SNPs.rds"))  # CMC eQTL
cmc_isoqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_isoQTL_SNPs.rds"))  # CMC isoQTL

### Annotation to CMC ###
# Set parameters
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
perm <- 10000

# Run annotation for CMC (Fromer et al.) DPC gene cis-eQTLs
results <- RunAnnotationForAllEQTLs(cmc_eqtl, "CMC_eQTL_SNP", results)

# Run annotation for CMC (Fromer et al.) DPC isoform cis-eQTLs
results <- RunAnnotationForAllEQTLs(cmc_isoqtl, "CMC_isoQTL_SNP", results)

# Write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_annotation_results_10000_permutations_CMC.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")


############
### GTEx ###
############
### Load public data: GTEx ###
gtex_wb_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))  # GTEx Whole Blood eQTL
gtex_wb_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))  # GTEx Whole Blood splice QTL
gtex_brains_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))  # GTEx All Brain regions eQTL
gtex_brains_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))  ## GTEx All Brain regions splice QTL

### Annotation to GTEx ###
# Set parameters
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
perm <- 10000

# Run annotation for GTEx whole blood gene cis-eQTLs and isoform cis-eQTL (sQTL)
results <- RunAnnotationForAllEQTLs(gtex_wb_eqtl, "GTEx_whole_blood_eQTL_SNP", results)
results <- RunAnnotationForAllEQTLs(gtex_wb_sqtl, "GTEx_whole_blood_sQTL_SNP", results)

# Run annotation for GTEx brain gene cis-eQTLs and isoform cis-eQTL (sQTL)
results <- RunAnnotationForAllEQTLs(gtex_brains_eqtl, "GTEx_brain_eQTL_SNP", results)
results <- RunAnnotationForAllEQTLs(gtex_brains_sqtl, "GTEx_brain_sQTL_SNP", results)

# Write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_annotation_results_10000_permutations_CMC.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

print(sessionInfo())

