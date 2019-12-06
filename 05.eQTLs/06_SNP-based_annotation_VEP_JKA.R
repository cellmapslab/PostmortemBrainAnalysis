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
  library(reshape2)
  library(rtracklayer)
})

# function to generate an eQTL GRanges object (with MAF bins)
# @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
# @param mdata Merging data (e.g. Background SNPs)
# @param ofile Name for the rds output file of the generated GRanges object 
# GenerateGrangesObject <- function(input, mdata, ofile){
#   # merge with MAF info
#   input <- merge(input, mdata, by = "SNP")
#   # generate GRanges object
#   input.gr <- GRanges(seqnames = input$CHR, 
#                       ranges = IRanges(start = as.numeric(as.character(input$POS)), 
#                                        end = as.numeric(as.character(input$POS))),
#                       snp_id = input$SNP, 
#                       bin = input$bin)
#   # remove duplicates
#   input.gr <- input.gr[!duplicated(input.gr)] 
#   # save GRanges object to an rds file
#   saveRDS(input.gr, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/", ofile))
#   # return finalized GRanges object
#   return(input.gr)
# }

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

### Public data ###
##VEP data
# library(rtracklayer)
# library(data.table) 
# library(reshape2)
# vep=fread("WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt", data.table = F)
# vep$SNP.Position =as.numeric(sapply(as.character(vep$SNP),function(x)unlist(strsplit(x,"_"))[2]))
# vep$SNP.Chr=as.character(sapply(as.character(vep$SNP),function(x)unlist(strsplit(x,"_"))[1]))
# vep.melt = melt(vep,id.vars= c("SNP", "SNP.Chr", "SNP.Position"))
# vepGR = makeGRangesFromDataFrame(
#   data.frame(
#     chr=vep.melt.1$SNP.Chr,
#     start=vep.melt.1$SNP.Position,
#     end=vep.melt.1$SNP.Position,
#     code=as.character(vep.melt.1$variable)),
#   keep.extra.columns=TRUE)
# 
# #lift down to hg19
# hg38ToHg19<- import.chain("hg38ToHg19.over.chain")
# vepGR_hg19 <- unlist(liftOver(vepGR, hg38ToHg19))#10719405
# seqlevelsStyle(vepGR_hg19) <- "NCBI"
# saveRDS(vepGR_hg19, "VEP_hg19_collapsed_MAF01.rds")
vepGR_hg19 <- readRDS("/storage/groups/ccm01/datasets/2019_VEP/VEP_hg19_collapsed_MAF01.rds")
seqlevels(vepGR_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(vepGR_hg19), c("X", "Y", "M"))

### LOAD eQTL results ###
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"
#ipath=""
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene.gr <- readRDS(paste0(ipath,"pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls.transcript.gr <- readRDS(paste0(ipath,"pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls.exon.gr <- readRDS(paste0(ipath,"pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls.gene.cases.gr <- readRDS(paste0(ipath,"pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls.transcript.cases.gr <- readRDS(paste0(ipath,"pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls.exon.cases.gr <- readRDS(paste0(ipath,"pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))

### Load Background data ###
background.gr <- readRDS(paste0(ipath,"background_SNPs.rds"))
#length(background.gr@ranges)#[1] 6798048
#background.gr.vep <- subsetByOverlaps(background.gr, vepGR_hg19)
#length(background.gr.vep@ranges)# 3382170


### Annotation ###
# First, set parameters
perm <- 10000

# Vep codes
states <- names(table(vepGR_hg19$code))

# 2.4.1 annotation to ChromHMM Blood (parallelized)
cl <- makeCluster(17)
registerDoSNOW(cl)

results2 <- foreach (i = seq_along(states), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                       results <- data.frame()
                       overlap <- paste0("VEP_", states[i], "_SNP_overlap")
                       pval <- paste0("VEP_", states[i], "_SNP_pVal")
                       or <- paste0("VEP_", states[i], "_SNP_or")
                       vep.selected <- vepGR_hg19[vepGR_hg19$code == states[i]]
                       results[c(overlap, pval, or), "Gene"] <- annotatePublicData(eqtls.gene.gr,  vep.selected, perm, background.gr)
                       results[c(overlap, pval, or), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr,  vep.selected, perm, background.gr)
                       results[c(overlap, pval, or), "Transcript"] <- annotatePublicData(eqtls.transcript.gr,  vep.selected, perm, background.gr)
                       results[c(overlap, pval, or), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr,  vep.selected, perm, background.gr)
                       results[c(overlap, pval, or), "Exon"] <- annotatePublicData(eqtls.exon.gr, vep.selected, perm, background.gr)
                       results[c(overlap, pval, or), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr,  vep.selected, perm, background.gr)
                       results
                     }
stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_VEP.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")


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





### CMC (Fromer et al.) gene eQTLs
cmc_gene <- data.table::fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc_gene <- cmc_gene[cmc_gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc_gene_gr <- GenomicRanges::GRanges(seqnames = cmc_gene$SNP_location, 
                                      snp_id = cmc_gene$SNP)  # generate GRanges object 
cmc_gene_gr <- cmc_gene_gr[!duplicated(cmc_gene_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
GenomeInfoDb::seqlevelsStyle(cmc_gene_gr) <- "NCBI"
saveRDS(cmc_gene_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/CMC_adjustedSVA_eQTL_SNPs_background.rds")
rm(cmc_gene, cmc_gene_gr)

