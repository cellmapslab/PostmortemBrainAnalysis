############################
### SNP-based Annotation ###
############################
# This R-Script runs a SNP-based annotation for the eQTL data against several public datasets (VEP, ChromHMM, GTEx, CMC, meQTL, GWAS)
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args inpath=<path/to/input-files> what=<data> outpath=<path/for/output>" 08_SNP-based_annotation.R &>> SNP-based_annotation.log
# or using Rscript: Rscript --slave 08_SNP-based_annotation.R inpath=<path/to/input-files> what=<data> outpath=<path/for/output> &>> SNP-based_annotation.log

# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(arules)
  library(reshape2)
  library(rtracklayer)
  library(foreach)
  library(doSNOW)
  library(S4Vectors)
})

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to run a SNP-based annotation for all eQTL datasets for a defined public dataset, e.g. VEP, ChromHMM, GTEx, CMC, meQTL GWAS.\n\nYou can invoke the script as follows:\nRscript 08_SNP-based_annotation.R inpath=\"/home/.../to/input-files\" what=\"dataset-name\" outpath=\"path/to/store/output\" &>> SNP-based_annotation.log \n\nCommandline parameters:\n\ninpath\t\tpath to all the input .rds files including the data as a GRanges object\nwhat\t\tname of the dataset (possible names: VEP, ChromHMM, GTEx, CMC, meQTL and GWAS)\noutpath\t\tpath to store all output files into\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  # function, to calculate the SNP annotation to a public dataset with multiple categories, but with respect to the MAF bins
  # @param own          Own cis-eQTL data set as GRanges object
  # @param public       A public cis-eQTL data set as GRanges object
  # @param nperm        Number of permutations (usually a 1000 times)
  # @param background   All background SNPs of the eQTL generation overlapping with the public dataset (changes for each public dataset)
  # @return             Vector including the number of overlapping SNPs (1 = overlap), the resampling p-value (2 = p_value) and the resampling odds ratio (3 = or)
  AnnotatePublicData <- function(own, public, nperm, background){
    #print(paste0("Number of overlapping background SNPs: ", length(background@ranges)))  # print number of overlapping background SNPs
    overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
    background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
    own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
    # resample nperm times
    resampling <- lapply(1:nperm, function(x)
      sum(sapply(1:11, function(y)
        # check that the background maf bin is not empty (happens sometimes, if the overlap is empty)
        ifelse(length(background_bins[[y]]) != 0, sum(
          # take n random samples of the background (n = own bin length except background bin length is smaller, than take all background genes of this MAF bin) 
          IRanges::overlapsAny(sample(background_bins[[y]], ifelse(length(background_bins[[y]]) >= own_bin_lengths[y], own_bin_lengths[y], length(background_bins[[y]]))), public)
        ), 0))))
    # IF 0 => actually smaller than 1e^(-5/6)
    p_value <- ifelse((sum(resampling >= overlap) / nperm) == 0, paste("<=", 1/nperm, sep = " "), sum(resampling >= overlap) / nperm)
    or <- overlap/mean(unlist(resampling))
    return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling p-value
  }
  
  # function, to calculate the annotation for one public dataset to all eQTL datasets
  # @param input  Public dataset (e.g. CMC eQTLs, GWAS AD-2019)
  # @param name   Name of the public dataset
  # @param res    Data frame in which the results should be stored
  # @return       the results data frame including the newly calculated results
  RunAnnotationForAllEQTLs <- function(input, name, results, output, nperm, background){
    #print(length(background@ranges))  # print number of overlapping background SNPs
    res <- data.frame()
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, input, nperm, background)
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, input, nperm, background)
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, input, nperm, background)
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, input, nperm, background)
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, input, nperm, background)
    res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, input, nperm, background)
    # Write single runs into a file
    write.table(res, file = paste0(output, "/pmb_cis-eQTL_annotation_results_10000_permutations_", name, ".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    results <- rbind(results, res)
    print(paste0("Done ", name))
    return(results)
  }
  
  ###########
  ### PMB ###
  ###########
  ### LOAD eQTL results ###
  # cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
  eqtls_gene_gr <- readRDS(paste0(inpath, "/01_pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
  eqtls_transcript_gr <- readRDS(paste0(inpath, "/01_pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
  eqtls_exon_gr <- readRDS(paste0(inpath, "/01_pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
  eqtls_gene_cases_gr <- readRDS(paste0(inpath, "/02_pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))
  eqtls_transcript_cases_gr <- readRDS(paste0(inpath, "/02_pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))
  eqtls_exon_cases_gr <- readRDS(paste0(inpath, "/02_pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))
  print("Done loading eQTL")
  
  ### Load Background data ###
  background_gr <- readRDS(paste0(inpath, "/background_SNPs.rds"))
  length(background_gr@ranges) #[1] 6,798,048
  print("Done loading background")

  ###################
  ### Public data ###
  ###################
  ### CMC ###
  ###########
  # Load CMC (Fromer et al.) gene eQTLs and isoQTLs
  cmc_eqtl <- readRDS(paste0(inpath, "/CMC_adjustedSVA_eQTL_SNPs_fdr0.05.rds"))
  cmc_isoqtl <- readRDS(paste0(inpath, "/CMC_adjustedSVA_isoQTL_SNPs_fdr0.05.rds"))
  print("Done loading CMC")
  
  ############
  ### GTEx ###
  ############
  ### Load GTEx Whole Blood eQTLs and sQTLs
  gtex_wb_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs_fdr0.05.rds"))
  gtex_wb_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs_fdr0.05.rds"))
  print("Done loading GTEx whole blood")
  ### GTEx Brain (all brain tissues) eQTLs and sQTLs
  gtex_brains_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_eQTL_SNPs_fdr0.05.rds"))
  gtex_brains_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_sQTL_SNPs_fdr0.05.rds"))
  print("Done loading GTEx all brains")
  ### GTEx Brain (only BA9 = Broadmann area 9) eQTLs and sQTLs
  gtex_ba9_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_eQTL_SNPs_fdr0.05.rds"))
  gtex_ba9_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_sQTL_SNPs_fdr0.05.rds"))
  print("Done loading GTEx BA9")
  
  # run annotation for the dataset (what = VEP/ChromHMM_Blood/ChromHMM_Brain/ChromHMM_DLPFC/CMC/GTEx_Blood/GTEx_Brain/GTEx_BA9/meQTL/GWAS)
  if (what == "VEP") {
    # Generate output path
    vep_outpath <- paste0(outpath, "/01_VEP")
    system(paste0("mkdir -p ", vep_outpath))
    
    ### Load VEP ###
    vep <- readRDS(paste0(inpath, "/VEP_hg19_collapsed_MAF01.rds"))
   
    # Get VEP states
    vep_counts <- table(vep$code)
    states <- names(vep_counts)
    print(states)
    
    # Save states with SNP counts into a file
    vep_counts <- as.data.frame(vep_counts)
    vep_counts$Var1 <- substr(vep_counts$Var1, 1, nchar(as.character(vep_counts$Var1)) - 2)
    write.table(vep_counts, file = paste0(vep_outpath, "/VEP_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")

    # get only background SNPs overlapping with VEP
    background <- subsetByOverlaps(background_gr, vep)
    
    # Get overlap of VEP to CMC, GTEx wb, all brains and BA9
    results <- data.frame()
    for (state in states) {
      vep_selected <- vep[vep$code == state]
      results[substr(state, 1, nchar(state) - 2), "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_eQTL_BA9"] <- sum(overlapsAny(gtex_ba9_eqtl, vep_selected))
      results[substr(state, 1, nchar(state) - 2), "GTEx_sQTL_BA9"] <- sum(overlapsAny(gtex_ba9_sqtl, vep_selected))
    }
    # Write results into a file
    write.table(results, file = paste0(vep_outpath, "/CMC_GTEx_SNP_overlap_to_VEP_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
    # Annotation to VEP (parallelized)
    perm <- 10000
    cl <- makeCluster(10)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_overlap")
                           pval <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_pVal")
                           or <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_or")
                           vep_selected <- vep[vep$code == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, vep_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, vep_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, vep_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, vep_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, vep_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, vep_selected, perm, background)
                           write.table(results, file = paste0(vep_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_", strsplit(states[i], "_d$")[[1]],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file and delete unnecessary data
    write.table(results2, file = paste0(vep_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(vep_outpath, vep, states, vep_counts, background, results, results2)
    
  } else if (what == "ChromHMM_Blood") {
    # Generate output path
    hmm_blood_outpath <- paste0(outpath, "/02_ChromHMM_Blood")
    system(paste0("mkdir -p ", hmm_blood_outpath))
    
    # Load ChromHMM Blood & T-cells
    hmm_blood <- readRDS(paste0(inpath, "/ChromHMM_Blood_and_T-cells.rds"))
    
    # Get ChromHMM Blood states
    hmm_blood_counts <- table(hmm_blood$name)
    states <- names(hmm_blood_counts)
    print(states)
    
    # Save states with SNP counts into a file
    hmm_blood_counts <- as.data.frame(hmm_blood_counts)
    write.table(hmm_blood_counts, file = paste0(hmm_blood_outpath, "/ChromHMM_Blood_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    # get only background SNPs overlapping with ChromHMM Blood
    background <- subsetByOverlaps(background_gr, hmm_blood)
    
    ### Get overlap to CMC, GTEx wb, all brains and BA9
    results <- data.frame()
    for (state in states) {
      hmm_blood_selected <- hmm_blood[hmm_blood$name == state]
      results[strsplit(state, "_d")[[1]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_BA9"] <- sum(overlapsAny(gtex_ba9_eqtl, hmm_blood_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_BA9"] <- sum(overlapsAny(gtex_ba9_sqtl, hmm_blood_selected))
    }
    # Write results into a file
    write.table(results, file = paste0(hmm_blood_outpath, "/CMC_GTEx_SNP_overlap_to_ChromHMM_blood_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
    # Annotation to ChromHMM Blood (parallelized)
    perm <- 10000
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_blood_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_blood_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_blood_", states[i], "_SNP_or")
                           hmm_blood_selected <- hmm_blood[hmm_blood$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_blood_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_blood_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_blood_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_blood_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_blood_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_blood_selected, perm, background)
                           if(states[i] == "8_ZNF/Rpts"){
                             write.table(results, file = paste0(hmm_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood_8_ZNF_Rpts.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           } else {
                             write.table(results, file = paste0(hmm_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood_", states[i], ".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           }
                           results
                         }
    stopCluster(cl)
    # Write results into a file and delete unnecessary data
    write.table(results2, file = paste0(hmm_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(hmm_blood_outpath, hmm_blood, states, hmm_blood_counts, background, results, results2)
    
  } else if (what == "ChromHMM_Brain") {
    # Generate output path
    hmm_brain_outpath <- paste0(outpath, "/02_ChromHMM_Brain")
    system(paste0("mkdir -p ", hmm_brain_outpath))
    
    # Load ChromHMM all Brains
    hmm_brain <- readRDS(paste0(inpath, "/ChromHMM_all_Brains.rds"))
    
    # Get ChromHMM Brain states
    hmm_brain_counts <- table(hmm_brain$name)
    states <- names(hmm_brain_counts)
    print(states)
    
    # Save states with SNP counts into a file
    hmm_brain_counts <- as.data.frame(hmm_brain_counts)
    write.table(hmm_brain_counts, file = paste0(hmm_brain_outpath, "/ChromHMM_Brain_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")

    # get only background SNPs overlapping with ChromHMM Blood
    background <- subsetByOverlaps(background_gr, hmm_brain)
    
    ### Get overlap to CMC, GTEx wb, all brains and BA9
    results <- data.frame()
    for (state in states) {
      hmm_brain_selected <- hmm_brain[hmm_brain$name == state]
      results[strsplit(state, "_d")[[1]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_BA9"] <- sum(overlapsAny(gtex_ba9_eqtl, hmm_brain_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_BA9"] <- sum(overlapsAny(gtex_ba9_sqtl, hmm_brain_selected))
    }
    # Write results into a file
    write.table(results, file = paste0(hmm_brain_outpath, "/CMC_GTEx_SNP_overlap_to_ChromHMM_brain_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
    # Annotation to ChromHMM Brain (parallelized)
    perm <- 10000
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_brain_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_brain_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_brain_", states[i], "_SNP_or")
                           hmm_brain_selected <- hmm_brain[hmm_brain$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_brain_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_brain_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_brain_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_brain_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_brain_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_brain_selected, perm, background)
                           if(states[i] == "8_ZNF/Rpts"){
                             write.table(results, file = paste0(hmm_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain_8_ZNF_Rpts.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           } else {
                             write.table(results, file = paste0(hmm_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain_", states[i],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           }
                           results
                         }
    stopCluster(cl)
    # Write results into a file and delete unnecessary data
    write.table(results2, file = paste0(hmm_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(hmm_brain_outpath, hmm_brain, states, hmm_brain_counts, background, results, results2)
    
  } else if (what == "ChromHMM_DLPFC") {
    # Generate output path
    hmm_dlpfc_outpath <- paste0(outpath, "/02_ChromHMM_DLPFC")
    system(paste0("mkdir -p ", hmm_dlpfc_outpath))
    
    # Load ChromHMM Brain DLPFC only
    hmm_dlpfc <- readRDS(paste0(inpath, "/ChromHMM_Brain_Dorsolateral_Prefrontal_Cortex.rds"))

    # Get ChromHMM DLPFC states
    hmm_dlpfc_counts <- table(hmm_dlpfc$name)
    states <- names(hmm_dlpfc_counts)
    print(states)
    
    # Save states with SNP counts into a file
    hmm_dlpfc_counts <- as.data.frame(hmm_dlpfc_counts)
    write.table(hmm_dlpfc_counts, file = paste0(hmm_dlpfc_outpath, "/ChromHMM_DLPFC_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")

    # get only background SNPs overlapping with ChromHMM DLPFC
    background <- subsetByOverlaps(background_gr, hmm_dlpfc)
    
    ### Get overlap to CMC, GTEx wb, all brains and BA9
    results <- data.frame()
    for (state in states) {
      hmm_dlpfc_selected <- hmm_dlpfc[hmm_dlpfc$name == state]
      results[strsplit(state, "_d")[[1]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_eQTL_BA9"] <- sum(overlapsAny(gtex_ba9_eqtl, hmm_dlpfc_selected))
      results[strsplit(state, "_d")[[1]], "GTEx_sQTL_BA9"] <- sum(overlapsAny(gtex_ba9_sqtl, hmm_dlpfc_selected))
    }
    # Write results into a file
    write.table(results, file = paste0(hmm_dlpfc_outpath, "/CMC_GTEx_SNP_overlap_to_ChromHMM_DLPFC.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
    # Annotation to ChromHMM Brain DLPFC (parallelized)
    perm <- 10000
    cl <- makeCluster(8)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_DLPFC_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_DLPFC_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_DLPFC_", states[i], "_SNP_or")
                           hmm_dlpfc_selected <- hmm_dlpfc[hmm_dlpfc$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_dlpfc_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_dlpfc_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_dlpfc_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_dlpfc_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_dlpfc_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_dlpfc_selected, perm, background)
                           if(states[i] == "8_ZNF/Rpts"){
                             write.table(results, file = paste0(hmm_dlpfc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_DLPFC_8_ZNF_Rpts.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           } else {
                             write.table(results, file = paste0(hmm_dlpfc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_DLPFC_", states[i],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           }
                           results
                         }
    stopCluster(cl)
    
    # Write results into a file and delete unnecessary data
    write.table(results2, file = paste0(hmm_dlpfc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_DLPFC_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(hmm_dlpfc_outpath, hmm_dlpfc, states, hmm_dlpfc_counts, background, results, results2)
    
  } else if (what == "GTEx_Blood") {
    # Generate output path
    gtex_blood_outpath <- paste0(outpath, "/03_GTEx_Blood")
    system(paste0("mkdir -p ", gtex_blood_outpath))
    
    # write number of SNPs in GTEx Blood into a file
    states <- cbind(rbind("GTEx Blood eQTL", "GTEx Blood sQTL"), rbind(length(gtex_wb_eqtl$snp_id), length(gtex_wb_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_blood_outpath, "/GTEx_Blood_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Get only background SNPs overlapping with the whole blood GTEx eQTL dataset (load data first)
    gtex_wb <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_wb)
    # Run annotation for GTEx whole blood gene cis-eQTLs (eQTL)
    results <- RunAnnotationForAllEQTLs(gtex_wb_eqtl, "GTEx_whole_blood_eQTL_SNP", results, gtex_blood_outpath, perm, background)
    rm(background, gtex_wb) # delete unnecessary data
    
    # Get only background SNPs overlapping with the whole blood GTEx eQTL dataset (load data first)
    gtex_wb <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_wb)
    # Run annotation for GTEx whole blood isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_wb_sqtl, "GTEx_whole_blood_sQTL_SNP", results, gtex_blood_outpath, perm, background)
    rm(background, gtex_wb) # delete unnecessary data
    
    # Write results into a file and remove unnecessary data
    write.table(results, file = paste0(gtex_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_whole_blood.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(gtex_blood_outpath, gtex_wb_eqtl, gtex_wb_sqtl, states, perm, results)
    
  } else if (what == "GTEx_Brain") {
    # Generate output path
    gtex_brain_outpath <- paste0(outpath, "/03_GTEx_Brain")
    system(paste0("mkdir -p ", gtex_brain_outpath))
    
    # write number of SNPs in GTEx Brain into a file
    states <- cbind(rbind("GTEx Brain eQTL", "GTEx Brain sQTL"), rbind(length(gtex_brains_eqtl$snp_id), length(gtex_brains_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_brain_outpath, "/GTEx_Brain_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Get only background SNPs overlapping with the brain GTEx eQTL dataset (load data first)
    gtex_brains <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_brains)
    # Run annotation for GTEx brain gene cis-eQTLs (eQTL)
    results <- RunAnnotationForAllEQTLs(gtex_brains_eqtl, "GTEx_all_brains_eQTL_SNP", results, gtex_brain_outpath, perm, background)
    rm(background, gtex_brains) # delete unnecessary data
    
    # Get only background SNPs overlapping with the brain GTEx sQTL dataset (load data first)
    gtex_brains <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_brains)
    # Run annotation for GTEx brain isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_brains_sqtl, "GTEx_all_brains_sQTL_SNP", results, gtex_brain_outpath, perm, background)
    rm(background, gtex_brains) # delete unnecessary data
    
    # Write results into a file and remove unnecessary data
    write.table(results, file = paste0(gtex_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_all_brains.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(gtex_brain_outpath, gtex_brains_eqtl, gtex_brains_sqtl, states, perm, results)
    
  } else if (what == "GTEx_BA9") {
    # Generate output path
    gtex_ba9_outpath <- paste0(outpath, "/03_GTEx_BA9")
    system(paste0("mkdir -p ", gtex_ba9_outpath))
    
    # write number of SNPs in GTEx BA9 into a file
    states <- cbind(rbind("GTEx Brain BA9 eQTL", "GTEx Brain BA9 sQTL"), rbind(length(gtex_ba9_eqtl$snp_id), length(gtex_ba9_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_ba9_outpath, "/GTEx_Brain_BA9_only_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Get only background SNPs overlapping with the brain BA9 GTEx eQTL dataset (load data first)
    gtex_ba9 <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_eQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_ba9)
    # Run annotation for GTEx brain gene cis-eQTLs (eQTL)
    results <- RunAnnotationForAllEQTLs(gtex_ba9_eqtl, "GTEx_brain_BA9_eQTL_SNP", results, gtex_ba9_outpath, perm, background)
    rm(background, gtex_ba9) # delete unnecessary data
    
    # Get only background SNPs overlapping with the brain BA9 GTEx sQTL dataset (load data first)
    gtex_ba9 <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_sQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, gtex_ba9)
    # Run annotation for GTEx brain isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_ba9_sqtl, "GTEx_brain_BA9_sQTL_SNP", results, gtex_ba9_outpath, perm, background)
    rm(background, gtex_ba9) # delete unnecessary data
    
    # Write results into a file and remove unnecessary data
    write.table(results, file = paste0(gtex_ba9_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_brain_BA9.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(gtex_ba9_outpath, gtex_ba9_eqtl, gtex_ba9_sqtl, states, perm, results)
    
  } else if (what == "CMC") {
    # Generate output path
    cmc_outpath <- paste0(outpath, "/04_CMC")
    system(paste0("mkdir -p ", cmc_outpath))
    
    # write number of SNPs in CMC into a file
    states <- cbind(rbind("CMC eQTL", "CMC isoQTL"), rbind(length(cmc_eqtl$snp_id), length(cmc_isoqtl$snp_id)))
    write.table(states, file = paste0(cmc_outpath, "/CMC_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to CMC ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Get only background SNPs overlapping with the CMC gene cis-eQTL dataset (load data first)
    cmc <- readRDS(paste0(inpath, "/CMC_adjustedSVA_eQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, cmc)
    
    # Run annotation for CMC (Fromer et al.) DPC gene cis-eQTLs
    results <- RunAnnotationForAllEQTLs(cmc_eqtl, "CMC_eQTL_SNP", results, cmc_outpath, perm, background)
    rm(background, cmc) # delete unnecessary data
    
    # Get only background SNPs overlapping with the CMC isoform-QTL (isoQTL) dataset (load data first)
    cmc <- readRDS(paste0(inpath, "/CMC_adjustedSVA_isoQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, cmc)
    
    # Run annotation for CMC (Fromer et al.) DPC isoQTL
    results <- RunAnnotationForAllEQTLs(cmc_isoqtl, "CMC_isoQTL_SNP", results, cmc_outpath, perm, background)

    # Write results into a file and remove unnecessary files
    write.table(results, file = paste0(cmc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_CMC.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(cmc_outpath, states, cmc_eqtl, cmc_isoqtl, results, perm, background, cmc)
    
  } else if (what == "meQTL") {
    # Generate output path
    meqtl_outpath <- paste0(outpath, "/05_meQTL")
    system(paste0("mkdir -p ", meqtl_outpath))
    
    # Load meQTL
    cortex_meqtl <- readRDS(paste0(inpath, "/Cortex_meQTL_SNPs_p0.05.rds"))

    # write number of SNPs in meQTL into a file
    states <- cbind("Cortex meQTL", length(cortex_meqtl$snp_id))
    write.table(states, file = paste0(meqtl_outpath, "/Cortex_meQTL_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")

    ### Annotation to meQTL ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Get only background SNPs overlapping with the cortex meQTL dataset (load data first)
    meqtl <- readRDS(paste0(inpath, "/Cortex_meQTL_SNPs.rds"))
    background <- subsetByOverlaps(background_gr, meqtl)
    
    # Run annotation for meQTL (Fromer et al.) DPC gene cis-eQTLs
    results <- RunAnnotationForAllEQTLs(cortex_meqtl, "Cortex_meQTL_SNP", results, meqtl_outpath, perm, background)

    # Write results into a file
    write.table(results, file = paste0(meqtl_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_cortex_meqtl.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(meqtl_outpath, cortex_meqtl, results, perm, background, meqtl)
    
  } else if (what == "GWAS") {
    # Generate output path
    gwas_outpath <- paste0(outpath, "/06_GWAS")
    system(paste0("mkdir -p ", gwas_outpath))
    
    # generate empty data frames to fill
    states_p0.05 <- data.frame()
    states_p5e8 <- data.frame()
    
    # get all GWAS RDS input file paths in the same order as the gwas and gwas_gw list
    gwas_names <- c("SCZ1-2012", "SCZ2-2014", "SCZ3-2020", "BIP-2019", "MDD1-2012", "MDD2-2018", "MDD3-2019", "CDG1-2013", "CDG2-2019", "ADHD-2019", "ASD-2019", "EA-2018", "BMI-2015_ENGAGE", "BMI-2015_GIANT", "BMI-2018", "BMI-2019", "AD-2019", "ALCDEP-2018", "AN-2019", "AUDIT-2019", "CAD-2017", "Cannabis-2018", "ED-2017", "HE-2014", "IN-2018", "OCD-2017", "PTSD-2019", "T2D-2017", "TS-2018")
    
    # load GWAS files and save them into a GRangesList
    gwas <- GRangesList()
    gwas_gw <- GRangesList()
    gwas_files <- list()  # save all unfiltered GWAS summary statistic file names with path
    for (i in seq_along(gwas_names)) {
      # load public data: GWAS summary statistic file
      sumstat <- readRDS(paste0(inpath, "/GWAS_sumstats_", gwas_names[[i]], "_p0.05.rds"))
      # Save GRanges objects of all GWAS (nominal p < 0.05) in the order they should be running
      gwas <- append(gwas, GRangesList(gwas = sumstat), after = length(gwas))
      # Save GRanges objects of all GWAS (genome-wide filtered: p < 5e-8) in the order they should be running
      sumstat <- sumstat[sumstat$p_value <= 5e-8, ]
      gwas_gw <- append(gwas_gw, GRangesList(gwas = sumstat), after = length(gwas_gw))
      # Add correct GWAS name to the list (did not work directly)
      if (i == 1) {
        names(gwas) <- gwas_names[[i]]
        names(gwas_gw) <- gwas_names[[i]]
      } else {
        names(gwas) <- c(names(gwas[1:(i-1)]), gwas_names[[i]])
        names(gwas_gw) <- c(names(gwas_gw[1:(i-1)]), gwas_names[[i]])
      }
      rm(sumstat)
      # get all GWAS RDS input file paths in the same order as the gwas and gwas_gw list
      gwas_files <- append(gwas_files, paste0(inpath, "/GWAS_sumstats_", gwas_names[[i]], ".rds"))
    }
    
    # save number of SNPs in each GWAS with nominal p (p < 0.05) into a file
    states <- as.data.frame(elementNROWS(gwas))
    write.table(states, file = paste0(gwas_outpath, "/GWAS_p0.05_sizes.tsv"), row.names = T, col.names = F, quote = F, sep = "\t")
    rm(states)
    
    # save number of SNPs in each GWAS with genome-wide p (p < 5e-8) into a file
    states <- as.data.frame(elementNROWS(gwas_gw))
    write.table(states, file = paste0(gwas_outpath, "/GWAS_p5e-8_sizes.tsv"), row.names = T, col.names = F, quote = F, sep = "\t")
    rm(states)
    
    ### Get overlap of all GWAS (p < 0.05) to CMC, GTEx wb, GTEx brain and GTEx BA9
    results <- data.frame()
    for (i in seq_along(gwas)) {
      gwas_selected <- gwas[[i]]
      print(names(gwas)[[i]])
      results[names(gwas)[[i]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, gwas_selected))
      results[names(gwas)[[i]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_eQTL_ba9"] <- sum(overlapsAny(gtex_ba9_eqtl, gwas_selected))
      results[names(gwas)[[i]], "GTEx_sQTL_ba9"] <- sum(overlapsAny(gtex_ba9_sqtl, gwas_selected))
    }
    ### write results into a file
    write.table(results, file = paste0(gwas_outpath, "/CMC_GTEx_SNP_overlap_to_all_GWAS_p0.05.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(results, gwas_selected)
    
    ### Get overlap of all GWAS (p < 5e-08) to CMC, GTEx wb, brain, ba9
    results <- data.frame()
    for (i in seq_along(gwas_gw)) {
      gwas_selected <- gwas_gw[[i]]
      print(names(gwas_gw)[[i]])
      results[names(gwas_gw)[[i]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_eQTL_ba9"] <- sum(overlapsAny(gtex_ba9_eqtl, gwas_selected))
      results[names(gwas_gw)[[i]], "GTEx_sQTL_ba9"] <- sum(overlapsAny(gtex_ba9_sqtl, gwas_selected))
    }
    ### write results into a file
    write.table(results, file = paste0(gwas_outpath, "/CMC_GTEx_SNP_overlap_to_all_GWAS_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(results, gwas_selected)
    
    ### Annotation to genome-wide GWAS hits ###
    # First, set parameters
    perm <- 10000
    # Annotation to GWAS summary statistics (parallelized)
    cl <- makeCluster(10)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(gwas_gw), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(names(gwas_gw)[i], "_SNP_overlap")
                           pval <- paste0(names(gwas_gw)[i], "_SNP_pVal")
                           or <- paste0(names(gwas_gw)[i], "_SNP_or")
                           gwas_selected <- gwas_gw[[i]]
                           gwas_full <- readRDS(gwas_files[[i]])
                           background <- subsetByOverlaps(background_gr, gwas_full)
                           rm(gwas_full)
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, gwas_selected, perm, background)
                           write.table(results, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas_gw)[i],"_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(results2)
    
    ### Annotation to nominal GWAS hits ###
    # First, set parameters
    perm <- 10000
    # Aannotation to GWAS summary statistics (parallelized)
    cl <- makeCluster(10)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(gwas), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(names(gwas)[i], "_SNP_overlap")
                           pval <- paste0(names(gwas)[i], "_SNP_pVal")
                           or <- paste0(names(gwas)[i], "_SNP_or")
                           gwas_selected <- gwas[[i]]
                           gwas_full <- readRDS(gwas_files[[i]])
                           background <- subsetByOverlaps(background_gr, gwas_full)
                           rm(gwas_full)
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, gwas_selected, perm, background)
                           write.table(results, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas)[i],"_p0.05.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p0.05.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
  } else {
    print("Parameter not defined!")
  }
  
  # Test frequency of eQTLs on chromosomes
  #eqtls_chr <- as.data.frame(table(seqnames(eqtls_transcript_cases_gr)))
  #eqtls_hmm <- as.data.frame(table(seqnames(hmm_blood)))
  #chromosomes <- merge(eqtls_chr, eqtls_hmm, by = "Var1")
  #print(wilcox.test(chromosomes$Freq.x, chromosomes$Freq.y, alternative = "two.sided"))
  print(sessionInfo())
}
