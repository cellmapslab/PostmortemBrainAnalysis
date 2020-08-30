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
  # @param background   All background genes coming from the gene expression (here: 14,048 unqiue genes)
  # @return             Vector including the number of overlapping SNPs (1 = overlap), the resampling p-value (2 = p_value) and the resampling odds ratio (3 = or)
  AnnotatePublicData <- function(own, public, nperm, background){
    # get only background SNPs overlapping with ChromHMM
    background <- subsetByOverlaps(background, public)
    #print(paste0("Number of overlapping background SNPs: ", length(background@ranges)))  # print number of overlapping background SNPs
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
    return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling p-value
  }
  
  # function, to calculate the annotation for one public dataset to all eQTL datasets
  # @param input  Public dataset (e.g. CMC eQTLs, GWAS AD-2019)
  # @param name   Name of the public dataset
  # @param res    Data frame in which the results should be stored
  # @return       the results data frame including the newly calculated results
  RunAnnotationForAllEQTLs <- function(input, name, results, output, nperm, background){
    # get only background SNPs overlapping with ChromHMM
    background <- subsetByOverlaps(background, input)
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
  cmc_eqtl <- readRDS(paste0(inpath, "/CMC_adjustedSVA_eQTL_SNPs.rds"))
  cmc_isoqtl <- readRDS(paste0(inpath, "/CMC_adjustedSVA_isoQTL_SNPs.rds"))
  print("Done loading CMC")
  
  ############
  ### GTEx ###
  ############
  ### Load GTEx Whole Blood eQTLs and sQTLs
  gtex_wb_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
  gtex_wb_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))
  print("Done loading GTEx whole blood")
  ### GTEx Brain (all brain tissues) eQTLs and sQTLs
  gtex_brains_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
  gtex_brains_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))
  print("Done loading GTEx all brains")
  ### GTEx Brain (only BA9 = Broadmann area 9) eQTLs and sQTLs
  gtex_ba9_eqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_eQTL_SNPs.rds"))
  gtex_ba9_sqtl <- readRDS(paste0(inpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_sQTL_SNPs.rds"))
  print("Done loading GTEx BA9")
  
  # run annotation for the dataset (what = VEP/ChromHMM/CMC/GTEx/meQTL/GWAS)
  if (what == "test") {
    print("\nsuccessfull\n\n")
    
  } else if (what == "VEP") {
    # Generate output path
    vep_outpath <- paste0(outpath, "/01_VEP")
    system(paste0("mkdir -p ", vep_outpath))
    
    ### Load VEP ###
    vep <- readRDS(paste0(inpath, "/VEP_hg19_collapsed_MAF01.rds"))
   
    # Vep codes
    vep_counts <- table(vep$code)
    states <- names(vep_counts)
    print(states)
    
    # Save states into a file
    vep_counts <- as.data.frame(vep_counts)
    vep_counts$Var1 <- substr(vep_counts$Var1, 1, nchar(as.character(vep_counts$Var1)) - 2)
    write.table(vep_counts, file = paste0(vep_outpath, "/VEP_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
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
    cl <- makeCluster(17)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_overlap")
                           pval <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_pVal")
                           or <- paste0(strsplit(states[i], "_d$")[[1]], "_SNP_or")
                           vep_selected <- vep[vep$code == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, vep_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, vep_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, vep_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, vep_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, vep_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, vep_selected, perm, background_gr)
                           write.table(results, file = paste0(vep_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_", strsplit(states[i], "_d$")[[1]],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(vep_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_VEP_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "ChromHMM_Blood") {
    # Generate output path
    hmm_blood_outpath <- paste0(outpath, "/02_ChromHMM_Blood")
    system(paste0("mkdir -p ", hmm_blood_outpath))
    
    # Load ChromHMM Blood & T-cells
    hmm_blood <- readRDS(paste0(inpath, "/ChromHMM_Blood_and_T-cells.rds"))
    
    # Get ChromHMM states
    hmm_blood_counts <- table(hmm_blood$name)
    states <- names(hmm_blood_counts)
    print(states)
    
    # Save states into a file
    hmm_blood_counts <- as.data.frame(hmm_blood_counts)
    write.table(hmm_blood_counts, file = paste0(hmm_blood_outpath, "/ChromHMM_Blood_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
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
    cl <- makeCluster(15)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_blood_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_blood_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_blood_", states[i], "_SNP_or")
                           hmm_blood_selected <- hmm_blood[hmm_blood$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_blood_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_blood_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_blood_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_blood_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_blood_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_blood_selected, perm, background_gr)
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(hmm_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "ChromHMM_Brain") {
    # Generate output path
    hmm_brain_outpath <- paste0(outpath, "/02_ChromHMM_Brain")
    system(paste0("mkdir -p ", hmm_brain_outpath))
    
    # Load ChromHMM all Brains
    hmm_brain <- readRDS(paste0(inpath, "/ChromHMM_all_Brains.rds"))
    
    # Get ChromHMM states
    hmm_brain_counts <- table(hmm_brain$name)
    states <- names(table(hmm_brain_counts))
    print(states)
    
    # Save states into a file
    hmm_brain_counts <- as.data.frame(hmm_brain_counts)
    write.table(hmm_brain_counts, file = paste0(hmm_brain_outpath, "/ChromHMM_Brain_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
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
    cl <- makeCluster(15)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_brain_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_brain_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_brain_", states[i], "_SNP_or")
                           hmm_brain_selected <- hmm_brain[hmm_brain$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_brain_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_brain_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_brain_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_brain_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_brain_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_brain_selected, perm, background_gr)
                           write.table(results, file = paste0(hmm_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_", states[i],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(hmm_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain_all_states.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "ChromHMM_DLPFC") {
    # Generate output path
    hmm_dlpfc_outpath <- paste0(outpath, "/02_ChromHMM_DLPFC")
    system(paste0("mkdir -p ", hmm_dlpfc_outpath))
    
    # Load ChromHMM Brain DLPFC only
    hmm_dlpfc <- readRDS(paste0(inpath, "/ChromHMM_Brain_Dorsolateral_Prefrontal_Cortex.rds"))

    # Get ChromHMM states
    hmm_dlpfc_counts <- table(hmm_dlpfc$name)
    states <- names(hmm_dlpfc_counts)
    print(states)
    
    # Save states into a file
    hmm_dlpfc_counts <- as.data.frame(hmm_dlpfc_counts)
    write.table(hmm_dlpfc_counts, file = paste0(hmm_dlpfc_outpath, "/ChromHMM_DLPFC_category_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
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
    write.table(results, file = paste0(hmm_dlpfc_outpath, "/CMC_GTEx_SNP_overlap_to_ChromHMM_brain_DLPFC_only.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
    # Annotation to ChromHMM Brain DLPFC (parallelized)
    perm <- 10000
    cl <- makeCluster(15)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(states), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0("ChromHMM_brain_DLPFC_", states[i], "_SNP_overlap")
                           pval <- paste0("ChromHMM_brain_DLPFC_", states[i], "_SNP_pVal")
                           or <- paste0("ChromHMM_brain_DLPFC_", states[i], "_SNP_or")
                           hmm_dlpfc_selected <- hmm_dlpfc[hmm_dlpfc$name == states[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, hmm_dlpfc_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, hmm_dlpfc_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, hmm_dlpfc_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, hmm_dlpfc_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, hmm_dlpfc_selected, perm, background_gr)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, hmm_dlpfc_selected, perm, background_gr)
                           write.table(results, file = paste0(hmm_dlpfc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_", states[i],".tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(hmm_dlpfc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_brain_DLPFC_only.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "GTEx_Blood") {
    # Generate output path
    gtex_blood_outpath <- paste0(outpath, "/03_GTEx_Blood")
    system(paste0("mkdir -p ", gtex_blood_outpath))
    
    # write number of SNPs in GTEx Blood
    states <- cbind(rbind("GTEx Blood eQTL", "GTEx Blood sQTL"), rbind(length(gtex_wb_eqtl$snp_id), length(gtex_wb_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_blood_outpath, "/GTEx_Blood_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Run annotation for GTEx whole blood gene cis-eQTLs and isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_wb_eqtl, "GTEx_whole_blood_eQTL_SNP", results, gtex_blood_outpath, perm, background_gr)
    results <- RunAnnotationForAllEQTLs(gtex_wb_sqtl, "GTEx_whole_blood_sQTL_SNP", results, gtex_blood_outpath, perm, background_gr)

    # Write results into a file
    write.table(results, file = paste0(gtex_blood_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_whole_blood.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "GTEx_Brain") {
    # Generate output path
    gtex_brain_outpath <- paste0(outpath, "/03_GTEx_Brain")
    system(paste0("mkdir -p ", gtex_brain_outpath))
    
    # write number of SNPs in GTEx Brain
    states <- cbind(rbind("GTEx Brain eQTL", "GTEx Brain sQTL"), rbind(length(gtex_brains_eqtl$snp_id), length(gtex_brains_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_brain_outpath, "/GTEx_Brain_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Run annotation for GTEx whole blood gene cis-eQTLs and isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_brains_eqtl, "GTEx_all_brains_eQTL_SNP", results, gtex_brain_outpath, perm, background_gr)
    results <- RunAnnotationForAllEQTLs(gtex_brains_sqtl, "GTEx_all_brains_sQTL_SNP", results, gtex_brain_outpath, perm, background_gr)
    
    # Write results into a file
    write.table(results, file = paste0(gtex_brain_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_all_brains.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "GTEx_BA9") {
    # Generate output path
    gtex_ba9_outpath <- paste0(outpath, "/03_GTEx_BA9")
    system(paste0("mkdir -p ", gtex_ba9_outpath))
    
    # write number of SNPs in GTEx BA9
    states <- cbind(rbind("GTEx Brain BA9 eQTL", "GTEx Brain BA9 sQTL"), rbind(length(gtex_ba9_eqtl$snp_id), length(gtex_ba9_sqtl$snp_id)))
    write.table(states, file = paste0(gtex_ba9_outpath, "/GTEx_Brain_BA9_only_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to GTEx ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Run annotation for GTEx whole blood gene cis-eQTLs and isoform cis-eQTL (sQTL)
    results <- RunAnnotationForAllEQTLs(gtex_ba9_eqtl, "GTEx_brain_BA9_eQTL_SNP", results, gtex_ba9_outpath, perm, background_gr)
    results <- RunAnnotationForAllEQTLs(gtex_ba9_sqtl, "GTEx_brain_BA9_sQTL_SNP", results, gtex_ba9_outpath, perm, background_gr)
    
    # Write results into a file
    write.table(results, file = paste0(gtex_ba9_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_GTEx_brain_BA9.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "CMC") {
    # Generate output path
    cmc_outpath <- paste0(outpath, "/04_CMC")
    system(paste0("mkdir -p ", cmc_outpath))
    
    # write number of SNPs in CMC
    states <- cbind(rbind("CMC eQTL", "CMC isoQTL"), rbind(length(cmc_eqtl$snp_id), length(cmc_isoqtl$snp_id)))
    write.table(states, file = paste0(cmc_outpath, "/CMC_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    ### Annotation to CMC ###
    # Set parameters
    results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
    perm <- 10000
    
    # Run annotation for CMC (Fromer et al.) DPC gene cis-eQTLs
    results <- RunAnnotationForAllEQTLs(cmc_eqtl, "CMC_eQTL_SNP", results, cmc_outpath, perm, background_gr)
    
    # Run annotation for CMC (Fromer et al.) DPC isoform cis-eQTLs
    results <- RunAnnotationForAllEQTLs(cmc_isoqtl, "CMC_isoQTL_SNP", results, cmc_outpath, per, background_gr)
    
    # Write results into a file
    write.table(results, file = paste0(cmc_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_CMC.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    
  } else if (what == "meQTL") {
    # Generate output path
    meqtl_outpath <- paste0(outpath, "/05_meQTL")
    system(paste0("mkdir -p ", meqtl_outpath))
    
    # Load meQTL
    cortex_meqtl_gr <- readRDS(paste0(inpath, "/Cortex_meQTL_SNPs.rds"))

    # write number of SNPs in meQTL
    states <- cbind("Cortex meQTL", length(cortex_meqtl_gr$snp_id))
    write.table(states, file = paste0(meqtl_outpath, "/Cortex_meQTL_sizes.tsv"), row.names = F, col.names = F, quote = F, sep = "\t")
    
    # Annotate meQTLs
    results <- data.frame()
    perm <- 10000
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene"] <- AnnotatePublicData(eqtls_gene_gr, cortex_meqtl_gr, perm, background_gr)
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases_gr, cortex_meqtl_gr, perm, background_gr)
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript"] <- AnnotatePublicData(eqtls_transcript_gr, cortex_meqtl_gr, perm, background_gr)
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases_gr, cortex_meqtl_gr, perm, background_gr)
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon"] <- AnnotatePublicData(eqtls_exon_gr, cortex_meqtl_gr, perm, background_gr)
    results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases_gr, cortex_meqtl_gr, perm, background_gr)
    
    # Write results into a file
    write.table(results, file = paste0(meqtl_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_meQTLs.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")

  } else if (what == "GWAS") {
    # Generate output path
    gwas_outpath <- paste0(outpath, "/06_GWAS")
    system(paste0("mkdir -p ", gwas_outpath))
    
    ### Load public data: GWAS summary statistics ###
    # datasets are already filtered for p-value <= 0.05
    ad <- readRDS(paste0(inpath, "/GWAS_sumstats_AD-2019_p0.05.rds"))  # GWAS Alzheimer's Disease (AD-2019)
    adhd <- readRDS(paste0(inpath, "/GWAS_sumstats_ADHD-2019_p0.05.rds"))  # GWAS Attention-Deficit Hyperactivity Disorder (ADHD-2017)
    alcdep <- readRDS(paste0(inpath, "/GWAS_sumstats_ALCDEP-2018_p0.05.rds"))  # GWAS Alcohol Dependence (ALCDEP-2018) <- similar to AUDIT-2018
    an <- readRDS(paste0(inpath, "/GWAS_sumstats_AN-2019_p0.05.rds"))  # GWAS Anorexia Nervosa (AN-2019) <- similar to ED-2017
    asd <- readRDS(paste0(inpath, "/GWAS_sumstats_ASD-2019_p0.05.rds"))  # GWAS Autism Spectrum Disorder (ASD-2017)
    audit <- readRDS(paste0(inpath, "/GWAS_sumstats_AUDIT-2019_p0.05.rds"))  # GWAS Alcohol Use Disorders Identification Test (AUDIT-2018) <- similar to ALCDEP-2018
    bip <- readRDS(paste0(inpath, "/GWAS_sumstats_BIP-2019_p0.05.rds"))  # GWAS Bipolar Disorder (BIP-2018)
    bmi5e <- readRDS(paste0(inpath, "/GWAS_sumstats_BMI-2015_ENGAGE_p0.05.rds"))  # GWAS Body Mass Index (BMI-2015_ENGAGE)
    bmi5g <- readRDS(paste0(inpath, "/GWAS_sumstats_BMI-2015_GIANT_p0.05.rds"))  # GWAS Body Mass Index (BMI-2015_GIANT)
    bmi8 <- readRDS(paste0(inpath, "/GWAS_sumstats_BMI-2018_p0.05.rds"))  # GWAS Body Mass Index (BMI-2018)
    bmi9 <- readRDS(paste0(inpath, "/GWAS_sumstats_BMI-2019_p0.05.rds"))  # GWAS Body Mass Index (BMI-2019)
    cad <- readRDS(paste0(inpath, "/GWAS_sumstats_CAD-2017_p0.05.rds"))  # GWAS Coronary Artery Disease (CAD-2017)
    cannabis <- readRDS(paste0(inpath, "/GWAS_sumstats_Cannabis-2018_p0.05.rds"))  # GWAS Cannabis Use (Cannabis-2018)
    cdg1 <- readRDS(paste0(inpath, "/GWAS_sumstats_CDG1-2013_p0.05.rds"))  # GWAS Cross-Disorder Group (CDG1-2013)
    cdg2 <- readRDS(paste0(inpath, "/GWAS_sumstats_CDG2-2019_p0.05.rds"))  # GWAS Cross-Disorder Group (CDG2-2019)
    ea <- readRDS(paste0(inpath, "/GWAS_sumstats_EA-2018_p0.05.rds"))  # GWAS Educational Attainment (EA-2018) <- similar to Intelligence-2018
    ed <- readRDS(paste0(inpath, "/GWAS_sumstats_ED-2017_p0.05.rds"))  # GWAS Eating Disorder (ED-2017) <- similar to AN-2019
    height <- readRDS(paste0(inpath, "/GWAS_sumstats_HE-2014_p0.05.rds"))  # GWAS Height (Height-2014)
    intelligence <- readRDS(paste0(inpath, "/GWAS_sumstats_IN-2018_p0.05.rds"))  # GWAS Intelligence (Intelligence-2018) <- similar to EA-2018
    mdd1 <- readRDS(paste0(inpath, "/GWAS_sumstats_MDD1-2012_p0.05.rds"))  # GWAS Major Depression Disorder  (MDD-2012)
    mdd2 <- readRDS(paste0(inpath, "/GWAS_sumstats_MDD2-2018_p0.05.rds"))  # GWAS Major Depression Disorder  (MDD-2018)
    mdd3 <- readRDS(paste0(inpath, "/GWAS_sumstats_MDD3-2019_p0.05.rds"))  # GWAS Major Depression Disorder (MDD-2019)
    ocd <- readRDS(paste0(inpath, "/GWAS_sumstats_OCD-2017_p0.05.rds"))  # GWAS Obsessive Compulsive Disorder (OCD-2017)
    ptsd <- readRDS(paste0(inpath, "/GWAS_sumstats_PTSD-2019_p0.05.rds"))  # GWAS Posttraumatic Stress Disorder (PTSD-2019)
    scz2 <- readRDS(paste0(inpath, "/GWAS_sumstats_SCZ2-2014_p0.05.rds"))  # GWAS Schizophrenia (SCZ2-2014)
    t2d <- readRDS(paste0(inpath, "/GWAS_sumstats_T2D-2017_p0.05.rds"))  # GWAS Type 2 Diabetes (T2D-2017)
    ts <- readRDS(paste0(inpath, "/GWAS_sumstats_TS-2018_p0.05.rds"))  # GWAS Tourette Syndrom (TS-2018)
    
    # filter GWAS datasets for genome-wide (gw) hits (p-value <= 5e-8)  
    ad_gw <- ad[ad$p_value <= 5e-8, ]
    adhd_gw <- adhd[adhd$p_value <= 5e-8, ]
    alcdep_gw <- alcdep[alcdep$p_value <= 5e-8, ]
    an_gw <- an[an$p_value <= 5e-8, ]
    asd_gw <- asd[asd$p_value <= 5e-8, ]
    audit_gw <- audit[audit$p_value <= 5e-8, ]
    bip_gw <- bip[bip$p_value <= 5e-8, ]
    bmi5e_gw <- bmi5e[bmi5e$p_value <= 5e-8, ]
    bmi5g_gw <- bmi5g[bmi5g$p_value <= 5e-8, ]
    bmi8_gw <- bmi8[bmi8$p_value <= 5e-8, ]
    bmi9_gw <- bmi9[bmi9$p_value <= 5e-8, ]
    cad_gw <- cad[cad$p_value <= 5e-8, ]
    cannabis_gw <- cannabis[cannabis$p_value <= 5e-8, ]
    cdg1_gw <- cdg1[cdg1$p_value <= 5e-8, ]
    cdg2_gw <- cdg2[cdg2$p_value <= 5e-8, ]
    ea_gw <- ea[ea$p_value <= 5e-8, ]
    ed_gw <- ed[ed$p_value <= 5e-8, ]
    height_gw <- height[height$p_value <= 5e-8, ]
    intelligence_gw <- intelligence[intelligence$p_value <= 5e-8, ]
    mdd1_gw <- mdd1[mdd1$p_value <= 5e-8, ]
    mdd2_gw <- mdd2[mdd2$p_value <= 5e-8, ]
    mdd3_gw <- mdd3[mdd3$p_value <= 5e-8, ]
    ocd_gw <- ocd[ocd$p_value <= 5e-8, ]
    ptsd_gw <- ptsd[ptsd$p_value <= 5e-8, ]
    scz2_gw <- scz2[scz2$p_value <= 5e-8, ]
    t2d_gw <- t2d[t2d$p_value <= 5e-8, ]
    ts_gw <- ts[ts$p_value <= 5e-8, ]
    
    # Save names of all GWAS (nominal p < 0.05) in the order they should be running
    gwas <- GRangesList("GWAS_SCZ-2014" = scz2, "GWAS_BIP-2018" = bip, "GWAS_MDD1-2012" = mdd1, "GWAS_MDD2-2018" = mdd2, "GWAS_MDD3-2019" = mdd3, "GWAS_CDG1-2013" = cdg1, "GWAS_CDG2-2019" = cdg2, "GWAS_ADHD-2017" = adhd, "GWAS_ASD-2017" = asd, "GWAS_EA-2018" = ea, "GWAS_BMI-2015_ENGAGE" = bmi5e, "GWAS_BMI-2015_GIANT" = bmi5g, "GWAS_BMI-2018" = bmi8, "GWAS_BMI-2019" = bmi9, "GWAS_AD-2019" = ad, "GWAS_ALCDEP-2018" = alcdep, "GWAS_AN-2019" = an, "GWAS_AUDIT-2018" = audit, "GWAS_CAD_2017" = cad, "GWAS_Cannabis-2018" = cannabis, "GWAS_ED-2017" = ed, "GWAS_Height-2014" = height, "GWAS_Intelligence-2018" = intelligence, "GWAS_OCD-2017" = ocd, "GWAS_PTSD-2019" = ptsd, "GWAS_T2D-2017" = t2d, "GWAS_TS-2018" = ts)
    rm(ad, adhd, alcdep, an, asd, audit, bip, bmi5e, bmi5g, bmi8, bmi9, cad, cannabis, cdg1, cdg2, ea, ed, height, intelligence, mdd1, mdd2, mdd3, ocd, ptsd, scz2, t2d, ts)
    # save number of SNPs in each GWAS
    states <- as.data.frame(elementNROWS(gwas))
    write.table(states, file = paste0(gwas_outpath, "/GWAS_p0.05_sizes.tsv"), row.names = T, col.names = F, quote = F, sep = "\t")
    
    # Save names of all GWAS (genome-wide filtered: p < 5e-8) in the order they should be running
    gwas_gw <- GRangesList("GWAS_SCZ-2014" = scz2_gw, "GWAS_BIP-2018" = bip_gw, "GWAS_MDD1-2012" = mdd1_gw, "GWAS_MDD2-2018" = mdd2_gw, "GWAS_MDD3-2019" = mdd3_gw, "GWAS_CDG1-2013" = cdg1_gw, "GWAS_CDG2-2019" = cdg2_gw, "GWAS_ADHD-2017" = adhd_gw, "GWAS_ASD-2017" = asd_gw, "GWAS_EA-2018" = ea_gw, "GWAS_BMI-2015_ENGAGE" = bmi5e_gw, "GWAS_BMI-2015_GIANT" = bmi5g_gw, "GWAS_BMI-2018" = bmi8_gw, "GWAS_BMI-2019" = bmi9_gw, "GWAS_AD-2019" = ad_gw, "GWAS_ALCDEP-2018" = alcdep_gw, "GWAS_AN-2019" = an_gw, "GWAS_AUDIT-2018" = audit_gw, "GWAS_CAD_2017" = cad_gw, "GWAS_Cannabis-2018" = cannabis_gw, "GWAS_ED-2017" = ed_gw, "GWAS_Height-2014" = height_gw, "GWAS_Intelligence-2018" = intelligence_gw, "GWAS_OCD-2017" = ocd_gw, "GWAS_PTSD-2019" = ptsd_gw, "GWAS_T2D-2017" = t2d_gw, "GWAS_TS-2018" = ts_gw)
    rm(ad_gw, adhd_gw, alcdep_gw, an_gw, asd_gw, audit_gw, bip_gw, bmi5e_gw, bmi5g_gw, bmi8_gw, bmi9_gw, cad_gw, cannabis_gw, cdg1_gw, cdg2_gw, ea_gw, ed_gw, height_gw, intelligence_gw, mdd1_gw, mdd2_gw, mdd3_gw, ocd_gw, ptsd_gw, scz2_gw, t2d_gw, ts_gw)
    # save number of SNPs in each GWAS
    states <- as.data.frame(elementNROWS(gwas_gw))
    write.table(states, file = paste0(gwas_outpath, "/GWAS_p5e-8_sizes.tsv"), row.names = T, col.names = F, quote = F, sep = "\t")
    
    ### Get overlap of all GWAS (p < 0.05) to CMC, GTEx wb, brain, ba9
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
    cl <- makeCluster(12)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(gwas_gw), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(names(gwas_gw)[i], "_SNP_overlap")
                           pval <- paste0(names(gwas_gw)[i], "_SNP_pVal")
                           or <- paste0(names(gwas_gw)[i], "_SNP_or")
                           gwas_selected <- gwas_gw[[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases, gwas_selected, perm, background)
                           write.table(results, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas_gw)[i],"_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
                           results
                         }
    stopCluster(cl)
    # Write results into a file
    write.table(results2, file = paste0(gwas_outpath, "/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p5e-8.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
    rm(results, results2)
    
    ### Annotation to nominal GWAS hits ###
    # First, set parameters
    perm <- 10000
    # Aannotation to GWAS summary statistics (parallelized)
    cl <- makeCluster(12)
    registerDoSNOW(cl)
    clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above
    results2 <- foreach (i = seq_along(gwas), .combine = rbind,
                         .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                           results <- data.frame()
                           overlap <- paste0(names(gwas)[i], "_SNP_overlap")
                           pval <- paste0(names(gwas)[i], "_SNP_pVal")
                           or <- paste0(names(gwas)[i], "_SNP_or")
                           gwas_selected <- gwas[[i]]
                           results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon, gwas_selected, perm, background)
                           results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases, gwas_selected, perm, background)
                           write.table(results, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas)[i],"_p0.05.tsv"), col.names = NA, quote = F, dec = ".", sep = "\t")
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
