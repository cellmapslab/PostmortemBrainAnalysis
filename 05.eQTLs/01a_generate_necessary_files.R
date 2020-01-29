###############################################################
### Generate input files for the preprcoessing of genotypes ###
###############################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args input1=<path/sample-file> input2=<path/sample-sheet> output1=<path/output-file1> output2=<path/output-file2>" 01a_generate_input_files.R

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to change the genotype probabilities (.gen.gz) file into the matrixEQTL input file format with only one probability value instead of three.\n\nYou can invoke the script as follows:\nR CMD BATCH --slave '--args geno=\"/home/.../genotype-file.gen.gz\" sample=\"/home/.../file.sample\" ids=\"/home/.../ids.txt\" output=\"/home/.../out-file.gen\"' 02b_prepare_geno.R\nCommandline parameters:\n\ngeno\t\tpath and filename of the merged genotype probabilites (ending in gen.gz)\nsample\t\tpath and filename of the sample file corresponding to the genotype file\nids\t\tpath and filename, which includes the ids of overlapping samples from the gene expression and the genotype dataset (here: 180)\noutput\t\tpath and filename of the output (input file for matrixEQTL)\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  
  ### 1. Copy paste the exclude sample file from 01.GWAS to 05.eQTLs 
  system(paste0("cp ", input1, " ", outpath))
  
  ### 2. Generate file including only samples for which also expression data is available
  # read pheno data
  pData <- read.table(input2, sep = "\t", header = TRUE, stringsAsFactors = T)
  length(pData$IID) # 209
  
  # remove older replicates
  replicates <- pData[(duplicated(pData$MHRIID) | duplicated(pData$MHRIID, fromLast = TRUE)) & pData$Batch == "A", "Filename"]
  pData <- pData[-(which(pData$Filename %in% replicates)), ] # remove the replicates from phenotype data
  length(pData$IID) # 204
  
  # filter out all samples with NA in IID
  pData <- pData[!is.na(pData$IID), ]
  length(pData$IID) # 202
  
  # filter out all samples with NA in C1
  pData <- pData[!is.na(pData$C1), ]
  length(pData$IID) # 171
  
  # filter out all samples with NA in Type_of_Death
  pData <- pData[!is.na(pData$Type_of_Death), ]
  length(pData$IID) # 169
  
  # save all important sample IIDs for filtering of the genotype data
  cat(sapply(sort(pData$IID), toString), file = output, sep = "\n")
}
