###############################################################
### Generate input files for the preprcoessing of genotypes ###
###############################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args input1=<path/sample-file> input2=<path/sample-sheet> output1=<path/output-file1> output2=<path/output-file2> 01a_generate_input_files.R
# or using Rscript: Rscript --slave 01a_generate_necessary_files.R input1=\"/home/.../input-file1.txt\" input2=\"/home/.../input2.txt\" output1=\"/home/.../out-file1.txt\" output2=\"/home/.../out-file2.txt\"

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate files necessary for the genotype preprocesing.\n\nYou can invoke the script as follows:\nRcript --slave 01a_generate_necessary_files.R input1=\"/home/.../input-file1.txt\" input2=\"/home/.../input2.txt\" output1=\"/home/.../out-file1.txt\" output2=\"/home/.../out-file2.txt\"\nCommandline parameters:\n\ninput1\t\tpath and filename of the merged genotype probabilites (ending in gen.gz)\ninput2\t\tpath and filename of the sample file corresponding to the genotype file\noutput1\t\tpath and filename of the output file including all sample IDs to exclude from further analysis as they are coming from other (unknown) cohorts\noutput2\t\tpath and filename of the output file including the 169 samples we want to include in all analysis (samples including genotype AND expression data)\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  
  ### 1. Generate file excluding the samples of other studies
  # read in sample file of the first chromosome
  samples <- read.table(input1, sep = " ", header = TRUE, stringsAsFactors = F)
  
  # get all samples including letters (as our IIDs are only numbers!!)
  exclude <- samples[grepl("[a-zA-Z]+", samples$ID_1), "ID_1"]
  
  # save data into a file
  write.table(exclude, file = output1, quote = FALSE, sep = " ", eol = " ", row.names = FALSE, col.names = F)
  
  
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
  
  # generate output table
  out <- as.data.frame(cbind(pData$FID, pData$IID))
  colnames(out) <- c("FID", "IID")
  
  # save all important sample IIDs for filtering of the genotype data
  write.table(out, file = output2, sep = "\t", quote = FALSE, row.names = FALSE)

  print(sessionInfo())
}

