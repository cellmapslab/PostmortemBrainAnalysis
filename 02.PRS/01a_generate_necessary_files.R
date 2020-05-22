###############################################################
### Generate input files for the preprcoessing of genotypes ###
###############################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args input1=<path/sample-file> input2=<path/sample-sheet> output1=<path/output-file1> output2=<path/output-file2> output3=<path/output-file3> output4=<path/output-file4>" 01a_generate_necessary_files.R
# or using Rscript: Rscript --slave 01a_generate_necessary_files.R input1=\"/home/.../input-file1.txt\" input2=\"/home/.../input2.txt\" output1=\"/home/.../out-file1.txt\" output2=\"/home/.../out-file2.txt\" output3=\"/home/.../out-file3.txt\" output4=\"/home/.../out-file4.txt

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate files necessary for the genotype preprocesing.\n\nYou can invoke the script as follows:\nRscript --slave 01a_generate_necessary_files.R input1=\"/home/.../input-file1.txt\" input2=\"/home/.../input2.txt\" output1=\"/home/.../out-file1.txt\" output2=\"/home/.../out-file2.txt\" output3=\"/home/.../out-file3.txt\" output4=\"/home/.../out-file4.txt\" output5=\"/home/.../out-file5.txt\"\nCommandline parameters:\n\ninput1\t\tpath and filename of the merged genotype probabilites (ending in gen.gz)\ninput2\t\tpath and filename of the sample file corresponding to the genotype file\noutput1\t\tpath and filename of the output file including all sample IDs to exclude from further analysis as they are coming from other (unknown) cohorts\noutput2\t\tpath and filename of the output file including the 169 samples we want to include in all analysis (samples including genotype AND expression data)\noutput3\t\tpath and filename of the output file including all covariates\noutput4\t\tpath and filename of the output file including the phenotypes\noutput5\t\tpath and filename of the output file including the covariates for LDpred (No FID, all covariates as continuous/numerical data, columns separated by whitespace instead of tab)\n\n")
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

  
  ### 3. Generate a file including all covariates for PRSice-2 and lassosum
  # get all important columns for the covariate file
  covars <- pData[, c("FID", "IID", "Age", "Sex", "Suicide", "Type_of_Death", "C1", "C2", "C3", "C4")]
  covars[, ] <- lapply(covars[, ], as.character) # change all factors to character
  write.table(covars, file = output3, sep = "\t", quote = FALSE, row.names = FALSE)


  ### 4. Save only FID, IID and DX (Diagnosis)
  pheno <- pData[, c("FID", "IID", "Dx")]
  pheno[, ] <- lapply(pheno[, ], as.character) # change all factors to character

  # change Diagnosis (Dx) from characters to numbers with all NON controls as one case (CTL = , SCZ/BPD/MDD = 2) 
  pheno[pheno$Dx == "CTL", "Dx"] <- "1"
  pheno[pheno$Dx %in% c("SCZ", "BPD", "MDD"), "Dx"] <- "2"

  # write phenotypes into file
  write.table(pheno, file = output4, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


  ### 5. Generate a file including all covariates for LDpred 
  # categorical covariates need to be changed to continuous/numerical ones, do not include FID column and whitespace between columns instead of tabs
  # get all important columns for the covariate file
  covars <- pData[, c("IID", "Age", "Sex", "Suicide", "Type_of_Death", "C1", "C2", "C3", "C4")]
  covars[, ] <- lapply(covars[, ], as.character) # change all factors to character
  
  # change sex from characters to numbers (M/Male = 1, F/Female = 2)
  covars[covars$Sex == "M", "Sex"] <- "1"
  covars[covars$Sex == "F", "Sex"] <- "2"
  
  # change suicide from characters to numbers (N/No = 1, Y/Yes = 2)
  covars[covars$Suicide == "N", "Suicide"] <- "1"
  covars[covars$Suicide == "Y", "Suicide"] <- "2"
  
  # change Type of Death from characters to numbers (Natural = 1, Non-violent = 2, Violent = 3) 
  covars[covars$Type_of_Death == "Natural", "Type_of_Death"] <- "1"
  covars[covars$Type_of_Death == "Non-violent", "Type_of_Death"] <- "2"
  covars[covars$Type_of_Death == "Violent", "Type_of_Death"] <- "3"
   
  # save all important sample IIDs for filtering of the genotype data
  write.table(covars, file = output5, sep = " ", quote = FALSE, row.names = FALSE)

  print(sessionInfo())
}
