# Script to generate a list of samples, we want to include in the analysis (removing all other samples)

# read pheno data
pData <- read.table("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2019-04-09_sample_sheet_KW.txt", sep="\t", header=TRUE, stringsAsFactors = T)
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
cat(sapply(sort(pData$IID), toString), file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/eqtl_samples.txt", sep="\n")
