# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)

# function to change the diagnosis column of the phenotypes in .fam from -9 to 1 for controls and 2 for cases (SCZ, BIP and MDD) using the IIDs to diagnosis file
pdata <- read.table("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/IIDs_and_diagnosis.tsv", sep = "\t", header = FALSE, stringsAsFactors = T)

ChangePheno <- function(filename){
  # load postmortem .fam file  
  pmb <- read.table(file = filename, header = FALSE, stringsAsFactors = FALSE)
  # remove last column
  pmb$V6 <- NULL
  # add pheno for matching IIDs
  pmb <- merge(pmb, pdata, by = "V1")
  # now write the new fam file into a file
  write.table(pmb, filename, col.names = FALSE, row.names = FALSE, sep = " ")
}

# Change the phenotype .fam file of the PMB genotypes
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/01_preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final.fam")

print(sessionInfo())
