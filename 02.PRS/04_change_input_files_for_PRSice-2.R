library(data.table)

# read in IID and corresponding phenotype information (1 = controls, 2 = cases -> SCZ, BIP and MDD)
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

### ADHD-2017 ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/ADHD-2017/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_adhd_filtered.fam")

### ASD ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/ASD/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_asd_filtered.fam")

### BIP-2018 ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/BIP-2018/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_filtered.fam")

### BIP+SCZvsCONT ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/BIP+SCZvsCONT/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_scz_vs_cont_filtered.fam")

### BIPvsSCZ ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/BIPvsSCZ/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_vs_scz_filtered.fam")

### BMI ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/BMI/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bmi_filtered.fam")

### Cross-Disorder ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/Cross_Disorder/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_cd_filtered.fam")

### EA #####
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/EA/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_ea_filtered.fam")

### MDD2 ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/MDD2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_mdd2_filtered.fam")

### SCZ2 ###
ChangePheno("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2/SCZ2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_scz2_filtered.fam")
