#!/bin/bash

tools_path="/storage/groups/ccm01/tools"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS"

###########################
### Run GWAS with PLINK ###
###########################
# 1. With Covariates
# --logistic uses either linear (without control data) or logistic regression (with control data)
$tools_path/PLINK/plink --bfile ${out_genotypes}/01_preprocessed_genotypes/1.4_GEX_filtered/pmb_info06_all_chr_maf0.01_merged_GEX_filtered --logistic hide-covar --freq --pheno /${out_genotypes}/pheno.txt --allow-no-sex --covar /${out_genotypes}/covar.txt --covar-name SEX,AGE,SUICIDE,TYPE_OF_DEATH,C1,C2,C3,C4 --out /${out_genotypes}/02_GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS 

# 2. Without Covariates
# --assoc does not allow any covariates (ignores the --covar parameter, if set)
$tools_path/PLINK/plink --bfile ${out_genotypes}/01_preprocessed_genotypes/1.4_GEX_filtered/pmb_info06_all_chr_maf0.01_merged_GEX_filtered --assoc --freq --pheno /${out_genotypes}/pheno.txt --allow-no-sex --covar /${out_genotypes}/covar.txt --covar-name SEX,AGE,SUICIDE,TYPE_OF_DEATH,C1,C2,C3,C4 --out /${out_genotypes}/02_GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS 

