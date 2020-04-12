#!/bin/bash

tools_path="/storage/groups/ccm01/tools"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS"
outpath_geno="${outpath}/01_preprocessed_genotypes/1.3_samples_filtered"
# Generate the output directory (if not existing)
outpath_gwas="${outpath}/02_GWAS_results"  # quality control output path
mkdir -p $outpath_gwas

###########################
### Run GWAS with PLINK ###
###########################
# 1. With Covariates
# --logistic uses either linear (without control data) or logistic regression (with control data)
$tools_path/PLINK/plink --bfile ${outpath_geno}/pmb_info06_all_chr_maf0.01_merged_samples_filtered --logistic hide-covar --freq --pheno /${outpath}/pheno.txt --allow-no-sex --covar /${outpath}/covar.txt --covar-name SEX,AGE,SUICIDE,TYPE_OF_DEATH,C1,C2,C3,C4 --out /${outpath_gwas}/pmb_info06_all_chr_maf0.01_merged_samples_filtered_GWAS_with_covar

# 2. Without Covariates
# --assoc does not allow any covariates (ignores the --covar parameter, if set)
$tools_path/PLINK/plink --bfile ${outpath_geno}/pmb_info06_all_chr_maf0.01_merged_samples_filtered --assoc --freq --pheno /${outpath}/pheno.txt --allow-no-sex --covar /${outpath}/covar.txt --covar-name SEX,AGE,SUICIDE,TYPE_OF_DEATH,C1,C2,C3,C4 --out /${outpath_gwas}/pmb_info06_all_chr_maf0.01_merged_samples_filtered_GWAS_without_covar

