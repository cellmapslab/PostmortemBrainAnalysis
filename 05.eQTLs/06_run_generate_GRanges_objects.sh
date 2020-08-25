#!/bin/bash
geno="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup"
maf_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.6_SNP_MAF_values"
frq="${maf_path}/pmb_info06_all_chr_maf0.05_merged_noDup_maf_information.frq"
eqtl_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/04_result_tables"
clump_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/05_clumping/03_result_tables"
cmc_path="/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL"
gwas_path="/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics"
gtex_path="/storage/groups/ccm01/datasets/2019_GTEx_v8"
vep="/storage/groups/ccm01/datasets/2019_Variant_Effect_Predictor/WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt"
chromhmm_path="/storage/groups/ccm01/datasets/2019_ChromHMM"
cortex="/storage/groups/ccm01/datasets/2019_cortex_meQTL/GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv"
out_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData_files"
# generate output folder
mkdir -p $maf_path $out_path
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"
tools_path="/storage/groups/ccm01/tools"

# 1. We first need to get MAF values for each SNP using PLINK
${tools_path}/PLINK/plink --bfile ${geno} --freq --out ${maf_path}/pmb_info06_all_chr_maf0.05_merged_noDup_maf_information

# generate GRanges objects needed for the functional annotation
Rscript ${script_path}/06_generate_GRanges_objects.R frq=\"${frq}\" geno=\"${geno}\" eqtls=\"${eqtl_path}\" clump=\"${clump_path}\" cmc=\"${cmc_path}\" gtex=\"${gtex_path}\" cortex=\"${cortex}\" vep=\"${vep}\" chromhmm=\"${chromhmm_path}\" gwas=\"${gwas_path}\" outpath=\"${out_path}\" &> ${out_path}/generate_GRanges_objects.log
