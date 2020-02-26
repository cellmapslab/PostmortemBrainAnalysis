#!/bin/bash
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
sample_file="${genotypes}/pmb_info06_chr1.sample"
sample_sheet="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2019-04-09_sample_sheet_KW.txt"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS"
outpath_geno="${outpath}/01_preprocessed_genotypes"
output1="${outpath_geno}/exclude_samples.txt"
output2="${outpath_geno}/samples_with_expression.txt"
output3="${outpath}/covar.txt"
output4="${outpath}/pheno.txt"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/01.GWAS"
tools_path="/storage/groups/ccm01/tools"

# 1. Generate all input files necessary for the preprocessing of the genotypes
R CMD BATCH --slave "--args input1=\"${sample_file}\" input2=\"${sample_sheet}\" output1=\"${output1}\" output2=\"${output2}\" output3=\"${output3}\" output4=\"${output4}\"" ${script_path}/01a_generate_necessary_files.R

# 2. Preprocess the genotypes: filter for MAF, HWE, merge chromosome files and remove duplicates
bash ${script_path}/01b_prepare_genotypes.sh -g ${genotypes} -e ${output1} -i ${output2} -o ${outpath_geno} -t ${tools_path}

