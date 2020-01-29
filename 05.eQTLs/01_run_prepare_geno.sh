#!/bin/bash
exclude_file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/01_preprocessed_genotypes/exclude_samples.txt"
sample_sheet="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2019-04-09_sample_sheet_KW.txt"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes"
output="${outpath}/eqtl_samples.txt"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
exclude="${outpath}/exclude_samples.txt"
gen_file="${outpath}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.gen.gz"
sample_file="${outpath}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.sample"
outpath2="${outpath}/1.5_probs_changed"
tools_path="/storage/groups/ccm01/tools"

# 1. Generate some necessary files for the eQTL analysis
R CMD BATCH --slave "--args input1=\"${exclude_file}\" input2=\"${sample_sheet}\" outpath=\"${outpath}\" output=\"${output}\"" ${script_path}/01a_generate_necessary_files.R

# 2. Preprocess the genotypes: filter for MAF, HWE, merge chromosome files and remove duplicates
bash ${script_path}/01b_prepare_genotypes.sh -g ${genotypes} -e ${exclude} -i ${output} -o ${outpath} -t ${tools_path}

# 3. Change the probabilities of the genotype file (three columns/numbers to one)
#R CMD BATCH --slave "--args geno=\"${gen_file}\" sample=\"${sample_file}\" ids=\"${output}\" output=\"${outpath2}\"" ${script_path}/01c_change_probs.R

