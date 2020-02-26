#!/bin/bash
exclude_file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/01_preprocessed_genotypes/exclude_samples.txt"
expr_samples_file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/01_preprocessed_genotypes/samples_with_expression.txt"
sample_sheet="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2019-04-09_sample_sheet_KW.txt"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/"
output="${outpath}/samples_with_expression.txt"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
exclude="${outpath}/exclude_samples.txt"
gen_file="${outpath}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.gen.gz"
sample_file="${outpath}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.sample"
outpath2="${outpath}/1.5_probs_changed"
tools_path="/storage/groups/ccm01/tools"

# 1. Copy paste some necessary files from the GWAS analysis to the eQTL analysis
cp ${exclude_file} $outpath
cp ${expr_samples_file} $outpath

# 2. Preprocess the genotypes: filter for MAF, HWE, merge chromosome files and remove duplicates
bash ${script_path}/01a_prepare_genotypes.sh -g ${genotypes} -e ${exclude} -i ${output} -o ${outpath} -t ${tools_path}

# 3. Change the probabilities of the genotype file (three columns/numbers to one)
R CMD BATCH --slave "--args geno=\"${gen_file}\" sample=\"${sample_file}\" ids=\"${output}\" output=\"${outpath2}\"" ${script_path}/01b_change_probs.R
cp ${sample_file} $outpath2
mv $outpath2/*.sample pmb_info06_all_chr_maf0.05_merged_noDup_matrixEQTL_transformed.sample

