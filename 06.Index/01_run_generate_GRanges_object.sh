#!/bin/bash
genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
sample_file="${genotypes}/pmb_info06_chr1.sample"
sample_sheet="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2019-04-09_sample_sheet_KW.txt"
inpath=""
# Generate the output directory (if not existing)
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/01_preprocessed_genotypes"  # preprocessed genotypes output path
mkdir -p $outpath_geno
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/06.index"

# 1. Generate all input files necessary for the preprocessing of the genotypes
Rscript ${script_path}/01a_generate_necessary_files.R input1=\"${sample_file}\" input2=\"${sample_sheet}\" output1=\"${output1}\" output2=\"${output2}\" &> ${outpath_geno}/generate_necessary_files.log

# 2. Preprocess the genotypes: filter for MAF, HWE, merge chromosome files, remove duplicates and change probabilities (three columns/numbers to one)
bash ${script_path}/01b_prepare_genotypes.sh -g ${genotypes} -e ${output1} -i ${output2} -o ${outpath_geno} -t ${tools_path} -s ${script_path} &> ${outpath_geno}/prepare_genotypes.log
