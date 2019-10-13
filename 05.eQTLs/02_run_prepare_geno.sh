#!/bin/bash
gen_file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz"
sample_file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample"
eqtl_samples="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/eqtl_samples.txt"
output_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_probs_changed"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# generate a file including only the IIDs for the 169 samples with gene expression data (and covariates) needed for the input file generation
R CMD BATCH ${script_path}/02a_get_eqtl_samples.R

# echo "Start matrixEQTL genotype input transformation"
R CMD BATCH --slave "--args geno=\"${gen_file}\" sample=\"${sample_file}\" ids=\"${eqtl_samples}\" output=\"${output_path}\"" ${script_path}/02b_prepare_geno.R

