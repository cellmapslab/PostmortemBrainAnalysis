#!/bin/bash
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS"
# Generate the output directory (if not existing)
outpath_figures="${outpath}/03_figures"  # figures output path
mkdir -p $outpath_figures
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/01.GWAS"

# 1. Generate manhattan plot with covariates
gwas_results="${outpath}/02_GWAS_results/pmb_info06_all_chr_maf0.01_merged_samples_filtered_GWAS_with_covar.assoc.logistic"
output="${outpath_figures}/manhattan_plot_GWAS_logistic_with_covariates"
Rscript ${script_path}/03a_manhattan_plots.R input=\"${gwas_results}\" output=\"${output}\" 

# 2. Generate manhattan plot without covariates
gwas_results="${outpath}/02_GWAS_results/pmb_info06_all_chr_maf0.01_merged_samples_filtered_GWAS_without_covar.assoc"
output="${outpath_figures}/manhattan_plot_GWAS_assoc_without_covariates"
Rscript ${script_path}/03a_manhattan_plots.R input=\"${gwas_results}\" output=\"${output}\"

