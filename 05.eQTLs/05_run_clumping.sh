#!/bin/bash
path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs"
genotypes="${path}/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup"
inpath="${path}/04_result_tables"  # preprocessed genotypes output path
# Generate the output directory (if not existing)
outpath="${path}/05_clumping"  # preprocessed genotypes output path
mkdir -p $outpath
outpath_input="${outpath}/01_input_files" # path to store all input files necessary for the clumping process
outpath_results="${outpath}/02_clumping_results" # path to store all clumped result files
outpath_tables="${outpath}/03_result_tables" # path to store clumping eQTL result tables
mkdir -p $outpath_input $outpath_results $outpath_tables
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# Clumping of eSNPs and generation of clumped result tables
Rscript ${script_path}/05_clumping.R inpath=\"${inpath}\" geno=\"${genotypes}\" outpath1=\"${outpath_input}\" outpath2=\"${outpath_results}\" outpath3=\"${outpath_tables}\" &> ${outpath}/clumping.log
