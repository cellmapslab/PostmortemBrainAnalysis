#!/bin/bash
gencode="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf"
path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs"
genotypes="${path}/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup"
snpspos="${geno}.bim"
genespos="${path}/02_matrixeqtl_input/gene_locations.txt"
transpos="${path}/02_matrixeqtl_input/transcript_locations.txt"
exonspos="${path}/02_matrixeqtl_input/exon_locations.txt"
inpath="${path}/03_matrixeqtl_output"
# Generate the output directory (if not existing)
outpath_tables="${path}/04_result_tables"  # preprocessed genotypes output path
outpath_clumping="${path}/05_clumping"  # preprocessed genotypes output path
mkdir -p $outpath_tables $outpath_clumping
outpath_clumping_input="${outpath_clumping}/01_input_files" # path to store all input files necessary for the clumping process
outpath_clumping_results="${outpath_clumping}/02_clumping_results" # path to store all clumped result files
outpath_clumping_tables="${outpath_clumping}/03_result_tables" # path to store clumping eQTL result tables
mkdir -p $outpath_clumping_input $outpath_clumping_results $outpath_clumping_tables
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"
tools_path="/storage/groups/ccm01/tools"

# 1. Generate eQTL result tables from the matrixEQTL output files
#Rscript ${script_path}/04a_generate_result_tables.R inpath=\"${inpath}\" snps=\"${snpspos}\" gtf=\"${gencode}\" gpos=\"${genespos}\" tpos=\"${transpos}\" epos=\"${exonspos}\" outpath=\"${outpath_tables}\" &> ${outpath_tables}/generate_result_tables.log

# 2. Clumping of eSNPs and generation of clumped result tables
Rscript ${script_path}/04b_clumping.R inpath=\"${outpath_tables}\" geno=\"${genotypes}\" outpath1=\"${outpath_clumping_input}\" outpath2=\"${outpath_clumping_results}\" outpath3=\"${outpath_clumping_tables}\" &> ${outpath_clumping}/clumping.log
