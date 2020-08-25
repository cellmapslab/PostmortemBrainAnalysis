#!/bin/bash
gencode="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf"
path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs"
snpspos="${path}/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.bim"
genespos="${path}/02_matrixeqtl_input/gene_locations.txt"
transpos="${path}/02_matrixeqtl_input/transcript_locations.txt"
exonspos="${path}/02_matrixeqtl_input/exon_locations.txt"
inpath="${path}/03_matrixeqtl_output"
# Generate the output directory (if not existing)
outpath="${path}/04_result_tables"  # preprocessed genotypes output path
mkdir -p $outpath
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# Generate eQTL result tables from the matrixEQTL output files
Rscript ${script_path}/04_generate_result_tables.R inpath=\"${inpath}\" snps=\"${snpspos}\" gtf=\"${gencode}\" gpos=\"${genespos}\" tpos=\"${transpos}\" epos=\"${exonspos}\" outpath=\"${outpath}\" &> ${outpath}/generate_result_tables.log
