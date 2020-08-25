#!/bin/bash
inpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData_files"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/08_SNP-based_annotation"

# Generate the output directory (if not existing)
mkdir -p ${outpath}
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# Run SNP-based annotation for ChromHMM, VEP, GWAS, GTEx, CMC and meQTL
# VEP
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" which=\"test\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation.log
