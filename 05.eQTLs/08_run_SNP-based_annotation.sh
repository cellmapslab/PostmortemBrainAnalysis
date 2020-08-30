#!/bin/bash
inpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData_files"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/08_SNP-based_annotation"

# Generate the output directory (if not existing)
mkdir -p ${outpath}
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# Run SNP-based annotation for ChromHMM, VEP, GWAS, GTEx, CMC and meQTL
# VEP
#Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"VEP\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_VEP.log
# ChromHMM Blood
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"ChromHMM_Blood\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_ChromHMM_Blood.log
# ChromHMM Brain
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"ChromHMM_Brain\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_ChromHMM_Brain.log
# ChromHMM_DLPFC
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"ChromHMM_DLPFC\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_ChromHMM_DLPFC.log
# GTEx_Blood
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"GTEx_Blood\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_GTEx_Blood.log
# GTEx_Brain
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"GTEx_Brain\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_GTEx_Brain.log
# GTEx_BA9
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"GTEx_BA9\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_GTEx_BA9.log
# CMC
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"CMC\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_CMC.log
# meQTL
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"meQTL\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_Cortex_meQTL.log
# GWAS
Rscript ${script_path}/08_SNP-based_annotation.R inpath=\"${inpath}\" what=\"GWAS\" outpath=\"${outpath}\" &> ${outpath}/SNP-based_annotation_GWAS.log
