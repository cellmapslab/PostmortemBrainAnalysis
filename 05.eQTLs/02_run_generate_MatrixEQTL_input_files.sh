#!/bin/bash
bim="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.bim"
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.5_probs_changed/pmb_info06_all_chr_maf0.05_merged_noDup_matrixEQTL_transformed.gen"
gene_expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData"
transcript_expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData"
exon_expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData"
# Generate the output directory (if not existing)
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input" # output path to store matrixeqtl input files
mkdir -p $outpath
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# 1. Generate all input files necessary for matrixeqtl
Rscript ${script_path}/02_generate_MatrixEQTL_input_files.R bim=\"${bim}\" snps=\"${snps}\" gene=\"${gene_expr}\" transcript=\"${transcript_expr}\" exon=\"${exon_expr}\" outpath=\"${outpath}\" &> ${outpath}/generate_MatrixEQTL_input.log
