#!/bin/bash
snvs="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/01_RData_file/Index-patients_rare_variants_filtered.rds"
eqts="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/04_FUMA/4.1_input_files/3_exon_eQTS_ADHD_BIP_MDD_SCZ.txt"
eqtl="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/04_result_tables/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"
gencode="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/02_overlap"

# Generate the output directory (if not existing)
mkdir -p ${outpath}
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Index"

# Generate all input files necessary for FUMA
Rscript ${script_path}/02_generate_gene_overlap_file.R snvs=\"$snvs\" eqts=\"$eqts\" eqtl=\"$eqtl\" gencode=\"$gencode\" outpath=\"$outpath\" &> ${outpath}/generate_gene_overlap_file.log
