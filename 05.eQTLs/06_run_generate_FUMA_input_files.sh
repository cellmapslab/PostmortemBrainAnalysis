#!/bin/bash
input_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/04_result_tables"
dea_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_files"
output_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_FUMA"
kegg_file="/storage/groups/ccm01/datasets/2020_KEGG/2020-06-14_KEGG_entrez_no_disease_pathways.rds"

# Generate the output directory (if not existing)
mkdir -p ${output_path}
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

# 1. Generate all input files necessary for FUMA
Rscript ${script_path}/06_generate_FUMA_input_files.R inpath=\"${input_path}\" dea=\"${dea_path}\" keggf=\"${kegg_file}\" outpath=\"${output_path}\" &> ${output_path}/generate_FUMA_input.log
