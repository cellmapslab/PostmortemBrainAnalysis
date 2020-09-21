#!/bin/bash
gencode="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf"
inpath="/storage/groups/ccm01/datasets/2018_IntegraMent_Mental_Disorders/Schizophrenia/Index-Patients"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/01_RData_file"  # preprocessed genotypes output path
# Generate the output directory (if not existing)
mkdir -p $outpath
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Index"

# Generate GRanges object and tables
Rscript ${script_path}/01_generate_GRanges_object.R inpath=\"$inpath\" gencode=\"$gencode\" outpath=\"${outpath}\" &> ${outpath}/generate_GRanges_object.log
