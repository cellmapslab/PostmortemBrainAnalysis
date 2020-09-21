#!/bin/bash
inpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/02_overlap"
surselect="/storage/groups/ccm01/datasets/2018_IntegraMent_Mental_Disorders/Schizophrenia/2017-09-27_S04380110_Covered_Agilent_Surselect_V5.bed"
gencode="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf"
ebkgd="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_FUMA/7.1_input_files/3_exon_bkgd.txt"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/03_FUMA"
outpath1="$outpath/3.1_input_files"
outpath2="$outpath/3.2_rare_variants_index"
outpath3="$outpath/3.3_overlap_rare_eeGenes_eQTS"
kegg="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_FUMA/7.1_input_files/4_kegg_entrez_no_disease_pathways.txt"

# Generate the output directory (if not existing)
mkdir -p ${outpath} ${outpath1} ${outpath2} ${outpath3}
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Index"

# Copy kegg input file from eQTL to Index folder and rename it
cp $kegg "$outpath1"
mv "$outpath1/4_kegg_entrez_no_disease_pathways.txt" "$outpath1/3_kegg_entrez_no_disease_pathways.txt"

# Generate all input files necessary for FUMA
Rscript ${script_path}/03_generate_FUMA_input_files.R inpath=\"$inpath\" surselect=\"$surselect\" gencode=\"$gencode\" ebkgd=\"$ebkgd\" outpath=\"$outpath1\" &> ${outpath1}/generate_FUMA_input_files.log
