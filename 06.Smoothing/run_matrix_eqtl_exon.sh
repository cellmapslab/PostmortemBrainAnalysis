#!/bin/bash

for value in {1..25} 
do
    file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/expression/01_pmb_exon_expr_CASE_CTL_chunk${value}.txt"
    cis="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_CASE_CTL_matrixEQTL_cis_results_exons_chunk${value}.txt"
    trans="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/01_pmb_CASE_CTL_matrixEQTL_trans_results_exons_chunk${value}.txt"
    qsub -e "." -o "." -hard -l job_mem=100G -q long_fed25 /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/wrapper_matrix_eqtl_exon_summaries.sh $file $cis $trans
done

