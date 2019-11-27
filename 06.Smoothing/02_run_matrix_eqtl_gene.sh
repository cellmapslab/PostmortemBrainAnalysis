#!/bin/bash
  
for value in {1..2}
do
    file="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/expression/gene/01_pmb_gene_probeset_expr_CASE_CTL_chunk${value}.txt"
    #qsub -e "." -o "." -hard -l job_mem=100G -q long_fed25 /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/02b_wrapper_matrix_eqtl_gene_summaries.sh $file $value
    bash /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/02b_wrapper_matrix_eqtl_gene_summaries.sh $file $value &
done

