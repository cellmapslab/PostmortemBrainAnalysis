#!/bin/bash

# load R version 3.5.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_makeenv.sh

# run matrixeqtl for gene based cis-/trans-eQTLs 
#R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_run_matrix_eqtl_gene_summaries.R
#wait 

for value in {1..54} 
do
    file="/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/input/expression/01_pmb_exon_expr_SCZ_CTL_chunk${value}.txt"
    cis="/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/01_pmb_SCZ_CTL_matrixEQTL_cis_results_exons_chunk${value}.txt"
    trans="/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/01_pmb_SCZ_CTL_matrixEQTL_trans_results_exons_chunk${value}.txt"
    qsub -e "." -o "." -hard -l job_mem=100G -q long_fed25 /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/wrapper_matrix_eqtl_exon_summaries.sh $file $cis $trans
done

