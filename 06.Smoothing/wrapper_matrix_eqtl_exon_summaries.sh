#! /bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_makeenv.sh

file=$1
cis=$2
trans=$3

R CMD BATCH --slave --no-save --no-restore "--args chunk=\"${file}\" cis=\"${cis}\" trans=\"${trans}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01b_matrix_eqtl_exon_summaries.R

