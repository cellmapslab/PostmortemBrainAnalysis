#! /bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/02_makeenv.sh

file=$1
val=$2

R CMD BATCH --slave --no-save --no-restore "--args chunk=\"${file}\" value=\"${val}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/02a_matrix_eqtl_exon_summaries.R 02a_matrix_eqtl_exon_summaries_${val}.Rout

