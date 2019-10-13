#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_makeenv.sh

# run matrixeqtl for gene based cis-/trans-eQTLs 
R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_matrix_eqtl_gene_summaries.R

# run matrixeqtl for transcript based cis-/trans-eQTLs 
R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_matrix_eqtl_transcript_summaries.R



