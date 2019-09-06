#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/04_makeenv.sh

# run matrixeqtl for gene based cis-/trans-eQTLs 
R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/04_run_matrix_eqtl_gene_summaries.R 
#wait

# run matrixeqtl for exon based cis-/trans-eQTLs
#R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/04_run_matrix_eqtl_transcript_summaries.R
#wait 

# run matrixeqtl for gene based cis-/trans-eQTLs 
#R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/04_run_matrix_eqtl_exon_summaries.R 

