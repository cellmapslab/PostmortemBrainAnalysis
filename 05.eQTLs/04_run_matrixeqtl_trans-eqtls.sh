#!/bin/bash

# run matrixeqtl for gene based trans-eQTLs 
R CMD BATCH 04_run_matrix_eqtl_gene_summaries.R

# run matrixeqtl for exon based trans-eQTLs
R CMD BATCH 04_run_matrix_eqtl_exon_summaries.R

