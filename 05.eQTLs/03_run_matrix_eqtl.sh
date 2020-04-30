#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03a_makeenv.sh

# get all important paths
input_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input"
# Generate the output directory (if not existing)
output_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output" # output path to store matrixeqtl input files
mkdir -p $output_path
genespos="${input_path}/gene_locations.txt"
transpos="${input_path}/transcript_locations.txt"
exonspos="${input_path}/exon_locations.txt"
snpspos="${input_path}/SNP_locations.txt"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs"

##################
### Gene-level ###
##################
# run matrixeqtl for gene based cis-eQTLs for CASE & CTL
snps="${input_path}/01_pmb_snps_CASE_CTL.txt"
exp="${input_path}/01_pmb_gene_expr_CASE_CTL.txt"
cov="${input_path}/01_pmb_cov_CASE_CTL.txt"
ciso="${output_path}/01_pmb_CASE_CTL_matrixEQTL_cis_output_genes.txt"
transo="${output_path}/01_pmb_CASE_CTL_matrixEQTL_trans_output_genes.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${genespos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/01_matrix_eqtl_gene_summaries_CASE_CTL.log
wait

# run matrixeqtl for gene based cis-eQTLs for CASE only
snps="${input_path}/02_pmb_snps_CASE_only.txt"
exp="${input_path}/02_pmb_gene_expr_CASE_only.txt"
cov="${input_path}/02_pmb_cov_CASE_only.txt"
ciso="${output_path}/02_pmb_CASE_only_matrixEQTL_cis_output_genes.txt"
transo="${output_path}/02_pmb_CASE_only_matrixEQTL_trans_output_genes.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${genespos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/02_matrix_eqtl_gene_summaries_CASE_only.log
wait

# run matrixeqtl for gene based cis-eQTLs for CTL only
snps="${input_path}/03_pmb_snps_CTL_only.txt"
exp="${input_path}/03_pmb_gene_expr_CTL_only.txt"
cov="${input_path}/03_pmb_cov_CTL_only.txt"
ciso="${output_path}/03_pmb_CTL_only_matrixEQTL_cis_output_genes.txt"
transo="${output_path}/03_pmb_CTL_only_matrixEQTL_trans_output_genes.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${genespos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/03_matrix_eqtl_gene_summaries_CTL_only.log
wait


########################
### Transcript-level ###
########################
# run matrixeqtl for transcript based cis-eQTLs for CASE & CTL
snps="${input_path}/01_pmb_snps_CASE_CTL.txt"
exp="${input_path}/01_pmb_transcript_expr_CASE_CTL.txt"
cov="${input_path}/01_pmb_cov_CASE_CTL.txt"
ciso="${output_path}/01_pmb_CASE_CTL_matrixEQTL_cis_output_transcripts.txt"
transo="${output_path}/01_pmb_CASE_CTL_matrixEQTL_trans_output_transcripts.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${transpos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/01_matrix_eqtl_transcript_summaries_CASE_CTL.log
wait

# run matrixeqtl for transcript based cis-eQTLs for CASE only
snps="${input_path}/02_pmb_snps_CASE_only.txt"
exp="${input_path}/02_pmb_transcript_expr_CASE_only.txt"
cov="${input_path}/02_pmb_cov_CASE_only.txt"
ciso="${output_path}/02_pmb_CASE_only_matrixEQTL_cis_output_transcripts.txt"
transo="${output_path}/02_pmb_CASE_only_matrixEQTL_trans_output_transcripts.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${transpos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/02_matrix_eqtl_transcript_summaries_CASE_only.log
wait

# run matrixeqtl for transcript based cis-eQTLs for CTL only
snps="${input_path}/03_pmb_snps_CTL_only.txt"
exp="${input_path}/03_pmb_transcript_expr_CTL_only.txt"
cov="${input_path}/03_pmb_cov_CTL_only.txt"
ciso="${output_path}/03_pmb_CTL_only_matrixEQTL_cis_output_transcripts.txt"
transo="${output_path}/03_pmb_CTL_only_matrixEQTL_trans_output_transcripts.txt"
Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${exp}\" expos=\"${transpos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/03_matrix_eqtl_transcript_summaries_CTL_only.log
wait


##################
### Exon-level ###
##################
# generate log file
touch "/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03_matrix_eqtl_exon_summaries.Rout"
# run matrixeqtl for exon chunk based cis-eQTLs for CASE & CTL
for value in {1..5} 
do
    snps="${input_path}/01_pmb_snps_CASE_CTL.txt"
    chunk="${input_path}/01_pmb_exon_expr_CASE_CTL_chunk${value}.txt"
    cov="${input_path}/01_pmb_cov_CASE_CTL.txt"
    ciso="${output_path}/01_pmb_CASE_CTL_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="${output_path}/01_pmb_CASE_CTL_matrixEQTL_trans_output_exons_chunk${value}.txt"
    Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${chunk}\" expos=\"${exonspos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/01_matrix_eqtl_exon_summaries_CASE_CTL_chunk${value}.log
done
wait

# run matrixeqtl for exon chunk based cis-eQTLs for CASE only
for value in {1..5} 
do
    snps="${input_path}/02_pmb_snps_CASE_only.txt"
    chunk="${input_path}/02_pmb_exon_expr_CASE_only_chunk${value}.txt"
    cov="${input_path}/02_pmb_cov_CASE_only.txt"
    ciso="${output_path}/02_pmb_CASE_only_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="${output_path}/02_pmb_CASE_only_matrixEQTL_trans_output_exons_chunk${value}.txt"
    Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${chunk}\" expos=\"${exonspos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/02_matrix_eqtl_exon_summaries_CASE_only_chunk${value}.log
done
wait

# run matrixeqtl for exon chunk based cis-eQTLs for CTL only
for value in {1..5} 
do
    snps="${input_path}/03_pmb_snps_CTL_only.txt"
    chunk="${input_path}/03_pmb_exon_expr_CTL_only_chunk${value}.txt"
    cov="${input_path}/03_pmb_cov_CTL_only.txt"
    ciso="${output_path}/03_pmb_CTL_only_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="${output_path}/03_pmb_CTL_only_matrixEQTL_trans_output_exons_chunk${value}.txt"
    Rscript ${script_path}/03b_matrix_eqtl.R snps=\"${snps}\" spos=\"${snpspos}\" exp=\"${chunk}\" expos=\"${exonspos}\" cov=\"${cov}\" ciso=\"${ciso}\" transo=\"${transo}\" &>> ${output_path}/03_matrix_eqtl_exon_summaries_CTL_only_chunk${value}.log
done

