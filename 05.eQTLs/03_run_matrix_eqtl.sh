#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03a_makeenv.sh

##################
### Gene-level ###
##################
# run matrixeqtl for gene based cis-eQTLs for CASE & CTL
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_snps_CASE_CTL.txt"
expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_gene_expr_CASE_CTL.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_cov_CASE_CTL.txt"
genespos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/gene_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_cis_output_genes.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_trans_output_genes.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${expr}\" cov=\"${cov}\" value=\"${value}\" pos=\"${genespos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R 03_matrix_eqtl_gene_summaries.Rout
wait

# run matrixeqtl for gene based cis-eQTLs for CASE only
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_snps_CASE_only.txt"
expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_gene_expr_CASE_only.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_cov_CASE_only.txt"
genespos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/gene_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_cis_output_genes.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_trans_output_genes.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${expr}\" cov=\"${cov}\" value=\"${value}\" pos=\"${genespos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_gene_summaries.Rout
wait

# run matrixeqtl for gene based cis-eQTLs for CTL only
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_snps_CTL_only.txt"
expr="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_gene_expr_CTL_only.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_cov_CTL_only.txt"
genespos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/gene_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_cis_output_genes.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_trans_output_genes.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${expr}\" cov=\"${cov}\" value=\"${value}\" pos=\"${genespos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_gene_summaries.Rout
wait

########################
### Transcript-level ###
########################
# run matrixeqtl for transcript based cis-eQTLs for CASE & CTL
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_snps_CASE_CTL.txt"
chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_transcript_expr_CASE_CTL.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_cov_CASE_CTL.txt"
transpos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/transcript_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_cis_output_transcripts.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_trans_output_transcripts.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" chunk=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${transpos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R 03_matrix_eqtl_transcript_summaries.Rout
wait

# run matrixeqtl for transcript based cis-eQTLs for CASE only
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_snps_CASE_only.txt"
chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_transcript_expr_CASE_only.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_cov_CASE_only.txt"
transpos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/transcript_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_cis_output_transcripts.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_trans_output_transcripts.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" chunk=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${transpos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_transcript_summaries.Rout
wait

# run matrixeqtl for transcript based cis-eQTLs for CTL only
snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_snps_CTL_only.txt"
chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_transcript_expr_CTL_only.txt"
cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_cov_CTL_only.txt"
transpos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/transcript_locations.txt"
ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_cis_output_transcripts.txt"
transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_trans_output_transcripts.txt"
R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" chunk=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${transpos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_transcript_summaries.Rout
wait

##################
### Exon-level ###
##################
# generate log file
touch "/storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03_matrix_eqtl_exon_summaries.Rout"
# run matrixeqtl for exon chunk based cis-eQTLs for CASE & CTL
for value in {1..5} 
do
    snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_snps_CASE_CTL.txt"
    chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_exon_expr_CASE_CTL_chunk${value}.txt"
    cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_cov_CASE_CTL.txt"
    exonspos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/exon_locations.txt"
    ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/01_pmb_CASE_CTL_matrixEQTL_trans_output_exons_chunk${value}.txt"
    R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${exonspos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_exon_summaries.Rout
done
wait

# run matrixeqtl for exon chunk based cis-eQTLs for CASE only
for value in {1..5} 
do
    snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_snps_CASE_only.txt"
    chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_exon_expr_CASE_only_chunk${value}.txt"
    cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_cov_CASE_only.txt"
    exonspos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/exon_locations.txt"
    ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/02_pmb_CASE_only_matrixEQTL_trans_output_exons_chunk${value}.txt"
    R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${exonspos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_exon_summaries.Rout
done
wait

# run matrixeqtl for exon chunk based cis-eQTLs for CTL only
for value in {1..5} 
do
    snps="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_snps_CTL_only.txt"
    chunk="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_exon_expr_CTL_only_chunk${value}.txt"
    cov="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_cov_CTL_only.txt"
    exonspos="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/exon_locations.txt"
    ciso="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_cis_output_exons_chunk${value}.txt"
    transo="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/03_matrixeqtl_output/03_pmb_CTL_only_matrixEQTL_trans_output_exons_chunk${value}.txt"
    R CMD BATCH --slave --no-save --no-restore "--args snps=\"${snps}\" expr=\"${chunk}\" cov=\"${cov}\" value=\"${value}\" pos=\"${exonspos}\" ciso=\"${ciso}\" transo=\"${transo}\"" /storage/groups/ccm01/code/PostmortemBrainAnalysis/05.eQTLs/03b_matrix_eqtl.R >> 03_matrix_eqtl_exon_summaries.Rout
done
