#!/bin/bash

#### Script to run PRSice-2 with filtered SNPs ####

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS/05_makeenv.sh
script_path="/storage/groups/ccm01/tools/PRSice-2"
sumstats_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs"
target="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final"
ld_reference="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_cohort/1.4_dups_removed/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final"
phenotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/phenotypes.txt"
covariates="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/covariates.txt"
smooth_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/"

####################
### run PRSice-2 ###
####################
### SCZ2 GWAS ###
# 36,988 cases and 113,075 controls: N = 150,063 samples

# 1. Run for filtered SNPs nearby genes (1Mb)
# concatenate eQTL chunk files
#cat <(tail -n+2 ${smooth_path}/input/eQTLs/gene/01_pmb_CASE_CTL_matrixEQTL_cis_results_1_genes_chunk1.txt) <(tail -n+2 ${smooth_path}/input/eQTLs/gene/01_pmb_CASE_CTL_matrixEQTL_cis_results_1_genes_chunk2.txt) > ${smooth_path}/input/eQTLs/01_pmb_CASE_CTL_matrixEQTL_cis_results_1_genes.txt
# get SNPs nearby genes (1Mb)
#cut -f 1 ${smooth_path}/input/eQTLs/01_pmb_CASE_CTL_matrixEQTL_cis_results_1_genes.txt | sort | uniq > ${smooth_path}/input/eQTLs/snps_nearby_genes_1mb.txt
#Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/SCZ2.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --all-score --clump-kb 1000 --clump-r2 0.2 --score avg --extract ${smooth_path}/input/eQTLs/snps_nearby_genes_1mb.txt --out ${smooth_path}/results/nearby_genes/pmb_nearby_genes_1mb_scz2_PRSice-2_PRS

# 2. Run for gene cis-eQTL hits cases vs. controls
#Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/SCZ2.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --all-score --clump-kb 1000 --clump-r2 0.2 --score avg --extract ${smooth_path}/input/eQTLs/snps_gene_cis-eqtl_hits.txt --out ${smooth_path}/results/gene_eqtl/pmb_gene_eSNPs_cases_vs_controls_scz2_PRSice-2_PRS

# 2. Run for gene cis-eQTL hits cases only
#Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/SCZ2.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --all-score --clump-kb 1000 --clump-r2 0.2 --score avg --extract ${smooth_path}/input/eQTLs/snps_gene_cis-eqtl_hits_cases.txt --out ${smooth_path}/results/gene_eqtl/pmb_gene_eSNPs_cases_only_scz2_PRSice-2_PRS

# GTEx whole blood cis-eQTLs
#Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/SCZ2.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --all-score --clump-kb 1000 --clump-r2 0.2 --score avg --extract ${smooth_path}/input/eQTLs/gtex_wb_gene_cis-eqtls.txt --out ${smooth_path}/results/gtex_wb/pmb_gtex_gene_eSNPs_scz2_PRSice-2_PRS

/storage/groups/ccm01/tools/PLINK/plink --bfile ${target} --score /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/eQTLs/ld_snps_nearby_genes_1mb.txt 

### BIP 2018 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
#Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/BIP_2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_filtered_PRSice-2_PRS &> ${out_path}/BIP-2018/PRSice-2
#wait




