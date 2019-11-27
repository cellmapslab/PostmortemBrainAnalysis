#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS/05_makeenv.sh
script_path="/storage/groups/ccm01/tools/PRSice-2"
out_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/PRSice-2"
sumstats_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.2_sumstats_changed_IDs"
target="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final"
ld_reference="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_cohort/1.4_dups_removed/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final"
phenotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/phenotypes.txt"
covariates="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/covariates.txt"

####################
### run PRSice-2 ###
####################
### SCZ2 GWAS ###
# 36,988 cases and 113,075 controls: N = 150,063 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/SCZ2.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/SCZ2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_scz2_filtered_PRSice-2_PRS &> ${out_path}/SCZ2/PRSice-2
# --no-regress outputs only PRS for all thresholds without the plots and the best score => better use --all-score to also get this file + plots + best PRS
wait

### MDD2 ###
# 135,458 cases and 344,901 controls: N = 480,359 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/MDD2_2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/MDD2/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_mdd2_filtered_PRSice-2_PRS &> ${out_path}/MDD2/PRSice-2
wait

### BIP 2018 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/BIP_2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_filtered_PRSice-2_PRS &> ${out_path}/BIP-2018/PRSice-2
wait

### Cross-Disorder ###
# 61,220 samples (33,332 cases and 27,888 controls) => PGC_OLD
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/Cross_disorder_hg19.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/Cross_Disorder/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_cd_filtered_PRSice-2_PRS &> ${out_path}/Cross_Disorder/PRSice-2
wait

### BIP+SCZvsCONT ###
# 20,129 BIP cases and 33,426 SCZ cases and 54,065 controls: N = 107,620 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/BIPSCZvsCONT.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/BIP+SCZvsCONT/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_scz_vs_cont_filtered_PRSice-2_PRS &> ${out_path}/BIP+SCZvsCONT/PRSice-2
wait

### BIPvsSCZ ###
# 15,270 BIP cases and 23,585 SCZ cases: N = 38,855 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/BIPvsSCZ.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/BIPvsSCZ/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bip_vs_scz_filtered_PRSice-2_PRS &> ${out_path}/BIPvsSCZ/PRSice-2
wait

### ADHD 2017 ###
# 20,183 cases and 35,191 controls: N = 55,374 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/ADHD_2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/ADHD-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_adhd_filtered_PRSice-2_PRS &> ${out_path}/ADHD-2017/PRSice-2
wait

### ASD ###
# 5,305 cases and 5305 pseudocontrols: N = 10,610 samples
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/ASD.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/ASD/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_asd_filtered_PRSice-2_PRS &> ${out_path}/ASD/PRSice-2
wait

### BMI ###
# 87,048 samples (no case/control here!)
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/BMI.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 effect_allele --A2 other_allele --stat effect --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/BMI/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_bmi_filtered_PRSice-2_PRS &> ${out_path}/BMI/PRSice-2
wait

### Educational Attainment (EA or IQ) ###
# 766,345 samples (no case/control, without 23andMe samples)
Rscript --slave ${script_path}/PRSice.R --prsice ${script_path}/PRSice_linux --base ${sumstats_path}/EA.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat Beta --pvalue Pval --snp MarkerName --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${out_path}/EA/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ea_filtered_PRSice-2_PRS &> ${out_path}/EA/PRSice-2
wait

