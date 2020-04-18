#!/bin/bash

# load R version 3.6.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS/04a_makeenv.sh

# get all important paths and necessary files
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS"
# Generate the output directory (if not existing already)
outpath_lassosum="${outpath}/05_lassosum"  # lassosum output path
mkdir -p $outpath_lassosum
sumstats_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/03_preprocessed_sumstats"
target="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/01_preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final"
ld_reference="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/02_preprocessed_cohort/2.4_dups_removed/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final"
phenotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/phenotypes.txt"
covariates="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/covariates.txt"

####################
### run lassosum ###
####################
### AD-2019 ###
#mkdir -p ${outpath_lassosum}/AD-2019  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat BETA --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/AD-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_AD-2019_PRSice-2_PRS &> ${outpath_lassosum}/AD-2019/PRSice-2
#wait
#echo "Done AD-2019"

### ADHD-2017 ###
# 20,183 cases and 35,191 controls: N = 55,374
mkdir -p ${outpath_lassosum}/ADHD-2017  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ADHD-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 55374 --thread 3 --out ${outpath_lassosum}/ADHD-2017 --gwas ADHD-2017 &> ${outpath_lassosum}/ADHD-2017/PRSice-2
wait
echo "Done ADHD-2017"

### ALCDEP-2018 ###
#mkdir -p ${outpath_lassosum}/ALCDEP-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ALCDEP-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat BETA --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/ALCDEP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ALCDEP-2018_PRSice-2_PRS &> ${outpath_lassosum}/ALCDEP-2018/PRSice-2
#wait
#echo "Done ALCDEP-2018"

### AN-2019 ###
#mkdir -p ${outpath_lassosum}/AN-2019  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AN-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 REF --A2 ALT --stat BETA --pvalue PVAL --snp ID --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/AN-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_AN-2019_PRSice-2_PRS &> ${outpath_lassosum}/AN-2019/PRSice-2
#wait
#echo "Done AN-2019"

### ASD-2019 ###
#mkdir -p ${outpath_lassosum}/ASD-2019  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ASD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/ASD-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ASD-2019_PRSice-2_PRS &> ${outpath_lassosum}/ASD-2019/PRSice-2
#wait
#echo "Done ASD-2019"

### AUDIT-2018 ###
#mkdir -p ${outpath_lassosum}/AUDIT-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AUDIT-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death --ignore-fid --thread 1 --A1 a_0 --A2 a_1 --stat beta_T --pvalue p_T --snp rsid --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/AUDIT-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_AUDIT-2018_PRSice-2_PRS &> ${outpath_lassosum}/AUDIT-2018/PRSice-2
#wait
#echo "Done AUDIT-2018"

### BIP-2018 ###
#mkdir -p ${outpath_lassosum}/BIP-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIP-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/BIP-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BIP-2018_PRSice-2_PRS &> ${outpath_lassosum}/BIP-2018/PRSice-2
#wait
#echo "Done BIP-2018"

### BIPvsSCZ-2018 ###
#mkdir -p ${outpath_lassosum}/BIPvsSCZ-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIPvsSCZ-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/BIPvsSCZ-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BIPvsSCZ-2018_PRSice-2_PRS &> ${outpath_lassosum}/BIPvsSCZ-2018/PRSice-2
#wait
#echo "Done BIPvsSCZ-2018"

### BIPSCZvsCONT-2018 ###
#mkdir -p ${outpath_lassosum}/BIPSCZvsCONT-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIPSCZvsCONT-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/BIPSCZvsCONT-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BIPSCZvsCONT-2018_PRSice-2_PRS &> ${outpath_lassosum}/BIPSCZvsCONT-2018/PRSice-2
#wait
#echo "Done BIPSCZvsCONT-2018"

### BMI-2015 ###
#mkdir -p ${outpath_lassosum}/BMI-2015  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BMI-2015.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 effect_allele --A2 other_allele --stat effect --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/BMI-2015/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_BMI-2015_PRSice-2_PRS &> ${outpath_lassosum}/BMI-2015/PRSice-2
#wait
#echo "Done BMI-2015"

### CAD-2017 ###
#mkdir -p ${outpath_lassosum}/CAD-2017  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CAD-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 effect_allele --A2 noneffect_allele --stat OR --pvalue p-value_gc --snp snptestid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/CAD-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_CAD-2017_PRSice-2_PRS &> ${outpath_lassosum}/CAD-2017/PRSice-2
#wait
#echo "Done CAD-2017"

### Cannabis-2018 ###
#mkdir -p ${outpath_lassosum}/Cannabis-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Cannabis-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat BETA --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/Cannabis-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_Cannabis-2018_PRSice-2_PRS &> ${outpath_lassosum}/Cannabis-2018/PRSice-2
#wait
#echo "Done Cannabis-2018"

### CDG1-2013 ###
#mkdir -p ${outpath_lassosum}/CDG1-2013  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CDG1-2013.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue pval --snp snpid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/CDG1-2013/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_CDG1-2013_PRSice-2_PRS &> ${outpath_lassosum}/CDG1-2013/PRSice-2
#wait
#echo "Done CDG1-2013"

### CDG2-2018 ###
#mkdir -p ${outpath_lassosum}/CDG2-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CDG2-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 REF --A2 ALT --stat BETA --pvalue PVAL --snp ID --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/CDG2-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_CDG2-2018_PRSice-2_PRS &> ${outpath_lassosum}/CDG2-2018/PRSice-2
#wait
#echo "Done CDG2-2018"

### EA-2018 (Educational attainment or IQ) ###
#mkdir -p ${outpath_lassosum}/EA-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/EA-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat Beta --pvalue Pval --snp MarkerName --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/EA-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_EA-2018_PRSice-2_PRS &> ${outpath_lassosum}/EA-2018/PRSice-2
#wait
#echo "Done EA-2018"

### ED-2017 ###
#mkdir -p ${outpath_lassosum}/ED-2017  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ED-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/ED-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_ED-2017_PRSice-2_PRS &> ${outpath_lassosum}/ED-2017/PRSice-2
#wait
#echo "Done ED-2017"

### Height-2014 ###
#mkdir -p ${outpath_lassosum}/Height-2014  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Height-2014.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 Allele1 --A2 Allele2 --stat b --pvalue p --snp MarkerName --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/Height-2014/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_Height-2014_PRSice-2_PRS &> ${outpath_lassosum}/Height-2014/PRSice-2
#wait
#echo "Done Height-2014"

### Intelligence-2018 ###
#mkdir -p ${outpath_lassosum}/Intelligence-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Intelligence-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat stdBeta --pvalue P --snp SNP --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/Intelligence-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_Intelligence-2018_PRSice-2_PRS &> ${outpath_lassosum}/Intelligence-2018/PRSice-2
#wait
#echo "Done Intelligence-2018"

### MDD1-2012 ###
#mkdir -p ${outpath_lassosum}/MDD1-2012  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD1-2012.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue pval --snp snpid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/MDD1-2012/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD1-2012_PRSice-2_PRS &> ${outpath_lassosum}/MDD1-2012/PRSice-2
#wait
#echo "Done MDD1-2012"

### MDD2-2018 ###
#mkdir -p ${outpath_lassosum}/MDD2-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD2-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/MDD2-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD2-2018_PRSice-2_PRS &> ${outpath_lassosum}/MDD2-2018/PRSice-2
#wait
#echo "Done MDD2-2018"

### MDD3-2019 ###
#mkdir -p ${outpath_lassosum}/MDD3-2019  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD3-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/MDD3-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_MDD3-2019_PRSice-2_PRS &> ${outpath_lassosum}/MDD3-2019/PRSice-2
#wait
#echo "Done MDD3-2019"

### OCD-2017 ###
#mkdir -p ${outpath_lassosum}/OCD-2017  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/OCD-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/OCD-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_OCD-2017_PRSice-2_PRS &> ${outpath_lassosum}/OCD-2017/PRSice-2
#wait
#echo "Done OCD-2017"

### PTSD-2019 ###
#mkdir -p ${outpath_lassosum}/PTSD-2019  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/PTSD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/PTSD-2019/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_PTSD-2019_PRSice-2_PRS &> ${outpath_lassosum}/PTSD-2019/PRSice-2
#wait
#echo "Done PTSD-2019"

### SCZ2-2014 ### 
#mkdir -p ${outpath_lassosum}/SCZ2-2014  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/SCZ2-2014.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 a1 --A2 a2 --stat or --pvalue p --snp snpid --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/SCZ2-2014/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_SCZ2-2014_PRSice-2_PRS &> ${outpath_lassosum}/SCZ2-2014/PRSice-2
#wait
#echo "Done SCZ2-2014"

### T2D-2017 ###
#mkdir -p ${outpath_lassosum}/T2D-2017  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/T2D-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 Allele1 --A2 Allele2 --stat Effect --pvalue P-value --snp Chr:Position --binary-target F --print-snp --beta --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/T2D-2017/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_T2D-2017_PRSice-2_PRS &> ${outpath_lassosum}/T2D-2017/PRSice-2
#wait
#echo "Done T2D-2017"

### TS-2018 ###
#mkdir -p ${outpath_lassosum}/TS-2018  # generate output path if not existing already
#Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/TS-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --cov-factor Sex,Suicide,Type_of_Death  --ignore-fid --thread 1 --A1 A1 --A2 A2 --stat OR --pvalue P --snp SNP --binary-target T --print-snp --clump-kb 1000 --clump-r2 0.2 --score avg --all-score --out ${outpath_lassosum}/TS-2018/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_TS-2018_PRSice-2_PRS &> ${outpath_lassosum}/TS-2018/PRSice-2
#wait
#echo "Done TS-2018"
# --no-regress outputs only PRS for all thresholds without the plots and the best score => better use --all-score to also get this file + plots + best PRS

