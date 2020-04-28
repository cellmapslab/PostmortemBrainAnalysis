#!/bin/bash

# Define neccessary paths
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"
tools_path="/storage/groups/ccm01/tools"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS"
# Generate the output directory (if not existing already)
outpath_ldpred="${outpath}/06_LDpred"  # LDpred output path
mkdir -p $outpath_ldpred
sumstats_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/03_preprocessed_sumstats"
target="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/01_preprocessed_genotypes/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final"
ld_reference="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/02_preprocessed_cohort/2.4_dups_removed/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final"
phenotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/phenotypes.txt"
covariates="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/covariates_ldpred.txt"

##################
### run LDpred ###
##################
# As LDpred takes a lot of memory and time, we will run the most important GWAS sumstats in following order: 
# SCZ2-2014, BIP-2018, MDD3-2019, CDG2-2018, ADHD-2017, ASD-2019, MDD2-2018, MDD1-2012, CDG1-2013 and then in alphabetical order 
# the not so important GWAS summary statistics

### SCZ2-2014 ### 
# 36,989 SCZ cases and 113,075 controls: N 150,064 samples
mkdir -p ${outpath_ldpred}/SCZ2-2014  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/SCZ2-2014.sumstats --nums 150064 --rsid snpid --ref a1 --alt a2 --pos bp --chr hg19chrc --pval p --eff or --out ${outpath_ldpred}/SCZ2-2014
wait
echo "Done SCZ2-2014"

### BIP-2018 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
mkdir -p ${outpath_lassosum}/BIP-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/BIP-2018.sumstats --nums 51710 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/BIP-2018
wait
echo "Done BIP-2018"

### MDD3-2019 ###
# 170,756 cases and 329,443 controls: N = 500,199 samples
mkdir -p ${outpath_lassosum}/MDD3-2019  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/MDD3-2019.sumstats --nums 500199 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/MDD3-2019
wait
echo "Done MDD3-2019"

### CDG2-2018 ###
# N = 438,997 samples
mkdir -p ${outpath_lassosum}/CDG2-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/CDG2-2018.sumstats --nums 438997 --rsid ID --ref REF --alt ALT --pos POS --chr CHROM --pval PVAL --eff BETA --beta yes --stderr SE --out ${outpath_ldpred}/CDG2-2018
wait
echo "Done CDG2-2018"

### ADHD-2017 ###
# 20,183 cases and 35,191 controls: N = 55,374 samples
mkdir -p ${outpath_lassosum}/ADHD-2017  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/ADHD-2017.sumstats --nums 55374 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/ADHD-2017
wait
echo "Done ADHD-2017"

### ASD-2019 ###
# 18,382 cases and 27,969 controls: N = 46,351 samples
mkdir -p ${outpath_lassosum}/ASD-2019  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/ASD-2019.sumstats --nums 46351 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/ASD-2019
wait
echo "Done ASD-2019"

### MDD2-2018 ###
# 135,458 cases and 344,901 control: N = 480,359 samples
mkdir -p ${outpath_lassosum}/MDD2-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/MDD2-2018.sumstats --nums 480359 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/MDD2-2018
wait
echo "Done MDD2-2018"

### MDD1-2012 ###
# 9,240 cases and 9,519 controls: N = 18,759 samples
mkdir -p ${outpath_lassosum}/MDD1-2012  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/MDD1-2012.sumstats --nums 18759 --rsid snpid --ref a1 --alt a2 --pos bp --chr hg18chr --pval pval --eff or --out ${outpath_ldpred}/MDD1-2012
wait
echo "Done MDD1-2012"

### CDG1-2013 ###
# 33,332 cases and 27,888 controls: N = 61,220 samples
mkdir -p ${outpath_lassosum}/CDG1-2013  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/CDG1-2013.sumstats --nums 61220 --rsid snpid --ref a1 --alt a2 --pos bp --chr hg18chr --pval pval --eff or --out ${outpath_ldpred}/CDG1-2013
wait
echo "Done CDG1-2013"

### AD-2019 ###
# 71,880 cases and 383,378 controls: N = 455,258 samples
mkdir -p ${outpath_lassosum}/AD-2019  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/AD-2019.sumstats --nums 455258 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff BETA --beta yes --stderr SE --out ${outpath_ldpred}/AD-2019
wait
echo "Done AD-2019"

### ALCDEP-2018 ###
# 11,476 cases and 23,080 controls: N = 34,556 samples
mkdir -p ${outpath_lassosum}/ALCDEP-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/ALCDEP-2018.sumstats --nums 34556 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff BETA --beta yes --stderr SE --out ${outpath_ldpred}/ALCDEP-2018
wait
echo "Done ALCDEP-2018"

### AN-2019 ###
# 16,992 cases and 55,525 controls: N = 72,517 samples
mkdir -p ${outpath_lassosum}/AN-2019  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/AN-2019.sumstats --nums 72517 --rsid ID --ref REF --alt ALT --pos POS --chr CHROM --pval PVAL --eff BETA --beta yes --stderr SE --out ${outpath_ldpred}/AN-2019
wait
echo "Done AN-2019"

### AUDIT-2018 ###
# NO bp positions given => NO LDpred
# N = 121,604 samples (no case/control distinction here!)
#mkdir -p ${outpath_lassosum}/AUDIT-2018  # generate output path if not existing already
#bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/AUDIT-2018.sumstats --nums 121604 --rsid rsid --ref a_0 --alt a_1 --pos bp --chr chr --pval p_T --eff beta_T --out ${outpath_ldpred}/AUDIT-2018
#wait
#echo "Done AUDIT-2018"

### BIPvsSCZ-2018 ###
# 15,270 BIP cases and 23,585 SCZ cases: N = 38,855 samples
mkdir -p ${outpath_lassosum}/BIPvsSCZ-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/BIPvsSCZ-2018.sumstats --nums 38855 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/BIPvsSCZ-2018
wait
echo "Done BIPvsSCZ-2018"

### BIPSCZvsCONT-2018 ###
# 20,129 BIP cases combined with 33,426 schizophrenia cases and 54,065 controls: N = 107,620 samples
mkdir -p ${outpath_lassosum}/BIPSCZvsCONT-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/BIPSCZvsCONT-2018.sumstats --nums 107620 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/BIPSCZvsCONT-2018
wait
echo "Done BIPSCZvsCONT-2018"

### BMI-2015 ###
# N = 87,048 samples
mkdir -p ${outpath_lassosum}/BMI-2015  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/BMI-2015.sumstats --nums 87048 --rsid SNP --ref effect_allele --alt other_allele --pos position --chr chromosome --pval P --eff effect --out ${outpath_ldpred}/BMI-2015
wait
echo "Done BMI-2015"

### CAD-2017 ###
# 10,801 cases (incl. angina) and 137,914 controls: N = 148,715 samples
mkdir -p ${outpath_lassosum}/CAD-2017  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/CAD-2017.sumstats --nums 148715 --rsid snptestid --ref effect_allele --alt noneffect_allele --pos bp_hg19 --chr chr --pval pvalue_gc --eff OR --out ${outpath_ldpred}/CAD-2017
wait
echo "Done CAD-2017"

### Cannabis-2018 ###
# 35,297 ICC and 22,683 23andMe and 126,785 UK-Biobank: N = 184,765 samples
mkdir -p ${outpath_lassosum}/Cannabis-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/Cannabis-2018.sumstats --nums 184765 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff BETA --beta yes --stderr SE --out ${outpath_ldpred}/Cannabis-2018
wait
echo "Done Cannabis-2018"

### EA-2018 (Educational attainment or IQ) ###
# N = 766,345 samples
mkdir -p ${outpath_lassosum}/EA-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/EA-2018.sumstats --nums 766345 --rsid MarkerName --ref A1 --alt A2 --pos POS --chr CHR --pval Pval --eff Beta --beta yes --stderr SE --out ${outpath_ldpred}/EA-2018
wait
echo "Done EA-2018"

### ED-2017 ###
# 3,495 cases and 10,982 controls: N = 14,477 samples
mkdir -p ${outpath_lassosum}/ED-2017  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/ED-2017.sumstats --nums 14477 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/ED-2017
wait
echo "Done ED-2017"

### Height-2014 ###
# N = 253,288 samples
mkdir -p ${outpath_lassosum}/Height-2014  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/Height-2014.sumstats --nums 253288 --rsid MarkerName --ref Allele1 --alt Allele2 --pos Position --chr Chromosome --pval p --eff b --beta yes --stderr SE --out ${outpath_ldpred}/Height-2014
wait
echo "Done Height-2014"

### Intelligence-2018 ###
# N = 269,867 samples
mkdir -p ${outpath_lassosum}/Intelligence-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/Intelligence-2018.sumstats --nums 269867 --rsid SNP --ref A1 --alt A2 --pos POS --chr CHR --pval P --eff stdBeta --beta yes --stderr SE --out ${outpath_ldpred}/Intelligence-2018
wait
echo "Done Intelligence-2018"

### OCD-2017 ###
# 2,688 cases and 7,037 controls: N = 9,725 samples
mkdir -p ${outpath_lassosum}/OCD-2017  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/OCD-2017.sumstats --nums 9725 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/OCD-2017
wait
echo "Done OCD-2017"

### PTSD-2019 ###
# 23,212 PTSD cases and 151,447 controls: N = 174,659 samples
mkdir -p ${outpath_lassosum}/PTSD-2019  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/PTSD-2019.sumstats --nums 174659 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/PTSD-2019
wait
echo "Done PTSD-2019"

### T2D-2017 ###
# 26,676 T2D cases and 132,532 controls: N = 159,208 samples
mkdir -p ${outpath_lassosum}/T2D-2017  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/T2D-2017.sumstats --nums 159208 --rsid Chr:Position --ref Allele1 --alt Allele2 --pos Position --chr Chr --pval P-value --eff Effect --beta yes --stderr StdErr --out ${outpath_ldpred}/T2D-2017
echo "Done T2D-2017"

### TS-2018 ###
# 4,819 TS cases and 9,488 controls: N = 14,307 samples
mkdir -p ${outpath_lassosum}/TS-2018  # generate output path if not existing already
bash ${script_path}/06a_LDpred.sh --genotypes ${target} --phenotypes ${phenotypes} --covariates ${covariates} --cohort ${ld_reference} --sumstats ${sumstats_path}/TS-2018.sumstats --nums 14307 --rsid SNP --ref A1 --alt A2 --pos BP --chr CHR --pval P --eff OR --out ${outpath_ldpred}/TS-2018
wait
echo "Done TS-2018"

