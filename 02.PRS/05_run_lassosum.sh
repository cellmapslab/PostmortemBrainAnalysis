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
# 71,880 cases and 383,378 controls: N = 455,258 samples
mkdir -p ${outpath_lassosum}/AD-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat BETA --pvalue P --chr CHR --bp BP --size 455258 --beta yes --thread 3 --out ${outpath_lassosum}/AD-2019 --gwas AD-2019 &> ${outpath_lassosum}/AD-2019/Lassosum.log
wait
echo "Done AD-2019"

### ADHD-2019 ###
# 20,183 cases and 35,191 controls: N = 55,374
mkdir -p ${outpath_lassosum}/ADHD-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ADHD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 55374 --beta no --thread 3 --out ${outpath_lassosum}/ADHD-2019 --gwas ADHD-2019 &> ${outpath_lassosum}/ADHD-2019/Lassosum.log
wait
echo "Done ADHD-2019"

### ALCDEP-2018 ###
# 11,476 cases and 23,080 controls: N = 34,556 samples
mkdir -p ${outpath_lassosum}/ALCDEP-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ALCDEP-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat BETA --pvalue P --chr CHR --bp BP --size 34556 --beta yes --thread 3 --out ${outpath_lassosum}/ALCDEP-2018 --gwas ALCDEP-2018 &> ${outpath_lassosum}/ALCDEP-2018/Lassosum.log
wait
echo "Done ALCDEP-2018"

### AN-2019 ###
# 16,992 cases and 55,525 controls: N = 72517 samples
mkdir -p ${outpath_lassosum}/AN-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AN-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 REF --A2 ALT --stat BETA --pvalue PVAL --chr CHROM --bp POS --size 72517 --beta yes --thread 3 --out ${outpath_lassosum}/AN-2019 --gwas AN-2019 &> ${outpath_lassosum}/AN-2019/Lassosum.log
wait
echo "Done AN-2019"

### ASD-2019 ###
# 18,382 cases and 27,969 controls: N = 46,351 samples
mkdir -p ${outpath_lassosum}/ASD-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ASD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 46351 --beta no --thread 3 --out ${outpath_lassosum}/ASD-2019 --gwas ASD-2019 &> ${outpath_lassosum}/ASD-2019/Lassosum.log
wait
echo "Done ASD-2019"

### AUDIT-2019 ###
# 121,604 (UK Biobank only, no 23AndMe)
mkdir -p ${outpath_lassosum}/AUDIT-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/AUDIT-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 a_0 --A2 a_1 --stat beta_T --pvalue p_T --chr chr --bp bp --size 121604 --beta yes --thread 3 --out ${outpath_lassosum}/AUDIT-2019 --gwas AUDIT-2019 &> ${outpath_lassosum}/AUDIT-2019/Lassosum.log
wait
echo "Done AUDIT-2019"

### BIP-2019 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
mkdir -p ${outpath_lassosum}/BIP-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIP-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 51710 --beta no --thread 3 --out ${outpath_lassosum}/BIP-2019 --gwas BIP-2019 &> ${outpath_lassosum}/BIP-2019/Lassosum.log
wait
echo "Done BIP-2019"

### BIPvsSCZ-2018 ###
# 15,270 BIP cases and 23,585 SCZ cases: N = 38,855 samples
mkdir -p ${outpath_lassosum}/BIPvsSCZ-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIPvsSCZ-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 38855 --beta no --thread 3 --out ${outpath_lassosum}/BIPvsSCZ-2018 --gwas BIPvsSCZ-2018 &> ${outpath_lassosum}/BIPvsSCZ-2018/Lassosum.log
wait
echo "Done BIPvsSCZ-2018"

### BIPSCZvsCONT-2018 ###
# 20,129 BIP cases combined with 33,426 schizophrenia cases and 54,065 controls: N = 107,620 samples
mkdir -p ${outpath_lassosum}/BIPSCZvsCONT-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BIPSCZvsCONT-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 107620 --beta no --thread 3 --out ${outpath_lassosum}/BIPSCZvsCONT-2018 --gwas BIPSCZvsCONT-2018 &> ${outpath_lassosum}/BIPSCZvsCONT-2018/Lassosum.log
wait
echo "Done BIPSCZvsCONT-2018"

### BMI-2015 (ENGAGE) ###
# N = 87,048 samples
mkdir -p ${outpath_lassosum}/BMI-2015_ENGAGE  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BMI-2015_ENGAGE.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 effect_allele --A2 other_allele --stat effect --pvalue P --chr chromosome --bp position --size 87048 --beta yes --thread 3 --out ${outpath_lassosum}/BMI-2015_ENGAGE --gwas BMI-2015_ENGAGE &> ${outpath_lassosum}/BMI-2015_ENGAGE/Lassosum.log
wait
echo "Done BMI-2015 (ENGAGE)"

### BMI-2015 (GIANT) ###
# N = 322,154 samples
mkdir -p ${outpath_lassosum}/BMI-2015_GIANT  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BMI-2015_GIANT.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat b --pvalue p --chr Chromosome --bp Position --size 322154 --beta yes --thread 3 --out ${outpath_lassosum}/BMI-2015_GIANT --gwas BMI-2015_GIANT &> ${outpath_lassosum}/BMI-2015_GIANT/Lassosum.log
wait
echo "Done BMI-2015 (GIANT)"

### BMI-2018 ###
# N = 681,275 samples
mkdir -p ${outpath_lassosum}/BMI-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BMI-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 Tested_Allele --A2 Other_Allele --stat BETA --pvalue P --chr CHR --bp POS --size 681275 --beta yes --thread 3 --out ${outpath_lassosum}/BMI-2018 --gwas BMI-2018 &> ${outpath_lassosum}/BMI-2018/Lassosum.log
wait
echo "Done BMI-2018"

### BMI-2019 ###
# N = 806,834 samples
mkdir -p ${outpath_lassosum}/BMI-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/BMI-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 Testes_Allele --A2 Other_Allele --stat BETA --pvalue P --chr CHR --bp POS --size 806834 --beta yes --thread 3 --out ${outpath_lassosum}/BMI-2019 --gwas BMI-2019 &> ${outpath_lassosum}/BMI-2019/Lassosum.log
wait
echo "Done BMI-2019"

### CAD-2017 ###
# 10,801 cases (incl. angina) and 137,914 controls: N = 148,715 samples
mkdir -p ${outpath_lassosum}/CAD-2017  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CAD-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 effect_allele --A2 noneffect_allele --stat OR --pvalue p-value_gc --chr chr --bp bp_hg19 --size 148715 --beta no --thread 3 --out ${outpath_lassosum}/CAD-2017 --gwas CAD-2017 &> ${outpath_lassosum}/CAD-2017/Lassosum.log
wait
echo "Done CAD-2017"

### Cannabis-2018 ###
# 35,297 ICC and 22,683 23andMe and 126,785 UK-Biobank: N = 184,765 samples
mkdir -p ${outpath_lassosum}/Cannabis-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Cannabis-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat BETA --pvalue P --chr CHR --bp BP --size 184765 --beta yes --thread 3 --out ${outpath_lassosum}/Cannabis-2018 --gwas Cannabis-2018 &> ${outpath_lassosum}/Cannabis-2018/Lassosum.log
wait
echo "Done Cannabis-2018"

### CDG1-2013 ###
# 33,332 cases and 27,888 controls: N = 61,220 samples
mkdir -p ${outpath_lassosum}/CDG1-2013  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CDG1-2013.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 a1 --A2 a2 --stat or --pvalue pval --chr hg18chr --bp bp --size 61220 --beta no --thread 3 --out ${outpath_lassosum}/CDG1-2013 --gwas CDG1-2013 &> ${outpath_lassosum}/CDG1-2013/Lassosum.log
wait
echo "Done CDG1-2013"

### CDG2-2018 ###
# N = 438,997 samples
mkdir -p ${outpath_lassosum}/CDG2-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/CDG2-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 REF --A2 ALT --stat BETA --pvalue PVAL --chr CHROM --bp POS --size 438997 --beta yes --thread 3 --out ${outpath_lassosum}/CDG2-2018 --gwas CDG2-2018 &> ${outpath_lassosum}/CDG2-2018/Lassosum.log
wait
echo "Done CDG2-2018"

### EA-2018 (Educational attainment or IQ) ###
# N = 766,345 samples
mkdir -p ${outpath_lassosum}/EA-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/EA-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat Beta --pvalue Pval --chr CHR --bp POS --size 766345 --beta yes --thread 3 --out ${outpath_lassosum}/EA-2018 --gwas EA-2018 &> ${outpath_lassosum}/EA-2018/Lassosum.log
wait
echo "Done EA-2018"

### ED-2017 ###
# 3,495 cases and 10,982 controls: N = 14,477 samples
mkdir -p ${outpath_lassosum}/ED-2017  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/ED-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 14477 --beta no --thread 3 --out ${outpath_lassosum}/ED-2017 --gwas ED-2017 &> ${outpath_lassosum}/ED-2017/Lassosum.log
wait
echo "Done ED-2017"

### Height-2014 ###
# N = 253,288 samples
mkdir -p ${outpath_lassosum}/Height-2014  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Height-2014.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 Allele1 --A2 Allele2 --stat b --pvalue p --chr Chromosome --bp Position --size 253288 --beta yes --thread 3 --out ${outpath_lassosum}/Height-2014 --gwas Height-2014 &> ${outpath_lassosum}/Height-2014/Lassosum.log
wait
echo "Done Height-2014"

### Intelligence-2018 ###
# N = 269,867 samples
mkdir -p ${outpath_lassosum}/Intelligence-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/Intelligence-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat stdBeta --pvalue P --chr CHR --bp POS --size 269867 --beta yes --thread 3 --out ${outpath_lassosum}/Intelligence-2018 --gwas Intelligence-2018 &> ${outpath_lassosum}/Intelligence-2018/Lassosum.log
wait
echo "Done Intelligence-2018"

### MDD1-2012 ###
# 9,240 cases and 9,519 controls: N = 18,759 samples
mkdir -p ${outpath_lassosum}/MDD1-2012  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD1-2012.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 a1 --A2 a2 --stat or --pvalue pval --chr hg18chr --bp bp --size 18759 --beta no --thread 3 --out ${outpath_lassosum}/MDD1-2012 --gwas MDD1-2012 &> ${outpath_lassosum}/MDD1-2012/Lassosum.log
wait
echo "Done MDD1-2012"

### MDD2-2018 ###
# 135,458 cases and 344,901 control: N = 480,359 samples
mkdir -p ${outpath_lassosum}/MDD2-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD2-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 480359 --beta no --thread 3 --out ${outpath_lassosum}/MDD2-2018 --gwas MDD2-2018 &> ${outpath_lassosum}/MDD2-2018/Lassosum.log
wait
echo "Done MDD2-2018"

### MDD3-2019 ###
# 45,396 cases and 97,250 controls: N = 142,646 samples
mkdir -p ${outpath_lassosum}/MDD3-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/MDD3-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 142646 --beta no --thread 3 --out ${outpath_lassosum}/MDD3-2019 --gwas MDD3-2019 &> ${outpath_lassosum}/MDD3-2019/Lassosum.log
wait
echo "Done MDD3-2019"

### OCD-2017 ###
# 2,688 cases and 7,037 controls: N = 9,725 samples
mkdir -p ${outpath_lassosum}/OCD-2017  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/OCD-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 9725 --beta no --thread 3 --out ${outpath_lassosum}/OCD-2017 --gwas OCD-2017 &> ${outpath_lassosum}/OCD-2017/Lassosum.log
wait
echo "Done OCD-2017"

### PTSD-2019 ###
# 23,212 PTSD cases and 151,447 controls: N = 174,659 samples (only European samples)
mkdir -p ${outpath_lassosum}/PTSD-2019  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/PTSD-2019.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 174659 --beta no --thread 3 --out ${outpath_lassosum}/PTSD-2019 --gwas PTSD-2019 &> ${outpath_lassosum}/PTSD-2019/Lassosum.log
wait
echo "Done PTSD-2019"

### SCZ2-2014 ### 
# 36,989 SCZ cases and 113,075 controls: N 150,064 samples
mkdir -p ${outpath_lassosum}/SCZ2-2014  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/SCZ2-2014.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 a1 --A2 a2 --stat or --pvalue p --chr hg19chrc --bp bp --size 150064 --beta no --thread 3 --out ${outpath_lassosum}/SCZ2-2014 --gwas SCZ2-2014 &> ${outpath_lassosum}/SCZ2-2014/Lassosum.log
wait
echo "Done SCZ2-2014"

### T2D-2017 ###
# 26,676 T2D cases and 132,532 controls: N = 159,208 samples
mkdir -p ${outpath_lassosum}/T2D-2017  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/T2D-2017.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 Allele1 --A2 Allele2 --stat Effect --pvalue P-value --chr Chr --bp Position --size 159208 --beta yes --thread 3 --out ${outpath_lassosum}/T2D-2017 --gwas T2D-2017 &> ${outpath_lassosum}/T2D-2017/Lassosum.log
echo "Done T2D-2017"

### TS-2018 ###
# 4,819 TS cases and 9,488 controls: N = 14,307 samples
mkdir -p ${outpath_lassosum}/TS-2018  # generate output path if not existing already
Rscript --slave ${script_path}/05a_lassosum.R --base ${sumstats_path}/TS-2018.sumstats --target ${target} --ld ${ld_reference} --pheno ${phenotypes} --cov ${covariates} --A1 A1 --A2 A2 --stat OR --pvalue P --chr CHR --bp BP --size 14307 --beta no --thread 3 --out ${outpath_lassosum}/TS-2018 --gwas TS-2018 &> ${outpath_lassosum}/TS-2018/Lassosum.log
wait
echo "Done TS-2018"

