#!/bin/bash

read_genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes"
tools_path="/storage/groups/ccm01/tools"

chromosome=($(seq 1 1 22))

##############################
### Prepare Genotype Files ###
##############################
## 1.1 Here, we will perform three different steps in one qctool v1.5 run:
##	1) Exclude samples with ZNP or MD in the IID (we assume that they come from other projects) 
## 	2) Filter for MAF > 1% (0.01) + HWE < 6 (-log10(1e-06))
## 	3) Remove NA in first column using the -omit-chromosome parameter
#for chr in ${chromosome[@]} ; do # go through all genotype files
#    $tools_path/qctool_v1.5/qctool -g ${read_genotypes}/pmb_info06_chr${chr}.gen.gz -s ${read_genotypes}/pmb_info06_chr${chr}.sample -hwe 6 -maf 0.05 1 -excl-samples ${out_genotypes}/exclude_samples.txt -filetype gen -omit-chromosome -og ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz -os ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.sample &>> ${out_genotypes}/1.1_MAF_HWE_filtered/1.1_MAF_0.05_HWE_6_filtered.log
#done
#echo "1. MAF & HWE filtering done"

## 1.2 Merge all chr files together to generate one .gen.gz file and copy the sample file of chr1
#command=""
#for chr in ${chromosome[@]} ; do
#	  command="${command}${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz "
#done
#cat $command > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz  # run cat
#cp ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr1_maf0.05.sample ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample  # copy sample file and change its name
#echo "2. Merging done"

## 1.3 Generate .bed/.bim/.fam files out of the merged chr files
#$tools_path/PLINK/plink --gen ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz --sample ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample --allow-no-sex --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged
#echo "3. Transformation done"

# 1.4 remove duplicates (use plink version 1 as it does not work for plink version 2 right now...)
cut -f 2 ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim | sort | uniq -d > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_IDstoRemove.txt
$tools_path/plink-1.07-x86_64/plink --bfile ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --noweb --exclude ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_IDstoRemove.txt --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup

# 1.5 we still have duplicated SNPs as PLINK just checks for two replicates not multiple ones
cut -f 2 ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup.bim | sort | uniq -d > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_IDstoRemove.txt 
$tools_path/plink-1.07-x86_64/plink --bfile ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup --noweb --exclude ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_IDstoRemove.txt --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final

# 1.6 AND again remove duplicates (as we have SNPs which can be found 4 times)
cut -f 2 ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final.bim | sort | uniq -d > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final_IDstoRemove.txt 
$tools_path/plink-1.07-x86_64/plink --bfile ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --noweb --exclude ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final_IDstoRemove.txt --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final2

# 1.7 last round, as SNP rs71074666 has 5 different allels
cut -f 2 ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final2.bim | sort | uniq -d > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final2_IDstoRemove.txt 
$tools_path/plink-1.07-x86_64/plink --bfile ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final2 --noweb --exclude ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final2_IDstoRemove.txt --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final

# PLINK v1.07 => important for removing duplicated (multilallelic SNPs) as this version of 
# PLINK removes all except for the first occassion of the multiallelic SNP. PLINK v1.90 
# removes ALL and we need one appearance of each SNP for the clumping as PLINK otherwise 
# gives a lot of warnings and can not do LD analysis for thousands of SNPs.

# The filtering for the 169 samples with gene expression data is done during the input generation step!

