#!/bin/bash

read_genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_genotypes"
tools_path="/storage/groups/ccm01/tools"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"

chromosome=($(seq 1 1 22))

##############################
### Prepare Genotype Files ###
##############################
# 1. Here, we will perform three different steps in one qctool v1.5 run:
#	1) Exclude samples with ZNP or MD in the IID (we assume that they come from other projects) 
# 	2) Filter for MAF > 1% (0.01) + HWE < 6 (-log10(1e-06))
# 	3) Remove NA in first column using the -omit-chromosome parameter
for chr in ${chromosome[@]} ; do # go through all genotype files
    $tools_path/qctool_v1.5/qctool -g ${read_genotypes}/pmb_info06_chr${chr}.gen.gz -s ${read_genotypes}/pmb_info06_chr${chr}.sample -hwe 6 -maf 0.01 1 -excl-samples ${out_genotypes}/exclude_samples.txt -filetype gen -omit-chromosome -og ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.01.gen.gz -os ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.01.sample &>> ${out_genotypes}/1.1_MAF_HWE_filtered/1.1_MAF_0.01_HWE_6_filtered.log
done
echo "1. MAF & HWE filtering done"

# 1.2 Reorder: recode SNP ID and remove Multiallelic SNPs
for chr in ${chromosome[@]} ; do
	Rscript ${script_path}/01_sub_prepare_genotypes.R ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.01 ${out_genotypes}/1.2_changed_IDs/pmb_info06_chr${chr}_maf0.01 ${tools_path} &>> ${out_genotypes}/1.2_changed_IDs/1.2_changed_IDs_MAF0.01.log
	rm ${out_genotypes}/1.2_changed_IDs/pmb_info06_chr${chr}_maf0.01_renamed.*
	gzip ${out_genotypes}/1.2_changed_IDs/pmb_info06_chr${chr}_maf0.01_renamed_noMultiAllelicSNPs.gen
done

# 1.3 Concatenate the .gen.gz files (once with MAF 0.05 and rs SNP IDs, once with MAF 0.01 and the changed SNP ID format)
command=""
for chr in ${chromosome[@]} ; do
	  command="${command}${out_genotypes}/1.2_changed_IDs/pmb_info06_chr${chr}_maf0.01_renamed_noMultiAllelicSNPs.gen.gz "
done
cat $command > ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged.gen.gz  # run cat
cp ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr1_maf0.01.sample ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged.sample  # copy sample file

# 1.4 convert gen.gz + sample to .bed/.bim/.fam files
$tools_path/PLINK/plink --gen ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged.gen.gz --sample ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged.sample --allow-no-sex --make-bed --out ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged

# 1.5 remove duplicates (use plink version 1 as it does not work for plink version 2 right now...)
$tools_path/PLINK/plink --bfile ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged --list-duplicate-vars suppress-first --out ${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged  # *.dupvar file is empty (NO duplicates to remove!)
$tools_path/PLINK/plink --bfile ${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_maf0.01_IDs_changed_merged --exclude ${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged.dupvar --make-bed  --out /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup

# 1.6 Check if dups are really removed (again)
cut -f 2 /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup.bim | sort | uniq -d > /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_IDstoRemove.txt
$tools_path/PLINK/plink --bfile /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup --exclude /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_IDstoRemove.txt --make-bed --out /${out_genotypes}/1.4_dups_removed/pmb_info06_all_chr_maf0.01_IDs_changed_merged_noDup_final

echo "Transform genotype data done"

