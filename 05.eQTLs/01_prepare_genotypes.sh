#!/bin/bash

read_genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes"
tools_path="/storage/groups/ccm01/tools"

chromosome=($(seq 1 1 22))

##############################
### Prepare Genotype Files ###
##############################
# 1.1 Here, we will perform three different steps in one qctool v1.5 run:
#	1) Exclude samples with ZNP or MD in the IID (we assume that they come from other projects) 
# 	2) Filter for MAF > 1% (0.01) + HWE < 6 (-log10(1e-06))
# 	3) Remove NA in first column using the -omit-chromosome parameter
for chr in ${chromosome[@]} ; do # go through all genotype files
    $tools_path/qctool_v1.5/qctool -g ${read_genotypes}/pmb_info06_chr${chr}.gen.gz -s ${read_genotypes}/pmb_info06_chr${chr}.sample -hwe 6 -maf 0.05 1 -excl-samples ${out_genotypes}/exclude_samples.txt -filetype gen -omit-chromosome -og ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz -os ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.sample &>> ${out_genotypes}/1.1_MAF_HWE_filtered/1.1_MAF_0.05_HWE_6_filtered.log
done
echo "1. MAF & HWE filtering done"

# 1.2 Merge all chr files together to generate one .gen.gz file and copy the sample file of chr1
command=""
for chr in ${chromosome[@]} ; do
	  command="${command}${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz "
done
cat $command > ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz  # run cat
cp ${out_genotypes}/1.1_MAF_HWE_filtered/pmb_info06_chr1_maf0.05.sample ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample  # copy sample file and change its name
echo "2. Merging done"

# 1.3 Generate .bed/.bim/.fam files out of the merged chr files
$tools_path/PLINK/plink --gen ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz --sample ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample --allow-no-sex --make-bed --out ${out_genotypes}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged
echo "3. Transformation done"

# The .bed/.bim/.fam file generation is NOT needed in the eQTL analysis!
# The filtering for the 169 samples with gene expression data is done during the input generation step!

