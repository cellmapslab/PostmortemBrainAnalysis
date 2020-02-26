#!/bin/bash

# Generate variable including the help text for the script
help="
This bash script is used to do the preprocessing of the genotypes.

The script runs only on bash!!! (not shell "sh")

Steps done by running this script:
1. Run qctool v1.5 to:
   - exclude samples with ZNP or MD in the IID (we assume that they come from other projects), 
   - filter for MAF > 5% (0.05) + HWE < 6 (-log10(1e-06)), 
   - remove NA in first column using the -omit-chromosome parameter
2. Merge all chr files together to generate one .gen.gz file and copy the sample file of chr1
3. Use PLINK to filter genotype data for samples with gene expression data (169 samples in total) and 
   generate .bed/.bim/.fam files out of the merged chr files
4. Use PLINK (v1.07) to remove duplicated SNPs and change back to Oxford file format (.gen.gz/.sample)

The last step (remove duplicates) needs to run multiple times as it just removes one duplicated SNP.
Several SNPs have multiple duplications included (e.g. SNP rs71074666 has 5 different allels). 

PLINK v1.07 is important for removing duplicated (multilallelic SNPs) as this version of PLINK removes 
all except for the first occassion of the multiallelic SNP. PLINK v1.90 removes ALL and we need one 
appearance of each SNP for the clumping as PLINK otherwise gives a lot of warnings and can not do LD 
analysis for thousands of SNPs.

Usage: bash 01a_prepare_genotypes.sh -g <path/to/genotype-files> -e <exclude-file> -i <include-file> -o <path/to/store/output-files> -t <path/to/tools>

-g  --genotypes   path to the genotype (chromosome based) files after imputation
-e  --exclude     file (name with path) including all sample IDs to exclude for further analyses steps (samples coming from other projects)
-i  --include     file (name with path) including all sample IDs for which expression data is also available and we want to keep for further analyses
-o  --out         path to store the output files into
-t  --tools       path to the neccessary tools

optional:
-h  --help        show this help text
"


# Define empty input variables
genotypes=  # path to genotype files
exclude=  # file with IDs to exclude
include=  # file with IDs to include (samples with expression data)
out_path=  # directory, to store the output files into, here: /storage/groups/groups_ccm/workspace/MPIP-Dean_postmortem_brain/01.GWAS/01_preprocessed_genotypes"
tools_path=  # path to the tools needed for the script 

# Function which prints the help variable
opt_h() {
    echo "$help"
}


# Parse command line parameters and save them into the defined variables from above
while [ "$1" != "" ]; do
    case $1 in
        -g | --genotypes )  shift
                            genotypes=$1
                            ;;
        -e | --exclude )    shift
                            exclude=$1
                            ;;
        -i | --include )    shift
                            include=$1
                            ;;
        -o | --out )        shift
                            out_path=$1
                            ;;
        -t | --tools )      shift
                            tools_path=$1
                            ;;
        -h | --help )       opt_h
                            exit
                            ;;
        * )                 opt_h
                            exit 1
    esac
    shift
done


# Check, if all necessary parameters are given
if [ $genotypes ] && [ $exclude ] && [ $include ] && [ $out_path ] && [ $tools_path ]; then

    chromosome=($(seq 1 1 22))

    ##############################
    ### Prepare Genotype Files ###
    ##############################
    # 1.1 Here, we will perform three different steps in one qctool v1.5 run:
    #	1) Exclude samples with ZNP or MD in the IID (we assume that they come from other projects) 
    # 	2) Filter for MAF > 5% (0.05) + HWE < 6 (-log10(1e-06))
    # 	3) Remove NA in first column using the -omit-chromosome parameter
    for chr in ${chromosome[@]} ; do # go through all genotype files
        $tools_path/qctool_v1.5/qctool -g ${genotypes}/pmb_info06_chr${chr}.gen.gz -s ${genotypes}/pmb_info06_chr${chr}.sample -hwe 6 -maf 0.05 1 -excl-samples ${exclude} -filetype gen -omit-chromosome -og ${out_path}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz -os ${out_path}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.sample &>> ${out_path}/1.1_MAF_HWE_filtered/1.1_MAF_0.05_HWE_6_filtered.log
    done
    echo "1. MAF & HWE filtering done"

    # 1.2 Merge all chr files together to generate one .gen.gz file and copy the sample file of chr1
    command=""
    for chr in ${chromosome[@]} ; do
	    command="${command}${out_path}/1.1_MAF_HWE_filtered/pmb_info06_chr${chr}_maf0.05.gen.gz "
    done
    cat $command > ${out_path}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz  # run cat
    cp ${out_path}/1.1_MAF_HWE_filtered/pmb_info06_chr1_maf0.05.sample ${out_path}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample  # copy sample file and change its name
    echo "2. Merging done"

    # 1.3 Filter genotype data for samples with gene expression data (169 samples) and generate .bed/.bim/.fam files out of the merged chr files
    $tools_path/PLINK/plink --gen ${out_path}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.gen.gz --sample ${out_path}/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.sample --keep ${include} --allow-no-sex --make-bed --out ${out_path}/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged
    echo "3. Transformation and filtering done"

    # Now, we need to remove all duplicated SNPs but keep the first one (needed for Matrix EQTL and the clumping)
    # 1.4 remove duplicates (use plink version 1 as it does not work for plink version 2 right now...)
    cut -f 2 ${out_path}/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged.bim | sort | uniq -d > ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_IDstoRemove.txt
    $tools_path/plink-1.07-x86_64/plink --bfile ${out_path}/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --noweb --exclude ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_IDstoRemove.txt --make-bed --out ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup

    # 1.5 we still have duplicated SNPs as PLINK just checks for two replicates not multiple ones
    cut -f 2 ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.bim | sort | uniq -d > ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup_IDstoRemove.txt
    $tools_path/plink-1.07-x86_64/plink --bfile ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --noweb --exclude ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup_IDstoRemove.txt --make-bed --out ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2

    # 1.6 AND again remove duplicates (as we have SNPs which can be found 4 times)
    cut -f 2 ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2.bim | sort | uniq -d > ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2_IDstoRemove.txt 
    $tools_path/plink-1.07-x86_64/plink --bfile ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2 --noweb --exclude ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2_IDstoRemove.txt --make-bed --out ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3

    # 1.7 last round, as SNP rs71074666 has 5 different allels
    cut -f 2 ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3.bim | sort | uniq -d > ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3_IDstoRemove.txt 
    $tools_path/plink-1.07-x86_64/plink --bfile ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3 --noweb --exclude ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3_IDstoRemove.txt --make-bed --out ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup

    # remove unnecessary files
    rm ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup2.*
    rm ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup3.*

    # 1.8 change .bed/.bim/.fam files back to Oxford Style (.gen.gz/.sample)
    $tools_path/PLINK/plink --bfile ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --recode oxford gen-gz --out ${out_path}/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup

else
    # If one of the parameters is empty or not specified print the help text
    opt_h
fi

