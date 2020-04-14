#!/bin/bash

# Generate variable including the help text for the script
help="
This bash script is used to do the preprocessing of the cohort genotypes.

The script runs only on bash!!! (not shell "sh")

Steps done by running this script:
1. Run qctool v1.5 to:
   - filter for MAF > 1% (0.01) + HWE < 6 (-log10(1e-06)), 
   - remove NA in first column using the -omit-chromosome parameter
2. Change SNP IDs to Chr:Pos and remove multiallelic SNPs
3. Merge all chr files together to generate one .gen.gz file and copy the sample file of chr1
4. Use PLINK to generate .bed/.bim/.fam files out of the merged chr files
5. Remove duplicated SNPs
6. Use PLINK to filter genotype data for samples with gene expression data (169 samples in total)

Usage: bash 01b_prepare_genotypes.sh -c <path/to/cohort-genotype-files> -o <path/to/store/output-files> -t <path/to/tools> -s <path/to/scripts>

-c  --cohort      path to the cohort genotype (chromosome based) files
-o  --out         path to store the output files into
-t  --tools       path to the neccessary tools
-s  --scripts     path to the script (code) files

optional:
-h  --help        show this help text
"


# Define empty input variables
cohort_path=  # path to genotype files
out_path=  # directory, to store the output files into, here: /storage/groups/groups_ccm/workspace/MPIP-Dean_postmortem_brain/01.GWAS/01_preprocessed_genotypes"
tools_path=  # path to the tools needed for the script 
scripts_path=  # path to the script (code) files

# Function which prints the help variable
opt_h() {
    echo "$help"
}


# Parse command line parameters and save them into the defined variables from above
while [ "$1" != "" ]; do
    case $1 in
        -c | --cohort )     shift
                            cohort_path=$1
                            ;;
        -o | --out )        shift
                            out_path=$1
                            ;;
        -t | --tools )      shift
                            tools_path=$1
                            ;;
        -s | --scripts )    shift
                            scripts_path=$1
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
if [ $cohort_path ] && [ $out_path ] && [ $tools_path ] && [ $scripts_path ]; then

chromosome=($(seq 1 1 22))

    # Generate the output directory (if not existing)
    maf_hwe_out_path="${out_path}/2.1_MAF_HWE_filtered"  # output path for maf and hwe filtered data
    recoded_out_path="${out_path}/2.2_changed_IDs"  # output path for recoded and multiallelic filtered SNP IDs
    merged_out_path="${out_path}/2.3_merged_chr"  # output path for merged files
    removed_dups_out_path="${out_path}/2.4_dups_removed"  # output path for removed duplicates
    mkdir -p $maf_hwe_out_path $recoded_out_path $merged_out_path $removed_dups_out_path

    #####################################
    ### Prepare Cohort Genotype Files ###
    #####################################
    # 1. Here, we will perform three different steps in one qctool v1.5 run:
    # 	1) Filter for MAF > 1% (0.01) + HWE < 6 (-log10(1e-06))
    # 	2) Remove NA in first column using the -omit-chromosome parameter
    for chr in ${chromosome[@]} ; do # go through all genotype files
        $tools_path/qctool_v1.5/qctool -g ${cohort_path}/03_updateMap_${chr}_hwe.out -s ${cohort_path}/03_updateMap_hwe.sample -hwe 6 -maf 0.01 1 -filetype gen -omit-chromosome -og ${maf_hwe_out_path}/03_updateMap_${chr}_maf0.01_hwe.out -os ${maf_hwe_out_path}/03_updateMap_${chr}_maf0.01_hwe.sample &>> ${maf_hwe_out_path}/1.1_MAF_0.01_HWE_6_filtered.log
    done
    echo "1. MAF & HWE filtering done"

    # 2. Reorder: recode SNP ID and remove Multiallelic SNPs
    for chr in ${chromosome[@]} ; do
        Rscript ${scripts_path}/02b_prepare_cohorts.R ${maf_hwe_out_path}/03_updateMap_${chr}_maf0.01_hwe ${recoded_out_path}/03_updateMap_${chr}_maf0.01_hwe ${tools_path} &>> ${recoded_out_path}/1.2_changed_IDs_MAF0.01.log
        rm ${recoded_out_path}/03_updateMap_${chr}_maf0.01_hwe_renamed.*
        gzip ${recoded_out_path}/03_updateMap_${chr}_maf0.01_hwe_renamed_noMultiAllelicSNPs.gen
    done
    echo "2. Transformation and Filtering done"

    # 3. Concatenate the .gen.gz files (with MAF 0.01 and the changed SNP ID format)
    command=""
    for chr in ${chromosome[@]} ; do
        command="${command}${recoded_out_path}/03_updateMap_${chr}_maf0.01_hwe_renamed_noMultiAllelicSNPs.gen.gz "
    done
    cat $command > ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged.gen.gz  # run cat
    cp ${maf_hwe_out_path}/03_updateMap_1_maf0.01_hwe.sample ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged.sample  # copy sample file
    echo "3. Merging done"

    # 4. convert gen.gz + sample to .bed/.bim/.fam files
    $tools_path/PLINK/plink --gen ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged.gen.gz --sample ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged.sample --allow-no-sex --make-bed --out ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged
    echo "4. Convertion to PLINK format done"

    # 5. remove duplicates (use plink version 1 as it does not work for plink version 2 right now...)
    $tools_path/PLINK/plink --bfile ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged --list-duplicate-vars suppress-first --out ${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged  # *.dupvar file is empty (NO duplicates to remove!)
    $tools_path/PLINK/plink --bfile ${merged_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged --exclude ${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged.dupvar --make-bed  --out /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup
    cut -f 2 /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup.bim | sort | uniq -d > /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_IDstoRemove.txt
    $tools_path/PLINK/plink --bfile /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup --exclude /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_IDstoRemove.txt --make-bed --out /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final
    echo "5. Removal of duplicates done"

    # 6. change .bed/.bim/.fam back to .gen/.sample
    $tools_path/PLINK/plink --bfile /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final --recode oxford --out /${removed_dups_out_path}/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_merged_noDup_final
    echo "6. Convertion back to oxford format done"

    echo "Preprocess cohort genotype data done"

else
    # If one of the parameters is empty or not specified print the help text
    opt_h
fi

