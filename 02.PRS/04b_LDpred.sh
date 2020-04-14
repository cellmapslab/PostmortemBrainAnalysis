#!/bin/bash

# Generate variable including the help text for the script
help="
This bash script is used to calculate polygenic risk scores (PRS) using LDpred.

The script runs only on bash!!! (not shell "sh")

1. Step: run coord_genotypes.py

2. Step: run LDpred.py

3. Step: run validate.py

Usage: LDpred.sh -g <path/to/get/plink-geno-files> -s <path/to/sumstats> -n <number-of-individuals> -o <path/to/store/output-files> -t <path/to/LDpred>

-g  --genotypes   path and file name prefix (without file ending) of the plink formatted genotype files (e.g. pmb_info06_all_chr_merged_scz2_filtered)
-c  --cohort      path and file name prefix (without file ending) of the plnk formatted cohort genotype files (e.g. 03_updateMap_all_chr_merged_scz2_filtered)
-s  --sumstats    path and file name of the GWAS summary statistics file (e.g. MDD2_2018_filtered.sumstats) 
-n  --num         number of individuals (sample size) used to generate the GWAS summary statistics (e.g. 20,000 cases and 30,000 controls: N = 50,000)
-f  --ssformat    format to use for the summary statistics file in the coord_genotypes.py script (e.g. STANDARD, BASIC, GIANT, GIANT2, DECODE, PGC, PGC_large)
-t  --tools       path to the neccessary tool/s for the pipeline (LDpred)
-o  --out         path to store the output files into

optional:
-r  --radius      LD radius to use, e.g. number of SNPs divided by 3000 results in a radius approximately to 2000 (2Mb) (default: (lines of sumstats - 1) / 3000)
-h  --help        show the help text
"


# Define empty input variables
genotypes=  # path and file name prefix of the plink formated genotype files
cohort=  # path and file name prefix of the plink formatted cohort files
sumstats=  # path and file name of the GWAS summary statistics file
out_path=  # directory, to store the output files into, here: /storage/groups/groups_ccm/workspace/MPIP-Dean_postmortem_brain/genotype_data/3_LDpred"
tools_path=  # path to the tools needed for the script 
nums=  # sample size (number of individuals) used to calculate the GWAS summary statistics
format=  # format to use for the summary statistics file in the coord_genotypes.py script
radius=  # ld radius (usually sample size/number of SNPs divided by 3000 =~ 2000)


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
        -c | --cohort )     shift
                            cohort=$1
                            ;;
        -s | --sumstats )   shift
                            sumstats=$1
                            ;;
        -n | --nums )       shift
                            nums=$1
                            ;;
        -f | --ssformat )   shift
                            format=$1
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
if [ $genotypes ] && [ $cohort ] && [ $sumstats ] && [ $nums ] && [ $out_path ] && [ $tools_path ]; then

    # Print the server name on which the job is running on the grid
    echo ">>> Start pipeline on $HOSTNAME <<<"

    # First, we need the number of SNPs
    if [[ -z "$radius" ]]; then
        # set radius to 1000 (as coord_genotypes.py removes more SNPs...)
        radius=1000
        # set radius as number of SNPs divided by 3000 (recommended in LDpred paper!)
        # radius=$(wc -l < $sumstats)
        # radius=$(((radius - 1) / 3000))
    fi
    #echo "Calculated radius for ${sumstats}: $radius"

    # Second, we need the genotype file name prefix (no file ending) without the path 
    out_file_name=${genotypes##*/}
    #echo $out_file_name
    
    # Third, we need to remove the output file of coord_genotypes.py, if it already exists (does not work otherwise...)
    if [[ -e ${out_path}/01_coord/${out_file_name}.coord_file ]]; then
        rm ${out_path}/01_coord/${out_file_name}.coord_file
    fi

    # 1. Run coord_genotypes.py
    /storage/groups/ccm01/tools/ldpred/ldpred/ldpred coord --gf $cohort --ssf ${sumstats} --N ${nums} --vgf ${genotypes} --ssf-format $format --out ${out_path}/01_coord/${out_file_name}.coord_file &> ${out_path}/01_coord/coord_genotypes.log
    echo "coord done"  # --ssf-format $format

    # 2. Run LDpred: change causal fraction back to PS = 1,0.3,0.1,0.03,0.01,0.003,0.001,0.0003,0.0001
    /storage/groups/ccm01/tools/ldpred/ldpred/ldpred gibbs --cf ${out_path}/01_coord/${out_file_name}.coord_file --ldr ${radius} --f 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.03 0.01 0.003 0.001 --ldf ${out_path}/02_LDpred/${out_file_name} --N ${nums} --out ${out_path}/02_LDpred/${out_file_name} &> ${out_path}/02_LDpred/LDpred.log
    echo "gibbs done"

    # 3. Run scoring (get score for each patient)
    /storage/groups/ccm01/tools/ldpred/ldpred/ldpred score --gf ${genotypes} --rf ${out_path}/02_LDpred/${out_file_name} --f 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.03 0.01 0.003 0.001 --pf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/phenotypes.txt --out ${out_path}/03_validate/${out_file_name}_PRS &> ${out_path}/03_validate/validate.log
    echo "score done"

    # 4. run LDpred for pruning and thresholding (P+T)
    /storage/groups/ccm01/tools/ldpred/ldpred/ldpred p+t --cf ${out_path}/01_coord/${out_file_name}.coord_file --ldr ${radius} --p 1 0.05 0.00000005 --out ${out_path}/04_P+T/${out_file_name} &> ${out_path}/04_P+T/LDpred_p+t.log
    echo "gibbs for P+T done"

    # 5. run scoring (get score for each patient) for the P+T
    /storage/groups/ccm01/tools/ldpred/ldpred/ldpred score --gf ${genotypes} --rf ${out_path}/04_P+T/${out_file_name} --p 1 0.05 0.00000005 --pf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/phenotypes.txt --out ${out_path}/04_P+T/${out_file_name}_PRS &> ${out_path}/04_P+T/validate_p+t.log
    echo "score for P+T done"

    # Print that the job is done
    echo "!!!Done with LDpred of ${out_file_name}!!!"

# If one of the parameters is empty or not specified print the help text
else
    opt_h
fi

