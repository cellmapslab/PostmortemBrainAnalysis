#!/bin/bash

# Generate variable including the help text for the script
help="
This bash script is used to calculate polygenic risk scores (PRS) using LDpred.

LDpred needs to be installed beforehand using: pip install --user ldpred
Check if LDpred installation works, running: ldpred --help

The script runs only on bash!!! (not shell "sh")

1. Step: run ldpred coord to generate a coordinated data file

2. Step: run ldpred gibbs to estimate Gibbs sampler (LDpred) SNP weights for the defined causal fractions

3. Step: run ldpred score to generate the polygenic risk scores (PRS) for the samples

4. Step run ldpred p+t to generate PRS using pruning and thresholding instead of LDpreds Gibbs sampling method

Usage: bash 06a_LDpred.sh -g <path/to/target-geno-files> -t <path/to/phenotypes.txt> -c <path/to/cohort-geno-files> -s <path/to/sumstats> -n <number-of-individuals> -r <rsID-colname> -a1 <ref-allele-colname> -a2 <alt-allele-colname> -i <position-colname> -m <chromosome-colname> -p <pvalue-colname> -e <effect-colname> -o <path/for/output-files> 

-g  --genotypes   path and file name prefix (without file ending) of the plink formatted genotype files (e.g. pmb_info06_all_chr_merged_scz2_filtered)
-t  --phenotypes  path and file name of the phenotypes file
-v  --covariates  path and file name of the covariates file
-c  --cohort      path and file name prefix (without file ending) of the plnk formatted cohort genotype files (e.g. 03_updateMap_all_chr_merged_scz2_filtered)
-s  --sumstats    path and file name of the GWAS summary statistics file (e.g. MDD2_2018_filtered.sumstats) 
-n  --nums        number of individuals (sample size) used to generate the GWAS summary statistics (e.g. 20,000 cases and 30,000 controls: N = 50,000)
-r  --rsid        name of rs ID column in the GWAS summary statistics file
-a1  --ref        name of reference allele (usually A1) column in the GWAS summary statistics file
-a2  --alt        name of alternate allele (usually A2) column in the GWAS summary statistics file
-i  --pos         name of position (base pairs = bp) column in the GWAS summary statistics file
-m  --chr         name of chromosome column in the GWAS summary statistics file
-p  --pval        name of p-value column in the GWAS summary statistics file
-e  --eff         name of effect column (either OR/odds ratio or beta) in the GWAS summary statistics file
-o  --out         path to store the output files into

optional:
-b  --beta        check, if the effect is beta, then the parameters --se and --z-from-se are set as well
-se --stderr      if --beta is set, you need to specifiy the name of Standard Error (usually SE) column in the GWAS summary statistics file
-r  --radius      LD radius to use, e.g. number of SNPs divided by 3000 results in a radius approximately to 2000 (2Mb) (default: (lines of sumstats - 1) / 3000)
-h  --help        show the help text
"

# Define empty input variables
genotypes=  # path and file name prefix of the plink formated genotype files
phenotypes=  # path and file name of the phenotypes file
covariates=  # path and file name of the covariates file
cohort=  # path and file name prefix of the plink formatted cohort files
sumstats=  # path and file name of the GWAS summary statistics file
nums=  # sample size (number of individuals) used to calculate the GWAS summary statistics on
rsid=  # name of rs ID column in the GWAS summary statistics file
ref=  # name of reference allele (usually A1) column in the GWAS summary statistics file
alt=  # name of alternate allele (usually A2)  column in the GWAS summary statistics file
pos=  # name of position (base pairs = bp) column in the GWAS summary statistics file
chr=  # name of chromosome column in the GWAS summary statistics file
pval=  # name of p-value column in the GWAS summary statistics file
effect=  # name of effect column (either OR/odds ratio or beta) in the GWAS summary statistics file
out_path=  # directory, to store the output files into, here: /storage/groups/groups_ccm/workspace/MPIP-Dean_postmortem_brain/genotype_data/3_LDpred"
beta=  # check if effect is beta, if yes: calculate Z-score from beta and standard error
se=  # standard error
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
        -t | --phenotypes ) shift
                            phenotypes=$1
                            ;;
        -v | --covariates ) shift
                            covariates=$1
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
        -r | --rsid )       shift
                            rsid=$1
                            ;;
        -a1 | --ref )       shift
                            ref=$1
                            ;;
        -a2 | --alt )       shift
                            alt=$1
                            ;;
        -i | --pos )        shift
                            pos=$1
                            ;;
        -m | --chr )        shift
                            chr=$1
                            ;;
        -p | --pval )       shift
                            pval=$1
                            ;;
        -e | --eff )        shift
                            effect=$1
                            ;;
        -o | --out )        shift
                            out_path=$1
                            ;;
        -b | --beta )       shift
                            beta=$1
                            ;;
        -se | --stderr )    shift
                            se=$1
                            ;;
        -d | --radius )     shift
                            radius=$1
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
if [ $genotypes ] && [ $phenotypes ] && [ $covariates ] && [ $cohort ] && [ $sumstats ] && [ $nums ] && [ $rsid ] && [ $ref ] && [ $alt ] && [ $pos ] && [ $chr ] && [ $pval ] && [ $effect ] && [ $out_path ] ; then

    # Print the server name on which the job is running on the grid
    echo ">>> Start pipeline on $HOSTNAME <<<"

    # Generate output folders (if not existing)
    coord_path="${out_path}/01_coord"  # output path for coordinated data file
    gibbs_path="${out_path}/02_gibbs"  # output path for different Gibbs sampler (LDpred) SNP weight files
    validate_path="${out_path}/03_score"  # output path for calculated polygenic risk scores (PRS)
    pt_path="${out_path}/04_P+T"  # output path for pruning and thresholding results
    mkdir -p $coord_path $gibbs_path $validate_path $pt_path

    # First, we need the number of SNPs
    if [[ -z "$radius" ]]; then
        # set radius to 1000 (as ldpred coord removes more SNPs...)
        radius=1000
        # set radius as number of SNPs divided by 3000 (recommended in LDpred paper!)
        # radius=$(wc -l < $sumstats)
        # radius=$(((radius - 1) / 3000))
    fi
    echo "Radius for ${sumstats} set to: $radius"

    # Second, we need the genotype file name prefix (no file ending) without the path 
    out_file_name=${genotypes##*/}
    
    # Third, we need to remove the output file of coord_genotypes.py, if it already exists (does not work otherwise...)
    if [[ -e ${coord_path}/${out_file_name}.coord_file ]]; then
        rm ${coord_path}/${out_file_name}.coord_file
    fi

    # 1. Run coord_genotypes.py
    if [[ "$beta" == "yes" ]] && [[ -n $se ]] ; then
        ldpred coord --gf $cohort --ssf ${sumstats} --N ${nums} --rs ${rsid} --A1 ${ref} --A2 ${alt} --pos ${pos} --chr ${chr} --pval ${pval} --eff ${effect} --se ${se} --z-from-se --maf 0 --vgf ${genotypes} --ssf-format CUSTOM --out ${coord_path}/${out_file_name}.coord_file &> ${coord_path}/coord_genotypes.log
    else
        ldpred coord --gf $cohort --ssf ${sumstats} --N ${nums} --rs ${rsid} --A1 ${ref} --A2 ${alt} --pos ${pos} --chr ${chr} --pval ${pval} --eff ${effect} --maf 0 --vgf ${genotypes} --ssf-format CUSTOM --out ${coord_path}/${out_file_name}.coord_file &> ${coord_path}/coord_genotypes.log
    fi
    echo "coord done"  # --ssf-format $format

    # 2. Run LDpred: test for 1 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.03 0.01 0.003 0.001 (takes forever...)
    ldpred gibbs --cf ${coord_path}/${out_file_name}.coord_file --ldr ${radius} --f 1 0.5 0.3 0.1 0.03 0.01 0.003 0.001 --ldf ${gibbs_path}/${out_file_name} --N ${nums} --out ${gibbs_path}/${out_file_name} &> ${gibbs_path}/LDpred.log
    echo "gibbs done"

    # 3. Run scoring (get score for each patient)
    # 3.1 Once without covariates
#    ldpred score --gf ${genotypes} --rf ${gibbs_path}/${out_file_name} --f 1 0.5 0.3 0.1 0.03 0.01 0.003 0.001 --pf ${phenotypes} --pf-format LSTANDARD --out ${out_path}/03_validate/${out_file_name}_PRS_without_cov &> ${out_path}/03_validate/validate.log
#    echo "score without covariates done"
    # 3.2 Once with covariates
#    ldpred score --gf ${genotypes} --rf ${gibbs_path}/${out_file_name} --f 1 0.5 0.3 0.1 0.03 0.01 0.003 0.001 --pf ${phenotypes} --pf-format LSTANDARD --cov-file ${covariates} --out ${out_path}/03_validate/${out_file_name}_PRS_with_cov &> ${out_path}/03_validate/validate.log
#    echo "score with covariates done"

    # 4. run LDpred for pruning and thresholding (P+T)
#    ldpred p+t --cf ${coord_path}/${out_file_name}.coord_file --ldr ${radius} --p 1 0.05 0.01 0.00000005 --out ${pt_path}/${out_file_name} &> ${pt_path}/LDpred_p+t.log
#    echo "gibbs for P+T done"

    # 5. run scoring (get score for each patient) for the P+T
    # 5.1 once without cov
#    ldpred score --gf ${genotypes} --rf ${pt_path}/${out_file_name} --p 1 0.05 0.01 0.00000005 --pf ${phenotypes} --pf-format LSTANDARD --out ${pt_path}/${out_file_name}_PRS_without_cov &> ${pt_path}/validate_p+t.log
#    echo "score for P+T without covariates done"
    # 5.2 Once with cov
#    ldpred score --gf ${genotypes} --rf ${pt_path}/${out_file_name} --p 1 0.05 0.01 0.00000005 --pf ${phenotypes} --pf-format LSTANDARD --cov-file ${covariates} --out ${pt_path}/${out_file_name}_PRS_with_cov &> ${pt_path}/validate_p+t.log
#    echo "score for P+T with covariates done"

    # Print that the job is done
    echo "!!!Done with LDpred of ${out_file_name}!!!"

# If one of the parameters is empty or not specified print the help text
else
    opt_h
fi

