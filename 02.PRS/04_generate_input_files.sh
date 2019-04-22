#!/bin/bash

read_genotypes="/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/2018-12-20_post_mortem_brain_imputation"
out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/genotype_data"
script_path="/storage/groups/ccm01/code/MPIP-Dean_postmortem_brain/PRS"
tools_path="/storage/groups/ccm01/tools"
sum_stats_path="/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/"

chromosome=($(seq 1 1 22))


###############################################################
### 3. Compare SNPs Between Genotype And Summary Statistics ###
###############################################################
# # Best method: run for each summary statistics file on it's own (comment rest out) some errors and problems with all executed after each other...
# # As the merged genotype file does not change, we will refer to its location and name here
geno_in="${out_genotypes}/1.3_merged_chr/pmb_info06_all_chr_IDs_changed_merged_noDup_final"

# # 3.1 SCZ2 
# # ---------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/SCZ2.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/SCZ2/00_input/SCZ2_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/SCZ2/00_input/pmb_info06_all_chr_merged_scz2_filtered"
# snps="${out_genotypes}/3_LDpred/SCZ2/00_input/pmb_info06_SCZ2_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/SCZ2/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # For SCZ2 we need to transform the summary statistic file back into the PGC small format (adding chr to chromsome numbers and change bp back to column 5)
# # 1) We add chr to the chromosome number and save it into a new file
# awk 'NR>1 {$1="chr"$1} 1' OFS="\t" $sumstats_out > ${sumstats_out}2
# # 2) Change order of columns back to: chr, snp, a1, a2, bp, info, or, se, p, ngt and save it into the original sumstat file
# awk 'BEGIN {OFS="\t";} {tmp=$3; $3=$4; $4=$5; $5=tmp; print $0}' OFS='\t' ${sumstats_out}2 > ${sumstats_out}
# # 3)Remove unneccessary file
# rm ${sumstats_out}2

# # 3.2 MDD2 + 23andMe
# # ---------------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/MDD2_2018.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/MDD2/00_input/MDD2_2018_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/MDD2/00_input/pmb_info06_all_chr_merged_mdd2_filtered"
# snps="${out_genotypes}/3_LDpred/MDD2/00_input/pmb_info06_MDD2_common_SNPs.txt"

# # Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/MDD2/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # 3.3 BIP2018
# # --------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/BIP_2018.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/BIP-2018/00_input/BIP_2018_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/BIP-2018/00_input/pmb_info06_all_chr_merged_bip_filtered"
# snps="${out_genotypes}/3_LDpred/BIP-2018/00_input/pmb_info06_BIP_2018_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/BIP-2018/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # 3.4 BIPvsSCZ
# # ---------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/BIPvsSCZ/00_input/BDvsSCZ_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/BIPvsSCZ/00_input/pmb_info06_all_chr_merged_bd_vs_scz_filtered"
# snps="${out_genotypes}/3_LDpred/BIPvsSCZ/00_input/pmb_info06_BIPvsSCZ_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/BIPvsSCZ/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # 3.5 BIP+SCZvsCONT
# # --------------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/BDSCZvsCONT.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/BIP+SCZvsCONT/00_input/BDSCZvsCONT_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/BIP+SCZvsCONT/00_input/pmb_info06_all_chr_merged_bd_scz_vs_cont_filtered"
# snps="${out_genotypes}/3_LDpred/BIP+SCZvsCONT/00_input/pmb_info06_BDSCZvsCONT_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/BIP+SCZvsCONT/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # 3.6 ADHD2017
# # ---------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/ADHD_2017.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/ADHD-2017/00_input/ADHD_2017_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/ADHD-2017/00_input/pmb_info06_all_chr_merged_adhd_filtered"
# snps="${out_genotypes}/3_LDpred/ADHD-2017/00_input/pmb_info06_ADHD_2017_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/ADHD-2017/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # For ADHD we need to transform the summary statistic file into one of the formats coord_genotypes.py (LDpred) can read in (here: PGC)
# # 1) We add chr to the chromosome number and save it into a new file
# awk 'NR>1 {$1="chr"$1} 1' OFS='\t' $sumstats_out > ${sumstats_out}2
# # 2) Add new column at the end of the file, which includes only zeros and save it into a new file
# awk 'BEGIN{FS=OFS="0"}{print $0 "\t" OFS value}' ${sumstats_out}2 > ${sumstats_out}3
# # 3) Change order of columns to: chr, snp, a1, a2, bp, info, or, se, p, ngt and save it into the original sumstat file
# awk -F " " '{print $1, $2, $4, $5, $3, $6, $7, $8, $9, $10}' OFS="\t" ${sumstats_out}3 > ${sumstats_out}
# # 4) Change header names directly in the file
# sed -i "1s/.*/hg19chr\tsnpid\ta1\ta2\tbp\tor\tse\tp\tngt/" ${sumstats_out}
# # 5)Remove unneccessary files (with ending 3)
# rm ${sumstats_out}2
# rm ${sumstats_out}3

# # 3.7 ASD
# # ----------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/ASD.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/ASD/00_input/ASD_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/ASD/00_input/pmb_info06_all_chr_merged_asd_filtered"
# snps="${out_genotypes}/3_LDpred/ASD/00_input/pmb_info06_ASD_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/ASD/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # 3.8 BMI
# # ----------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/BMI.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/BMI/00_input/BMI_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/BMI/00_input/pmb_info06_all_chr_merged_bmi_filtered"
# snps="${out_genotypes}/3_LDpred/BMI/00_input/pmb_info06_BMI_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/BMI/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # Transform BMI sumstats into GIANT input format for LDpred
# # 1. Remove column 1 and 3 (chr and pos)
# awk '{$1=$3=""; print $0}' OFS="\t" $sumstats_out > ${sumstats_out}2
# # 2. Add column including nothing
# awk 'BEGIN{FS=OFS="NA"}{print $0 "\t" OFS value}' ${sumstats_out}2 > ${sumstats_out}3
# # 3. Reorder columns to: snpid, a1, a2, freq, effect, se, p, n
# awk 'BEGIN {OFS="\t";} {tmp=$4; $4=$8; $8=$7; $7=$6; $6=$5; $5=tmp; print $0}' ${sumstats_out}3 > ${sumstats_out}
# # 4) Change header names directly in the file
# sed -i "1s/.*/MarkerName\tAllele1\tAllele2\tb\teffect\tSE\tp\tN/" ${sumstats_out}
# # 5. Remove unneccessary files
# rm ${sumstats_out}2
# rm ${sumstats_out}3

# # 3.9 Height
# # -------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/Height.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/Height/00_input/Height_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/Height/00_input/pmb_info06_all_chr_merged_height_filtered"
# snps="${out_genotypes}/3_LDpred/Height/00_input/pmb_info06_Height_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/Height/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # Here, we need to remove the added chromosome and position columns again, save it into a summary statistics file number 2 and use it for the LDpred run
# cut -f 2,4,5,6,7,8,9,10 $sumstats_out > ${sumstats_out}2
# mv ${sumstats_out}2 $sumstats_out

# # 3.10 Educational Attainment (EA or IQ)
# # ---------------------------------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/EA.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/EA/00_input/EA_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/EA/00_input/pmb_info06_all_chr_merged_ea_filtered"
# snps="${out_genotypes}/3_LDpred/EA/00_input/pmb_info06_EA_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/BMI/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}

# # Transform EA sumstats into GIANT input format for LDpred
# # 1. Remove column 1 and 3 (chr and pos)
# awk '{$1=$3=""; print $0}' OFS="\t" $sumstats_out > ${sumstats_out}2
# # 2. Add column including number of individuals (n = 766,345)
# awk 'NR>1 {FS=OFS="766345"}{print $0 "\t" OFS value}' ${sumstats_out}2 > ${sumstats_out}3
# # 3. Reorder columns to: snpid, a1, a2, freq, effect, se, p, n
# awk 'BEGIN {OFS="\t";} {tmp=$4; $4=$8; $8=$7; $7=$6; $6=$5; $5=tmp; print $0}' OFS='\t' ${sumstats_out}3 > ${sumstats_out}
# # 4) Change header names directly in the file
# sed -i "1s/.*/MarkerName\tAllele1\tAllele2\tFreq.Allele1\tb\tSE\tp\tN/" ${sumstats_out}
# # 5. Remove unneccessary files
# rm ${sumstats_out}2
# rm ${sumstats_out}3

# # 3.11 Cross-Disorder
# # --------------------
# sumstats_in="${out_genotypes}/2.2_sumstats_changed_IDs/Cross_disorder_hg19.sumstats"
# sumstats_out="${out_genotypes}/3_LDpred/Cross_Disorder/00_input/Cross_disorder_filtered.sumstats"
# geno_out="${out_genotypes}/3_LDpred/Cross_Disorder/00_input/pmb_info06_all_chr_merged_cross_disorder_filtered"
# snps="${out_genotypes}/3_LDpred/Cross_Disorder/00_input/pmb_info06_cross_disorder_common_SNPs.txt"

# # 1) Run R-Script with the given parameters
# Rscript ${script_path}/01_overlaps_prepare_genotypes.R ${sumstats_in} ${sumstats_out} ${geno_in}.bim ${snps} &>> ${out_genotypes}/3_LDpred/Cross_Disorder/00_input/01_overlaps_prepare_genotypes.log

# # 2) Run PLINK to get overlapping SNPs for the genotype file + store data in .bed/.bim/.fam format
# $tools_path/PLINK/plink --bfile ${geno_in} --extract ${snps} --make-bed --out ${geno_out}


