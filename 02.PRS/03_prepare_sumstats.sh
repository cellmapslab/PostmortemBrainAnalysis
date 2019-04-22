#!/bin/bash

out_genotypes="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"
tools_path="/storage/groups/ccm01/tools"
sum_stats_path="/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/"

chromosome=($(seq 1 1 22))


#######################################
### Prepare Summary Statistic Files ### 
#######################################
# 1. Reorder and format summary statistic files
# Summary Statistic Formats for LDpred: 
# "Standard", "BASIC", "PGC", "PGC_large", "GIANT", "GIANT2" und "DECODE"
#######################################

# 2.1.1 SCZ2 
# ------------
# 1) SCZ2 summary statistic file, change columns to order: chr - SNP ID - pos (atm: SNP ID - chr - pos) and remove chr before chromosome number
awk 'BEGIN {OFS="\t";} {tmp=$5; $5=$4; $4=$3; $3=tmp; print $0}' OFS='\t' ${sum_stats_path}/SCZ2/scz2.snp.results.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results.txt
# 2) Remove chr before chromosome number
awk 'NR>1 {gsub("chr", "", $1)} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results2.txt
# 3) Delete old SCZ file
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results.txt
# 4) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the SCZ2 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results2.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/SCZ2.sumstats
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/scz2.snp.results2.txt  # remove old SCZ2 file
# 5) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/SCZ2.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/SCZ2.sumstats

# 2.1.2 MDD2 + 23andMe
# ---------------------
# 1) Remove last four columns
awk 'NF{NF-=4};1' OFS="\t" ${sum_stats_path}/MDD2/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe
# 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the MDD2 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/MDD2_2018.sumstats
# 3) Remove old MDD file
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe
# 4) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/MDD2_2018.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/MDD2_2018.sumstats

# 2.1.3 BIP2018
# --------------
# 1) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the BIP2018 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${sum_stats_path}/BIP-2018/daner_PGC_BIP32b_mds7a_0416a > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BIP_2018.sumstats
# 2) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BIP_2018.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/BIP_2018.sumstats

# 2.1.4 BIPvsSCZ
# ---------------
# 1) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the BIPvsSCZ summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${sum_stats_path}/BIP+SCZ/BDvsSCZ.sumstats > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BDvsSCZ.sumstats
# 2) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BDvsSCZ.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats
# Remove last two columns (some do not have a value there and plink does not like it...)
awk 'NR>1 {$12=$13=""} 1' OFS="\t" ${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats2
mv ${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats2 ${out_genotypes}/2.2_sumstats_changed_IDs/BDvsSCZ.sumstats

# 2.1.5 BIP+SCZvsCONT
# --------------------
# 1) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the SCZ2 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${sum_stats_path}/BIP+SCZ/BDSCZvsCONT.sumstats > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BDSCZvsCONT.sumstats
# 2) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BDSCZvsCONT.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/BDSCZvsCONT.sumstats

# 2.1.6 ADHD2017
# ---------------
# Keep in mind: ADHD summary statistic file has sometimes randomly "m"'s behind the SNP IDs!!! 
# As we change the SNP IDs completely it is no problem, but here is an example code for removing the m's: 
# awk '{gsub("m", "", $2); print}' OFS="\t" input-file.txt
# 1) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the ADHD2017 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${sum_stats_path}/ADHD-2017/adhd_jul2017 > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ADHD_2017.sumstats
# 2) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ADHD_2017.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/ADHD_2017.sumstats

# 2.1.7 ASD
# ----------
# 1) ASD summary statistic file has no zeros on floating point numbers (e.g. .123 instead of 0.123), so we need to add 0's before the dot in the columns 6-11
awk '{ if ($6 ~ /^\./) $6="0"$6; if ($7 ~ /^\./) $7="0"$7; if ($8 ~ /^\./) $8="0"$8; if ($9 ~ /^\./) $9="0"$9; if ($10 ~ /^\./) $10="0"$10; if ($11 ~ /^\./) $11="0"$11; print $0}' OFS="\t" ${sum_stats_path}/ASD/daner_pgc_asd_euro_all_25Mar2015 > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_asd_euro_all_25Mar2015
# 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the SCZ2 summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_asd_euro_all_25Mar2015 > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ASD.sumstats
# 3) Remove old ASD file
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/daner_pgc_asd_euro_all_25Mar2015
# 4) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ASD.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/ASD.sumstats

# 2.1.8 BMI
# ----------
# 1) BMI summary statistic file, change columns to order: chr - SNP ID - pos (atm: SNP ID - chr - pos)
awk 'BEGIN {OFS="\t";} {tmp=$1; $1=$2; $2=tmp; print $0;}' OFS="\t" ${sum_stats_path}/BMI/ENGAGE1000G_BMI.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ENGAGE1000G_BMI.txt
# 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the BMI summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ENGAGE1000G_BMI.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BMI.sumstats
# 3) Remove old BMI file
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/ENGAGE1000G_BMI.txt 
# 4) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/BMI.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/BMI.sumstats

# # 2.1.9 Height
# # -------------
# # 1) As Height GWAS summary statistic has hg18 (NCBI 36) SNP locations instead of hg19 (GRCh37), we need to lift them over => run Rscript
# Rscript "/storage/groups/ccm01/code/MPIP-Dean_postmortem_brain/PRS/convert_hg18_to_hg19.R"
# # 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the Height summary statistic file
# awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats2
# # 3) Change height summary statistics file name back to origin
# mv ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats2 ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats  
# # 4) Change all SNP IDs in the summary statistic files to following example format: chr10_10001753_C_T & leave first line untouched
# awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/Height.sumstats

# 2.1.10 Educational Attainment (EA or also known as IQ)
# -------------------------------------------------------
# 1) EA summary statistics file, change columns to order: CHR - MarkerName - POS (atm: MarkerName - CHR - POS)
awk 'BEGIN {OFS="\t";} {tmp=$1; $1=$2; $2=tmp; print $0;}' OFS="\t" ${sum_stats_path}/EA/GWAS_EA_excl23andMe.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/GWAS_EA_excl23andMe.txt
# 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the EA summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/GWAS_EA_excl23andMe.txt > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/EA.sumstats
# 3) Remove old BMI file
rm ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/GWAS_EA_excl23andMe.txt
# 4) Change all SNP IDs in the summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/EA.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/EA.sumstats

echo "Transform summary statistic files done"

# 2.1.11 Cross-Disorder
# ----------------------
# 1) As the Cross-Disorder GWAS summary statistic has hg18 (NCBI 36) SNP locations instead of hg19 (GRCh37), we need to lift them over => run Rscript
Rscript "/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS/03_convert_hg18_to_hg19.R"
# 2) Remove all deletions (D), insertions (I) and other strange alleles (<CN0>, <CN2>, <INS..:>) in the Height summary statistic file
awk '{ if (!index($2, "D") && !index($2, "I") && !index($2, "<")) { print $0; } }' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats > ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats2
# 3) Change cross-disorder summary statistics file name back to origin
mv ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats2 ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats  
# 4) Change all SNP IDs in the summary statistics file to following example format: chr10:10001753 & leave first line untouched
awk 'NR>1 {$2=$1":"$5} 1' OFS="\t" ${out_genotypes}/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats > ${out_genotypes}/2.2_sumstats_changed_IDs/Cross_disorder_hg19.sumstats


