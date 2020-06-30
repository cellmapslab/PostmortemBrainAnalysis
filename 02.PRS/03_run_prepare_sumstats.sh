#!/bin/bash
sum_stats_path="/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS"
# Generate the output directories (if not existing already)
outpath_sumstats="${outpath}/03_preprocessed_sumstats"  # preprocessed sumstats output path
mkdir -p $outpath_sumstats
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"
tools_path="/storage/groups/ccm01/tools"

#######################################
### Prepare Summary Statistic Files ### 
#######################################
# Summary statistic formats  for LDpred v1.0.11: (own definition: CUSTOM)
# LDPRED
# CHR POS SNP_ID REF ALT REF_FRQ PVAL BETA SE N
# chr1 1020428 rs6687776 C T 0.85083 0.0587 -0.0100048507289348 0.0100 8000
# 
# STANDARD
# chr pos ref alt reffrq info rs pval effalt
# chr1 1020428 C T 0.85083 0.98732 rs6687776 0.0587 -0.0100048507289348
# 
# BASIC
# hg19chrc snpid a1 a2 bp or p
# chr1 rs4951859 C G 729679 0.97853 0.2083
# 
# PGC
# CHR SNP BP A1 A2 FRQ_A_30232 FRQ_U_40578 INFO OR SE P ngt Direction HetISqt HetChiSq HetDf HetPVa
# ...
# 
# GIANT
# MarkerName Allele1 Allele2 Freq.Allele1.HapMapCEU p N
# rs10 a c 0.0333 0.8826 78380 
#######################################

###############
### AD-2019 ###
###############
sumstats="${sum_stats_path}/AD-2019/AD_sumstats_Jansenetal_2019sept.txt"
output_sumstats="${outpath_sumstats}/AD-2019.sumstats"
# Change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$6="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### ADHD-2019 ###
#################
sumstats="${sum_stats_path}/ADHD-2019/adhd_jul2017"
output_sumstats="${outpath_sumstats}/ADHD-2019.sumstats"
# Change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###################
### ALCDEP-2018 ###
###################
sumstats="${sum_stats_path}/ALCDEP-2018/pgc_alcdep.trans_fe_unrel_geno.aug2018_release.txt"
output_sumstats="${outpath_sumstats}/ALCDEP-2018.sumstats"
# Change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###############
### AN-2019 ###
###############
sumstats="${sum_stats_path}/AN-2019/pgcAN2.2019-07.vcf.tsv"
output_sumstats="${outpath_sumstats}/AN-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
sed -n '/^ *[^#]/p' ${sumstats} | awk 'NR==1 {$1=$1}; NR>1 {$3="chr"$1":"$2} 1' OFS="\t" > ${output_sumstats}

################
### ASD-2019 ###
################
sumstats="${sum_stats_path}/ASD-2019/iPSYCH-PGC_ASD_Nov2017"
output_sumstats="${outpath_sumstats}/ASD-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

##################
### AUDIT-2019 ###
##################
sumstats="${sum_stats_path}/AUDIT-2019/AUDIT_UKB_2018_AJP_with_positions.txt"
output_sumstats="${outpath_sumstats}/AUDIT-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {split($2,a,"_"); $2="chr"a[1]} 1' OFS="\t" ${sumstats} > ${output_sumstats}

################
### BIP-2019 ###
################
sumstats="${sum_stats_path}/BIP-2019/daner_PGC_BIP32b_mds7a_0416a"
output_sumstats="${outpath_sumstats}/BIP-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

####################
### BIP+SCZ-2018 ###
####################
# 1. BIPvsSCZ
# ------------
sumstats="${sum_stats_path}/BIP+SCZ-2018/BDvsSCZ.sumstats"
output_sumstats="${outpath_sumstats}/BIPvsSCZ-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

# 2. BIP+SCZvsCONT
# -----------------
sumstats="${sum_stats_path}/BIP+SCZ-2018/BDSCZvsCONT.sumstats"
output_sumstats="${outpath_sumstats}/BIPSCZvsCONT-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

########################
### BMI-2015 (ENGAGE) ###
########################
sumstats="${sum_stats_path}/BMI-2015_ENGAGE/ENGAGE1000G_BMI.txt"
output_sumstats="${outpath_sumstats}/BMI-2015_ENGAGE.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t"  ${sumstats} > ${output_sumstats}

########################
### BMI-2015 (GIANT) ###
########################
sumstats="${sum_stats_path}/BMI-2015_GIANT/SNP_gwas_mc_merge_nogc_with_locations.tbl.uniq"
output_sumstats="${outpath_sumstats}/BMI-2015_GIANT.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t"  ${sumstats} > ${output_sumstats}

################
### BMI-2018 ###
################
sumstats="${sum_stats_path}/BMI-2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"
output_sumstats="${outpath_sumstats}/BMI-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$3="chr"$1":"$2} 1' OFS="\t"  ${sumstats} > ${output_sumstats}

################
### BMI-2019 ###
################
sumstats="${sum_stats_path}/BMI-2019/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt"
output_sumstats="${outpath_sumstats}/BMI-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$3="chr"$1":"$2} 1' OFS="\t"  ${sumstats} > ${output_sumstats}

################
### CAD-2017 ###
################
sumstats="${sum_stats_path}/CAD-2017/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.hg19.OR.txt"
output_sumstats="${outpath_sumstats}/CAD-2017.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$3":"$4} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#####################
### Cannabis-2018 ###
#####################
sumstats="${sum_stats_path}/Cannabis-2018/Cannabis_ICC_23andmetop_UKB_het.txt"
output_sumstats="${outpath_sumstats}/Cannabis-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### CDG1-2013 ###
#################
sumstats="${sum_stats_path}/CDG1-2013/pgc.cross.full.2013-03.hg19.txt"
output_sumstats="${outpath_sumstats}/CDG1-2013.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### CDG2-2019 ###
#################
sumstats="${sum_stats_path}/CDG2-2019/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt"
output_sumstats="${outpath_sumstats}/CDG2-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$3="chr"$1":"$2} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###############
### EA-2018 ###
###############
sumstats="${sum_stats_path}/EA-2018/GWAS_EA_excl23andMe.txt"
output_sumstats="${outpath_sumstats}/EA-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###############
### ED-2017 ###
###############
sumstats="${sum_stats_path}/ED-2017/pgc.ed.freeze1.summarystatistics.July2017.txt"
output_sumstats="${outpath_sumstats}/ED-2017.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###################
### Height-2014 ###
###################
sumstats="${sum_stats_path}/Height-2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_hg19_with_locations.txt"
output_sumstats="${outpath_sumstats}/Height-2014.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#########################
### Intelligence-2018 ###
#########################
sumstats="${sum_stats_path}/Intelligence-2018/SavageJansen_2018_intelligence_metaanalysis.txt"
output_sumstats="${outpath_sumstats}/Intelligence-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$3":"$4} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### MDD1-2012 ###
#################
sumstats="${sum_stats_path}/MDD1-2012/pgc.mdd.full.2012-04.hg19.txt"
output_sumstats="${outpath_sumstats}/MDD1-2012.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### MDD2-2018 ###
#################
sumstats="${sum_stats_path}/MDD2-2018/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe"
output_sumstats="${outpath_sumstats}/MDD2-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### MDD3-2019 ###
#################
sumstats="${sum_stats_path}/MDD3-2019/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB"
output_sumstats="${outpath_sumstats}/MDD3-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

################
### OCD-2017 ###
################
sumstats="${sum_stats_path}/OCD-2017/ocd_aug2017"
output_sumstats="${outpath_sumstats}/OCD-2017.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### PTSD-2019 ###
#################
sumstats="${sum_stats_path}/PTSD-2019/pts_all_freeze2_overall.results"
output_sumstats="${outpath_sumstats}/PTSD-2019.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2="chr"$1":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

#################
### SCZ2-2014 ### 
#################
sumstats="${sum_stats_path}/SCZ2-2014/scz2.snp.results.txt"
output_sumstats="${outpath_sumstats}/SCZ2-2014.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$2=$1":"$5} 1' OFS="\t" ${sumstats} > ${output_sumstats}

################
### T2D-2017 ###
################
sumstats="${sum_stats_path}/T2D-2017/METAANALYSIS_DIAGRAM_SE1_BP_POS.txt"
output_sumstats="${outpath_sumstats}/T2D-2017.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$1} 1' OFS="\t" ${sumstats} > ${output_sumstats}

###############
### TS-2018 ###
###############
sumstats="${sum_stats_path}/TS-2018/TS_Oct2018"
output_sumstats="${outpath_sumstats}/TS-2018.sumstats"
# Remove comment lines and change all SNP IDs in the GWAS summary statistic file to following example format: chr:pos, e.g. chr10:10001753 & leave first line untouched
awk 'NR==1 {$1=$1}; NR>1 {$1="chr"$2":"$3} 1' OFS="\t" ${sumstats} > ${output_sumstats}

