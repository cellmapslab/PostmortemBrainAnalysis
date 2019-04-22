#!/bin/bash
# Define neccessary paths
script_path="/storage/groups/ccm01/code/MPIP-Dean_postmortem_brain/PRS"
tools_path="/storage/groups/ccm01/tools"
out_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/genotype_data/3_LDpred"

### SCZ2 ###
# 36,988 cases and 113,075 controls: N = 150,063 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/SCZ2/00_input/pmb_info06_all_chr_merged_scz2_filtered_dup_removed -s ${out_path}/SCZ2/00_input/SCZ2_filtered.sumstats -n 150063 -f PGC -o ${out_path}/SCZ2 -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/SCZ2/00_input/pmb_info06_all_chr_merged_scz2_filtered -s ${out_path}/SCZ2/00_input/SCZ2_filtered.sumstats -n 150063 -f PGC -o ${out_path}/SCZ2 -t $tools_path &

### MDD2 ###
# 135,458 cases and 344,901 controls: N = 480,359 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/MDD2/00_input/pmb_info06_all_chr_merged_mdd2_filtered_dup_removed -s ${out_path}/MDD2/00_input/MDD2_2018_filtered.sumstats -n 480359 -f PGC_large -o ${out_path}/MDD2 -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/MDD2/00_input/pmb_info06_all_chr_merged_mdd2_filtered -s ${out_path}/MDD2/00_input/MDD2_2018_filtered.sumstats -n 480359 -f PGC_large -o ${out_path}/MDD2 -t $tools_path 

# ### BIP 2018 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/BIP-2018/00_input/pmb_info06_all_chr_merged_bip_filtered_dup_removed -s ${out_path}/BIP-2018/00_input/BIP_2018_filtered.sumstats -n 51710 -f PGC_large -o ${out_path}/BIP-2018 -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/BIP-2018/00_input/pmb_info06_all_chr_merged_bip_filtered -s ${out_path}/BIP-2018/00_input/BIP_2018_filtered.sumstats -n 51710 -f PGC_large -o ${out_path}/BIP-2018 -t $tools_path &

### BIP+SCZvsCONT ###
# 20,129 BIP cases and 33,426 SCZ cases and 54,065 controls: N = 107,620 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/BIP+SCZvsCONT/00_input/pmb_info06_all_chr_merged_bd_scz_vs_cont_filtered_dup_removed -s ${out_path}/BIP+SCZvsCONT/00_input/BDSCZvsCONT_filtered.sumstats -n 107620 -f PGC_large -o ${out_path}/BIP+SCZvsCONT -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/BIP+SCZvsCONT/00_input/pmb_info06_all_chr_merged_bd_scz_vs_cont_filtered -s ${out_path}/BIP+SCZvsCONT/00_input/BDSCZvsCONT_filtered.sumstats -n 107620 -f PGC_large -o ${out_path}/BIP+SCZvsCONT -t $tools_path 

### BIPvsSCZ ###
# 15,270 BIP cases and 23,585 SCZ cases: N = 38,855 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/BIPvsSCZ/00_input/pmb_info06_all_chr_merged_bd_vs_scz_filtered_dup_removed -s ${out_path}/BIPvsSCZ/00_input/BDvsSCZ_filtered.sumstats -n 38855 -f PGC_large -o ${out_path}/BIPvsSCZ -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/BIPvsSCZ/00_input/pmb_info06_all_chr_merged_bd_vs_scz_filtered -s ${out_path}/BIPvsSCZ/00_input/BDvsSCZ_filtered.sumstats -n 38855 -f PGC_large -o ${out_path}/BIPvsSCZ -t $tools_path &

### ADHD 2017 ###
# 20,183 cases and 35,191 controls: N = 55,374 samples
# Run LDpred with the transformed summary statistics file number 2
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/ADHD-2017/00_input/pmb_info06_all_chr_merged_adhd_filtered_dup_removed -s ${out_path}/ADHD-2017/00_input/ADHD_2017_filtered2.sumstats -n 55374 -f PGC -o ${out_path}/ADHD-2017 -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/ADHD-2017/00_input/pmb_info06_all_chr_merged_adhd_filtered -s ${out_path}/ADHD-2017/00_input/ADHD_2017_filtered.sumstats -n 55374 -f PGC -o ${out_path}/ADHD-2017 -t $tools_path 

### ASD ###
# 5,305 cases and 5305 pseudocontrols: N = 10,610 samples
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/ASD/00_input/pmb_info06_all_chr_merged_asd_filtered_dup_removed -s ${out_path}/ASD/00_input/ASD_filtered.sumstats -n 10610 -f PGC_large -o ${out_path}/ASD -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/ASD/00_input/pmb_info06_all_chr_merged_asd_filtered -s ${out_path}/ASD/00_input/ASD_filtered.sumstats -n 10610 -f PGC_large -o ${out_path}/ASD -t $tools_path &

### BMI ###
# 87,048 samples (no case/control here!)
# qsub -e "." -o "." -hard -l job_mem=50G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/BMI/00_input/pmb_info06_all_chr_merged_bmi_filtered_dup_removed -s ${out_path}/BMI/00_input/BMI_filtered.sumstats -n 87048 -f GIANT -o ${out_path}/BMI -t $tools_path &  # run all in parallel
bash ${script_path}/LDpred.sh -g ${out_path}/BMI/00_input/pmb_info06_all_chr_merged_bmi_filtered -s ${out_path}/BMI/00_input/BMI_filtered.sumstats -n 87048 -f GIANT -o ${out_path}/BMI -t $tools_path

### Height ###
# 253,288 samples (no case/control here!)
# qsub -e "." -o "." -hard -l job_mem=32G -q long_fed25 ${script_path}/LDpred.sh -g ${out_path}/Height/00_input/pmb_info06_all_chr_merged_height_filtered_dup_removed -s ${out_path}/Height/00_input/Height_filtered.sumstats -n 253288 -f GIANT -o ${out_path}/Height -t $tools_path &  # run all in parallel
# bash ${script_path}/LDpred.sh -g ${out_path}/Height/00_input/pmb_info06_all_chr_merged_height_filtered -s ${out_path}/Height/00_input/Height_filtered.sumstats -n 253288 -f GIANT -o ${out_path}/Height -t $tools_path 

### Educational Attainment (EA or IQ) ###
# 766,345 samples (no case/control, without 23andMe samples)
bash ${script_path}/LDpred.sh -g ${out_path}/EA/00_input/pmb_info06_all_chr_merged_ea_filtered -s ${out_path}/EA/00_input/EA_filtered.sumstats -n 766345 -f GIANT -o ${out_path}/EA -t $tools_path &

### Cross-Disorder ###
# 61,220 samples (33,332 cases and 27,888 controls)
bash ${script_path}/LDpred.sh -g ${out_path}/Cross_Disorder/00_input/pmb_info06_all_chr_merged_cross_disorder_filtered -s ${out_path}/Cross_Disorder/00_input/Cross_disorder_filtered.sumstats -n 61220 -f PGC -o ${out_path}/Cross_Disorder -t $tools_path


