#!/bin/bash
# Define neccessary paths
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"
tools_path="/storage/groups/ccm01/tools"
out_path="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred"


### ADHD-2017 ###
# 20,183 cases and 35,191 controls: N = 55,374 samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/ADHD-2017/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_adhd_filtered -g ${out_path}/ADHD-2017/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_adhd_filtered -s ${out_path}/ADHD-2017/00_input/ADHD_2017_filtered.sumstats -n 55374 -f PGC_OLD -o ${out_path}/ADHD-2017 -t $tools_path

### ASD-2019 ###
# 18,382 cases and 27,969 controls: N = 46,351
bash ${script_path}/05_LDpred.sh -c ${out_path}/ASD/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_asd_filtered -g ${out_path}/ASD/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_asd_filtered -s ${out_path}/ASD/00_input/ASD_filtered.sumstats -n 46351 -f PGC -o ${out_path}/ASD -t $tools_path
#wait 

### BMI-2015 ###
# 87,048 samples (no case/control here!)
bash ${script_path}/05_LDpred.sh -c ${out_path}/BMI/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bmi_filtered -g ${out_path}/BMI/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bmi_filtered -s ${out_path}/BMI/00_input/BMI_filtered.sumstats -n 87048 -f GIANT -o ${out_path}/BMI -t $tools_path

# ### BIP 2018 ###
# 20,352 cases and 31,358 controls: N = 51,710 samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/BIP-2018/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_filtered -g ${out_path}/BIP-2018/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_filtered -s ${out_path}/BIP-2018/00_input/BIP_2018_filtered.sumstats -n 51710 -f PGC -o ${out_path}/BIP-2018 -t $tools_path
#wait

### Cross-Disorder ###
# 61,220 samples (33,332 cases and 27,888 controls) => PGC_OLD
bash ${script_path}/05_LDpred.sh -c ${out_path}/Cross_Disorder/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_cd_filtered -g ${out_path}/Cross_Disorder/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_cd_filtered -s ${out_path}/Cross_Disorder/00_input/CD_filtered.sumstats -n 61220 -f PGC_OLD -o ${out_path}/Cross_Disorder -t $tools_path

### BIP+SCZvsCONT ###
# 20,129 BIP cases and 33,426 SCZ cases and 54,065 controls: N = 107,620 samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/BIP+SCZvsCONT/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_scz_vs_cont_filtered -g ${out_path}/BIP+SCZvsCONT/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_scz_vs_cont_filtered -s ${out_path}/BIP+SCZvsCONT/00_input/BIPSCZvsCONT_filtered.sumstats -n 107620 -f PGC -o ${out_path}/BIP+SCZvsCONT -t $tools_path

### BIPvsSCZ ###
# 15,270 BIP cases and 23,585 SCZ cases: N = 38,855 samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/BIPvsSCZ/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_vs_scz_filtered -g ${out_path}/BIPvsSCZ/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_bip_vs_scz_filtered -s ${out_path}/BIPvsSCZ/00_input/BIPvsSCZ_filtered.sumstats -n 38855 -f PGC -o ${out_path}/BIPvsSCZ -t $tools_path
#wait 

### Educational Attainment (EA or IQ) ###
# 766,345 samples (no case/control, without 23andMe samples)
bash ${script_path}/05_LDpred.sh -c ${out_path}/EA/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_ea_filtered -g ${out_path}/EA/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_ea_filtered -s ${out_path}/EA/00_input/EA_filtered.sumstats -n 766345 -f GIANT -o ${out_path}/EA -t $tools_path
#wait 

### Height ###
# 253,288 samples (no case/control here!)
bash ${script_path}/LDpred.sh -g ${out_path}/Height/00_input/pmb_info06_all_chr_merged_height_filtered -s ${out_path}/Height/00_input/Height_filtered.sumstats -n 253288 -f GIANT -o ${out_path}/Height -t $tools_path 

### MDD1-2012 ###
# xxx cases and xxx controls: N = xxx samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/MDD2/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_mdd2_filtered -g ${out_path}/MDD2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_mdd2_filtered -s ${out_path}/MDD2/00_input/MDD2_2018_filtered.sumstats -n 480359 -f PGC -o ${out_path}/MDD2 -t $tools_path

### MDD2-2018 ###
# xxx cases and xxx controls: N = xxx samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/MDD2/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_mdd2_filtered -g ${out_path}/MDD2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_mdd2_filtered -s ${out_path}/MDD2/00_input/MDD2_2018_filtered.sumstats -n 480359 -f PGC -o ${out_path}/MDD2 -t $tools_path

### MDD3-2019 ###
# xxx cases and xxx controls: N = xxx samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/MDD2/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_mdd2_filtered -g ${out_path}/MDD2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_mdd2_filtered -s ${out_path}/MDD2/00_input/MDD2_2018_filtered.sumstats -n 480359 -f PGC -o ${out_path}/MDD2 -t $tools_path

### SCZ2-2014 ###
# xxx cases and xxx controls: N = xxx samples
bash ${script_path}/05_LDpred.sh -c ${out_path}/SCZ2/00_input/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered -g ${out_path}/SCZ2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_scz2_filtered -s ${out_path}/SCZ2/00_input/SCZ2_filtered.sumstats -n 150063 -f PGC_OLD -o ${out_path}/SCZ2 -t $tools_path

