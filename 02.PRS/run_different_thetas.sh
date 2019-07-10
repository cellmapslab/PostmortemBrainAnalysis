# 1. Run LDpred: change causal fraction from 0.1 to 0.9 in 01 steps (0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9)
python3 /storage/groups/ccm01/tools/ldpred-master/LDpred.py gibbs --cf /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/01_coord/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered.coord_file --ldr 1000 --f 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 --ldf /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered --N 150063 --out /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered &> /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/02_LDpred/LDpred.log
echo "LDpred.py done"

# 2. Run Validation (get score for each patient)
python3 /storage/groups/ccm01/tools/ldpred-master/LDpred.py score --gf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/00_input/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_scz2_filtered --rf /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/02_LDpred/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered --f 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 --pf /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/phenotypes.txt --out /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_PRS &> /storage/scratch/users/karolina.worf/MPIP/02.PRS/LDpred/SCZ2/03_validate/validate.log
echo "validate.py done"

