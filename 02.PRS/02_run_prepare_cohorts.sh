#!/bin/bash
cohort="/storage/groups/ccm01/datasets/2018_cohort_GSK_MDD"
outpath="/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS"
# Generate the output directory (if not existing)
outpath_cohort="${outpath}/02_preprocessed_cohort"  # cohort output path
mkdir -p $outpath_cohort
tools_path="/storage/groups/ccm01/tools"
script_path="/storage/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS"

# Preprocess the cohort genotypes: filter for MAF, HWE, merge chromosome files and remove duplicates
bash ${script_path}/02a_prepare_cohorts.sh -c ${cohort} -o ${outpath_cohort} -t ${tools_path} -s ${script_path}

