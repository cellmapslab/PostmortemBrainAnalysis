#!/bin/bash

#SBATCH -o PRS-figures.o
#SBATCH -e PRS-figures.e
#SBATCH -J PRS-figures
#SBATCH -p cpu_p
#SBATCH --mem=10000 # total memory in MB
#SBATCH -t 03:00:00 # format dd-hh:mm:ss
#SBATCH --nice=10000  # adjusts scheduling priority

# load R version 3.6.1 and its packages
source /lustre/groups/ccm01/code/PostmortemBrainAnalysis/02.PRS/04a_makeenv.sh

Rscript -e "rmarkdown::render('07_PRS_Analysis.Rmd')"

