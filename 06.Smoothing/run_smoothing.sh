#!/bin/bash

# load R version 3.5.1 and its packages
source /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/01_makeenv.sh

# run smoothing for all cis-eqtls + filtered trans-eQTLs
R CMD BATCH /storage/groups/ccm01/code/PostmortemBrainAnalysis/06.Smoothing/02.smoothing.R

