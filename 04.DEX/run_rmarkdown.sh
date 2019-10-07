#!/bin/bash

# run Rmarkdown
/home/icb/karolina.worf/R-3.6.1/bin/exec/R --slave --no-save --no-restore -e "rmakrdown::render('/storage/groups/ccm01/code/PostmortemBrainAnalysis/04.DEX/01_Preprocessing.Rmd', encoding = 'UTF-8')"

