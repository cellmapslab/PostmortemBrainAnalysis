################################
### Generate GRanges objects ###
################################
# This R-Script generates GRanges objects for the eQTL data as well as several public datasets (CMC, GTEx, VEP)
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args frq=<path/maf-file> geno=<path/to/geno-files> eqtls=<path/to/eQTL-files> clump=<path/to/clumped-eQTL-files> cmc=<path/to/CMC-files> gtex=<path/to/gtex-v8> cortex=<path/to/cortex-files> gwas=<path/to/gwas-sumstats> outpath=<path/for/output>" 05_generate_GRanges_objects.R &>> generate_GRanges_objects.log
# or using Rscript: Rscript --slave 04a_generate_result_table.R frq=<path/maf-file> geno=<path/to/geno-files> eqtls=<path/to/eQTL-files> clump=<path/to/clumped-eQTL-files> cmc=<path/to/CMC-files> gtex=<path/to/gtex-v8> cortex=<path/to/cortex-files> gwas=<path/to/gwas-sumstats> outpath=<path/for/output> &>> generate_GRanges_objects.log

# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")
