# Use packages from own library (R version 3.6)
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# First, install all necessary packages as well as lassosum (only done once in the beginning)
#install.packages(c("RcppArmadillo", "data.table", "Matrix"), dependencies = TRUE)
# Then, download lassosum in the terminal using following command:
# wget https://github.com/tshmak/lassosum/releases/download/v0.4.4/lassosum_0.4.4.tar.gz
#install.packages("/home/icb/karolina.worf/R/lassosum_0.4.4.tar.gz", repos = NULL, type = "source")

# Load packages
library(data.table)
library(Matrix)
library(lassosum)
library(rcompanion)
library(parallel)
library(optparse)

option_list = list(
  make_option(c("-b", "--base"), type = "character", default = NULL, 
              help = "path and name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-t", "--target"), type = "character", default = NULL, 
              help = "path and name of target genotype file", metavar = "character"),
  make_option(c("-l", "--ld"), type = "character", default = NULL, 
              help = "path and name of LD reference genotype file", metavar = "character"),
  make_option(c("-y", "--pheno"), type = "character", default = NULL, 
              help = "path and name of phenotype file from target samples", metavar = "character"),
  make_option(c("-i", "--cov"), type = "character", default = NULL, 
              help = "path and name of covariate file for target samples", metavar = "character"),
  make_option(c("-d", "--thread"), type = "numeric", default = 2, 
              help = "number of threads to use for calculation (default = %default)", metavar = "character"),
  make_option(c("-r", "--A1"), type = "character", default = NULL, 
              help = "reference allele (A1) column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-a", "--A2"), type = "character", default = NULL, 
              help = "alternate allele (A2) column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-x", "--chr"), type = "character", default = NULL, 
              help = "chromosome column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-n", "--size"), type = "numeric", default = NULL, 
              help = "number of samples used for the GWAS summary statistics (size)", metavar = "character"),
  make_option(c("-e", "--bp"), type = "character", default = NULL, 
              help = "base pair (bp) column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-q", "--stat"), type = "character", default = NULL, 
              help = "statistics (beta or odds ratio = OR) column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-p", "--pvalue"), type = "character", default = NULL, 
              help = "p-value column name of GWAS summary statistics file", metavar = "character"),
  make_option(c("-g", "--gwas"), type = "character", default = NULL, 
              help = "name of GWAS", metavar = "character"),
  make_option(c("-o", "--out"), type = "character", default = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/05_lassosum", 
              help = "path to store the output files into (default: %default)", metavar = "character"),
  make_option(c("-f", "--beta"), type = "character", default = "no", 
              help = "variable to check if statistic is beta or odds ratio, if statistic is beta \"yes\" if not \"no\" (default: \"no\")", metavar = "character")
); 

opt_parser = OptionParser(option_list = option_list, add_help_option = FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$base)){
  print_help(opt_parser)
  stop("At least one argument for GWAS summary statistics file must be supplied (--base).", call.=FALSE)
}
if (is.null(opt$target)){
  print_help(opt_parser)
  stop("At least one argument for target genotype file must be supplied (--target).", call.=FALSE)
}
if (is.null(opt$ld)){
  print_help(opt_parser)
  stop("At least one argument for LD reference genotype file must be supplied (--ld).", call.=FALSE)
}
if (is.null(opt$pheno)){
  print_help(opt_parser)
  stop("At least one argument for phenotype file must be supplied (--pheno).", call.=FALSE)
}
if (is.null(opt$cov)){
  print_help(opt_parser)
  stop("At least one argument for covariates file must be supplied (--cov).", call.=FALSE)
}
if (is.null(opt$A1)){
  print_help(opt_parser)
  stop("At least one argument for the reference allele column name of GWAS sumstats file, e.g. A1 or ref, must be supplied (--A1).", call.=FALSE)
}
if (is.null(opt$A2)){
  print_help(opt_parser)
  stop("At least one argument for the alternate allele column name of GWAS sumstats file, e.g. A2 or alt, must be supplied (--A2).", call.=FALSE)
}
if (is.null(opt$chr)){
  print_help(opt_parser)
  stop("At least one argument for chromosome column name of GWAS sumstats file must be supplied (--chr).", call.=FALSE)
}
if (is.null(opt$size)){
  print_help(opt_parser)
  stop("At least one argument for the number of samples used in GWAS sumstats file must be supplied (--size).", call.=FALSE)
}
if (is.null(opt$bp)){
  print_help(opt_parser)
  stop("At least one argument for the position column name of GWAS sumstats file, e.g. bp, must be supplied (--bp).", call.=FALSE)
}
if (is.null(opt$stat)){
  print_help(opt_parser)
  stop("At least one argument for the statistics column name of GWAS sumstats file, e.g. BETA or OR, must be supplied (--stat).", call.=FALSE)
}
if (is.null(opt$pvalue)){
  print_help(opt_parser)
  stop("At least one argument for the p-value column name of GWAS summary statistics file must be supplied (--pvalue).", call.=FALSE)
}
if (is.null(opt$gwas)){
  print_help(opt_parser)
  stop("At least one argument for the name of the GWAS must be supplied (--gwas).", call.=FALSE)
}

# get all GWAS column names
pvalue <- opt$pvalue
chr <- opt$chr
ref_allele <- opt$A1
alt_allele <- opt$A2
bp <- opt$bp
stat <- opt$stat

# load GWAS summary statistics file
sumstats <- fread(opt$base, data.table = FALSE)
# remove rows with p-value = 0 as it causes problems in the transformation
sumstats <- sumstats[!sumstats[, pvalue] == 0, ]

### Run lassosum (parallelized) ###
# Get LD region file (Other alternatives available: type ?lassosum.pipeline for more details)
LDblocks <- "EUR.hg19" # This will use LD regions as defined in Berisa and Pickrell (2015) for the European population and the hg19 genome
# check if statistic is beta or odds ratio (or)
if (opt$beta == "yes") {
  # Calculate SNP-wise correlations from p-values (n = sample size) without specifying sign
  cor <- p2cor(p = sumstats[, pvalue], n = opt$size)
  # Run lassosum (parallelized with number of threads, default = 2)
  cl <- makeCluster(opt$thread, type = "FORK")
  out <- lassosum.pipeline(cor = cor, chr = sumstats[, chr], pos = sumstats[, bp], 
                           A1 = sumstats[, ref_allele], A2 = sumstats[, alt_allele], # A2 is not required but advised
                           ref.bfile = opt$ld, test.bfile = opt$target, 
                           LDblocks = LDblocks, cluster = cl)
} else {
  # Calculate SNP-wise correlations from p-values (n = sample size) using odds ratio (OR)
  cor <- p2cor(p = sumstats[, pvalue], n = opt$size, sign = log(sumstats[, stat]))
  # Run lassosum (parallelized with number of threads, default = 2)
  cl <- makeCluster(opt$thread, type = "FORK")
  out <- lassosum.pipeline(cor = cor, chr = sumstats[, chr], pos = sumstats[, bp], 
                           A1 = sumstats[, ref_allele], A2 = sumstats[, alt_allele], # A2 is not required but advised
                           ref.bfile = opt$ld, test.bfile = opt$target, 
                           LDblocks = LDblocks, cluster = cl)
}
save(out, file = paste0(opt$out, "/lassosum_pipeline_result_", opt$gwas, ".RData"))

### Generate PRS ###
# generate PRS with covariates 
val <- validate(out, pheno = opt$pheno, cov = opt$cov)
# save output as .RData file
#save(val, file = paste0(opt$out, "/lassosum_PRS_with_covars_", opt$gwas, ".RData"))
# save PRS per sample in a .best file
write.table(val$results.table, file = paste0(opt$out, "/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_", opt$gwas, "_lassosum_PRS_with_cov.best"), sep = "\t", quote = FALSE)
# save plot comparing all four s
pdf(file = paste0(opt$out, "/plot_lassosum_PRS_with_covars_SCZ2-2014.pdf"))
plot(val)
dev.off()
# generate PRS without covariates
val2 <- validate(out, pheno = opt$pheno)
# save output as .RData file
#save(val2, file = paste0(opt$out, "/lassosum_PRS_without_covars_", opt$gwas, ".RData"))
# save PRS per sample in a .best file
write.table(val2$results.table, file = paste0(opt$out, "/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_", opt$gwas, "_lassosum_PRS_without_cov.best"), sep = "\t", quote = FALSE)
# save plot comparing all four s values
pdf(file = paste0(opt$out, "/plot_lassosum_PRS_without_covars_", opt$gwas, ".pdf"))
plot(val2)
dev.off()
