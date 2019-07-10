.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

library(MatrixEQTL)
library(data.table)

# Define function
run_matrix_eqtl <- function(SNP_file, expression_file, cov_file, outfile_cis, outfile_trans, cis, trans){
  # 1. Set general parameters
  # Use linear model (other models: modelANOVA, modelLINEAR, or modelLINEAR_CROSS)
  useModel <- modelLINEAR
  # We need to set up the p-value cutoff and the error model (in this case assuming independent errors)
  # Only associations significant at this level will be saved
  pvOutputThreshold_cis = cis  # cis-eQTLS cutoff: 0.05 = 5e-2
  pvOutputThreshold_tra = trans  # trans-eQTLs cutoff: 0.01 = 1e-2 (0 means no trans-eQTLs)
  # Error covariance matrix: set to numeric() for identity
  errorCovariance <- numeric()
  # Distance for local (cis) gene-SNP pairs
  cisDist = 1e6  # cis window of 1Mb = 1e6, 100kb = 1e5
  
  # 2. Set Data Up
  # Now we need to set up the snp and gene expression data in the special format required by the MatrixEQTL package
  # 2.1 SNP Data
  snps = SlicedData$new()
  snps$fileDelimiter <- "\t"     # the TAB character
  snps$fileOmitCharacters <- "NA" # denote missing values;
  snps$fileSkipRows <- 1          # one row of column labels
  snps$fileSkipColumns <- 1       # one column of row labels
  snps$fileSliceSize <- 100000     # read file in pieces of 100,000 rows
  snps$LoadFile(SNP_file)
  # 2.2 Gene Expression Data
  gene <- SlicedData$new()
  gene$fileDelimiter <- "\t"      # the TAB character
  gene$fileOmitCharacters <- "NA" # denote missing values;
  gene$fileSkipRows <- 1          # one row of column labels
  gene$fileSkipColumns <- 1      # one column of row labels
  gene$fileSliceSize <- 2000      # read file in pieces of 2,000 rows
  gene$LoadFile(expression_file)
  # 2.3 Covariates
  cvrt <- SlicedData$new()
  cvrt$fileDelimiter = "\t"      # the TAB character
  cvrt$fileOmitCharacters = "NA" # denote missing values;
  cvrt$fileSkipRows = 1          # one row of column labels
  cvrt$fileSkipColumns = 1       # one column of row labels
  if(length(cov_file) > 0){
    cvrt$LoadFile(cov_file)  # read file if given
  }
  
  # 3. Run Matrix_eQTL
  # Gene Expr x SNP
  me = Matrix_eQTL_main(
    snps = snps, 
    gene = gene, 
    cvrt = cvrt,
    output_file_name = outfile_trans,
    pvOutputThreshold = pvOutputThreshold_tra , #pvOutputThreshold_tra only cis
    useModel = useModel, 
    errorCovariance = errorCovariance, 
    verbose = TRUE, 
    output_file_name.cis = outfile_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos, 
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  return(me)
}

# Load SNP and gene positions
# MatrixEQTL does not work without giving the tool the SNP and gene position. 
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/SNP_locations.txt", header = T, stringsAsFactors = FALSE, sep = "\t") # snp chr pos
head(snpspos)
genepos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/gene_locations.txt", header = T, stringsAsFactors = FALSE, sep = "\t") # geneid chr start end
genepos <- genepos[!is.na(genepos$start), ]  # remove NA's
head(genepos)

# 1. Calculate eQTLS For Complete Data (SCZ + AFF + CTL)
# Genotype (SNP) file name
SNP_file_name_all <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_snps_SCZ_CTL.txt"
# Gene expression file name
expression_file_name_all <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_gene_expr_SCZ_CTL.txt"
# Covariates file name
covariates_file_name_all <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_cov_SCZ_CTL.txt"
# Output file name
output_file_name_cis_all <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/01_pmb_SCZ_CTL_matrixEQTL_cis_results_genes.txt"
output_file_name_trans_all <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/01_pmb_SCZ_CTL_matrixEQTL_trans_results_genes.txt"
# Run matrixEQTL only for cis and save results into an RData file
me.all <- run_matrix_eqtl(SNP_file = SNP_file_name_all, expression_file = expression_file_name_all, cov_file = covariates_file_name_all, outfile_cis = output_file_name_cis_all, outfile_trans = output_file_name_trans_all, cis = 1, trans = 5e-5)
#save(me.all, file = "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/01_pmb_SCZ_CTL_matrixEQTL.RData")

# 2. Only SCZ Data
# Genotype (SNP) file name
SNP_file_name_scz <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_snps_SCZ.txt"
# Gene expression file name
expression_file_name_scz <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_gene_expr_SCZ.txt"
# Covariates file name
covariates_file_name_scz <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_cov_SCZ.txt"
# Output file name
output_file_name_cis_scz <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/02_pmb_SCZ_matrixEQTL_cis_results_genes.txt"
output_file_name_trans_scz <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/02_pmb_SCZ_matrixEQTL_trans_results_genes.txt"
# Run matrixEQTL only for cis-eQTLs and save results into an RData file
me.scz <- run_matrix_eqtl(SNP_file = SNP_file_name_scz, expression_file = expression_file_name_scz, cov_file = covariates_file_name_scz, outfile_cis = output_file_name_cis_scz, outfile_trans = output_file_name_trans_scz, cis = 1, trans = 5e-5)
#save(me.scz, file = "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/02_pmb_SCZ_matrixEQTL.RData")

# 3. Only CTL Data
# Genotype (SNP) file name
SNP_file_name_ctl <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_snps_CTL.txt"
# Gene expression file name
expression_file_name_ctl <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_gene_expr_CTL.txt"
# Covariates file name
covariates_file_name_ctl <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_cov_CTL.txt"
# Output file name
output_file_name_cis_ctl <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/03_pmb_CTL_matrixEQTL_cis_results_genes.txt"
output_file_name_trans_ctl <- "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/03_pmb_CTL_matrixEQTL_trans_results_genes.txt"
# Run matrixEQTL only for cis and save results into an RData file
me.ctl <- run_matrix_eqtl(SNP_file = SNP_file_name_ctl, expression_file = expression_file_name_ctl, cov_file = covariates_file_name_ctl, outfile_cis = output_file_name_cis_ctl, outfile_trans = output_file_name_trans_ctl, cis = 1, trans = 5e-5)
#save(me.ctl, file = "/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/results/03_pmb_CTL_matrixEQTL.RData")

