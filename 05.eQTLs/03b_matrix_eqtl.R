################################
### Script to run MatrixEQTL ###
################################
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(MatrixEQTL)
library(data.table)

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args=(commandArgs(trailingOnly=TRUE))	

# go through the list of arguments
if(length(args)==0){
    cat("\nScript to run MatrixEQTL\n\nYou can invoke the script as follows:\nRscript 03b_matrix_eqtl.R snps=\"/home/.../SNP.file\" spos=\"/home/.../SNP-positions.file\" exp=\"/home/.../expression.file\" expos=\"/home/.../expression-position.file\" cov=\"/home/.../covariates.file\" ciso=\"/home/.../cis_eqtl_output.file\" transo=\"/home/.../trans_eqtl_output.file\"\n\nCommandline parameters:\n\nsnps\t\tpath and file name of the SNP genotype data\nspos\t\tpath and file name of the positions for the SNP genotype data\nexp\t\tpath and file name of the expression data\nexpos\t\tpath and file name of the positions for the expression data\ncov\t\tpath and file name of the covariates\nciso\t\tpath and output file name for the cis-eQTLs\ntranso\t\tpath and output file name for the trans-eQTLs\n\n")
} else {
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
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
            noFDRsaveMemory = FALSE
        )
        return(me)
    }

    # Load SNP and gene positions
    # MatrixEQTL does not work without giving the tool the SNP and gene position. 
    snpspos <- read.table(spos, header = T, stringsAsFactors = FALSE, sep = "\t") # snp chr pos
    head(snpspos)
    # load exon location file
    genepos <- read.table(expos, header = T, stringsAsFactors = FALSE, sep = "\t") # geneid chr start end
    genepos <- genepos[!is.na(genepos$start), ]  # remove NA's
    head(genepos)

    # Run matrixEQTL For cis-eQTLs
    me.all <- run_matrix_eqtl(SNP_file = snps, expression_file = exp, cov_file = cov, outfile_cis = ciso, outfile_trans = transo, cis = 1, trans = 0)
}
