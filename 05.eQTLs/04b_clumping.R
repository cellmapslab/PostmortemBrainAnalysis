####################################################################
### Generate eQTL result tables from the MatrixEQTL output files ###
####################################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args input=<path/sample-file> input2=<path/sample-sheet> output1=<path/output-file1> output2=<path/output-file2> 04b_clumping.R
# or using Rscript: Rscript --slave 04b_clumping.R input1=\"/home/.../input-file1.txt\" input2=\"/home/.../input2.txt\" output1=\"/home/.../out-file1.txt\" output2=\"/home/.../out-file2.txt\"

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly=TRUE))	

# go through the list of arguments
if(length(args) == 0){
    cat("\nScript to perform clumping of the SNPs from the eQTL result tables (also called eSNPs).\nThe script performs clumping and generates result tables for the clumped eSNPs.\n\nYou can invoke the script as follows:\nRscript 04a_generate_result_tables.R inpath=\"/home/.../to/input-files\" snps=\"/home/.../SNP.file\" gtf=\"/home/.../gencode.gtf\" gpos=\"/home/.../gene-positions.file\" tpos=\"/home/.../transcript-positions.file\" epos=\"/home/.../exons-positions.file\" outpath=\"/home/.../for/output-files\" \n\nCommandline parameters:\n\ninpath\t\tpath, where the input files are stored\nsnps\t\tpath and file name of the SNP genotype data\ngtf\t\tpath and file name of the gencode .gtf file\ngpos\t\tpath and file name of the positions for all genes\ntpos\t\tpath and file name of the positions for all transcripts\nepos\t\tpath and file name of the positions for all exons\noutpath\t\tpath to store the output files\n\n")
} else {
    for(i in 1:length(args)){
        eval(parse(text = args[[i]]))
    }
    ### Functions ###
    # function to get all unique SNPs of the eQTL results
    GetUniqueSnps <- function(eqtls, file){
        unique_snps <- as.data.frame(unique(eqtls$`SNP ID`))
        colnames(unique_snps) <- "SNP"
        fwrite(unique_snps, file = paste0(outpath1, file), col.names = F, sep = "\n")
        return(unique_snps)
    }

    # function to get the chromosome and minimal p-value for all unique SNPs of the eQTL results
    GetChrAndPvalueOfUniqueSnps <- function(eqtls, snps, file){
        min_pval_snps <- eqtls[ , .SD[which.min(`P-value`)], by = `SNP ID`]  # get minimal p-value for each SNP
        unique_snps <- merge(snps, min_pval_snps[, c("SNP ID", "SNP Location", "P-value")], by.x = "SNP", by.y = "SNP ID")
        unique_snps$CHR <- sapply(unique_snps$`SNP Location`, function(x) strsplit(x, ":")[[1]][1])
        unique_snps$P <- unique_snps$`P-value`
        unique_snps <- unique_snps[, c("SNP", "CHR", "P")]
        print(head(unique_snps))
        print(dim(unique_snps))
        fwrite(unique_snps, file = paste0(outpath1, file), sep = "\t")
    }

    # function to load the clumping results, reduce eQTL results to top LD SNPs and save the resulting eQTLs into an extra file
    GetClumpedSnpsEqtls <- function(input, eqtls, output){
        clumped_snps <- fread(paste0(outpath2, input), header = T, select = "SNP", stringsAsFactors = F)
        eqtls <- eqtls[eqtls$`SNP ID` %in% clumped_snps$SNP, ]
        print(dim(eqtls))
        write.table(eqtls, file = paste0(outpath3, output), row.names = F, quote = F, dec = ".", sep = "\t")
    }

    ### 1. LOAD eQTL results ###
    # 1.1 gene, transcript and exon cis-eQTLs of case vs. control data (ALL)
    eqtls_gene <- fread(paste0(inpath, "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    eqtls_transcript <- fread(paste0(inpath, "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    eqtls_exon <- fread(paste0(inpath, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    # 1.2 gene, transcript and exon cis-eQTLs of case data only
    eqtls_gene_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_gene_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    eqtls_transcript_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_transcript_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    eqtls_exon_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_exon_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    # 1.3 gene, transcript and exon cis-eQTLs of control data only
    eqtls_gene_controls <- fread(paste0(inpath, "/03_pmb_cis-eQTL_gene_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    eqtls_transcript_controls <- fread(paste0(inpath, "/03_pmb_cis-eQTL_transcript_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    eqtls_exon_controls <- fread(paste0(inpath, "/03_pmb_cis-eQTL_exon_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    print("1. Done loading eQTL results")
    
    # 2. Generate list of all SNP rs IDs, which we want to do clumping for
    # 2.1 gene, transcript and exon cis-eQTLs of case vs. control data
    gene_snps <- GetUniqueSnps(eqtls_gene, file = "/01_pmb_cis-eQTL_genes_SNP_IDs.txt")
    transcript_snps <- GetUniqueSnps(eqtls_transcript, file = "/01_pmb_cis-eQTL_transcripts_SNP_IDs.txt")
    exon_snps <- GetUniqueSnps(eqtls_exon, file = "/01_pmb_cis-eQTL_exons_SNP_IDs.txt")
    # 2.2 gene, transcript and exon cis-eQTLs of case data only
    gene_cases_snps <- GetUniqueSnps(eqtls_gene_cases, file = "/02_pmb_cis-eQTL_genes_cases_only_SNP_IDs.txt")
    transcript_cases_snps <- GetUniqueSnps(eqtls_transcript_cases, file = "/02_pmb_cis-eQTL_transcripts_cases_only_SNP_IDs.txt")
    exon_cases_snps <- GetUniqueSnps(eqtls_exon_cases, file = "/02_pmb_cis-eQTL_exons_cases_only_SNP_IDs.txt")
    # 2.3 gene, transcript and exon cis-eQTLs of control data only
    gene_controls_snps <- GetUniqueSnps(eqtls_gene_controls, file = "/03_pmb_cis-eQTL_genes_ctls_only_SNP_IDs.txt")
    transcript_controls_snps <- GetUniqueSnps(eqtls_transcript_controls, file = "/03_pmb_cis-eQTL_transcripts_ctls_only_SNP_IDs.txt")
    exon_controls_snps <- GetUniqueSnps(eqtls_exon_controls, file = "/03_pmb_cis-eQTL_exons_ctls_only_SNP_IDs.txt")
    print("2. Done generating eSNP files for clumping")
    
    # 3. Generate list of all SNP rs IDs, which we want to do clumping for including their corrresponding chromosome and minimal eQTL p-value
    # 3.1 gene, transcript and exon cis-eQTLs of case vs. control data
    GetChrAndPvalueOfUniqueSnps(eqtls_gene, gene_snps, file = "/01_pmb_cis-eQTL_genes_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_transcript, transcript_snps, file = "/01_pmb_cis-eQTL_transcripts_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_exon, exon_snps, file = "/01_pmb_cis-eQTL_exons_clumping.txt")
    # 3.2 gene, transcript and exon cis-eQTLs of case data only
    GetChrAndPvalueOfUniqueSnps(eqtls_gene_cases, gene_cases_snps, file = "/02_pmb_cis-eQTL_genes_cases_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_transcript_cases, transcript_cases_snps, file = "/02_pmb_cis-eQTL_transcripts_cases_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_exon_cases, exon_cases_snps, file = "/02_pmb_cis-eQTL_exons_cases_clumping.txt")
    # 3.3 gene, transcript and exon cis-eQTLs of control data only
    GetChrAndPvalueOfUniqueSnps(eqtls_gene_controls, gene_controls_snps, file = "/03_pmb_cis-eQTL_genes_ctls_only_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_transcript_controls, transcript_controls_snps, file = "/03_pmb_cis-eQTL_transcripts_ctls_only_clumping.txt")
    GetChrAndPvalueOfUniqueSnps(eqtls_exon_controls, exon_controls_snps, file = "/03_pmb_cis-eQTL_exons_ctls_only_clumping.txt")
    print("3. Done generating eSNP files including chromosome and minimal p-value")    
    
    # remove unnecessary data
    rm(gene_snps, transcript_snps, exon_snps, gene_cases_snps, transcript_cases_snps, exon_cases_snps, gene_controls_snps, transcript_controls_snps, exon_controls_snps)

    # 3. Run clumping using PLINK
    # 3.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/01_pmb_cis-eQTL_genes_SNP_IDs.txt --clump ", outpath1, "/01_pmb_cis-eQTL_genes_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/01_pmb_cis-eQTL_genes"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/01_pmb_cis-eQTL_transcripts_SNP_IDs.txt --clump ", outpath1, "/01_pmb_cis-eQTL_transcripts_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/01_pmb_cis-eQTL_transcripts"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/01_pmb_cis-eQTL_exons_SNP_IDs.txt --clump ", outpath1, "/01_pmb_cis-eQTL_exons_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/01_pmb_cis-eQTL_exons"))
    # 3.2 cis-eQTLs of cases only in genes, transcripts and exons
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/02_pmb_cis-eQTL_genes_cases_only_SNP_IDs.txt --clump ", outpath1, "/02_pmb_cis-eQTL_genes_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/02_pmb_cis-eQTL_genes_cases_only"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/02_pmb_cis-eQTL_transcripts_cases_only_SNP_IDs.txt --clump ", outpath1, "/02_pmb_cis-eQTL_transcripts_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/02_pmb_cis-eQTL_transcripts_cases_only"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/02_pmb_cis-eQTL_exons_cases_only_SNP_IDs.txt --clump ", outpath1, "/02_pmb_cis-eQTL_exons_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/02_pmb_cis-eQTL_exons_cases_only"))
    # 3.3 gene, transcript and exon cis-eQTLs of control data only
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/03_pmb_cis-eQTL_genes_ctls_only_SNP_IDs.txt --clump ", outpath1, "/03_pmb_cis-eQTL_genes_ctls_only_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/03_pmb_cis-eQTL_genes_controls_only"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/03_pmb_cis-eQTL_transcripts_ctls_only_SNP_IDs.txt --clump ", outpath1, "/03_pmb_cis-eQTL_transcripts_ctls_only_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/03_pmb_cis-eQTL_transcripts_controls_only"))
    system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile ", geno, " --extract ", outpath1, "/03_pmb_cis-eQTL_exons_ctls_only_SNP_IDs.txt --clump ", outpath1, "/03_pmb_cis-eQTL_exons_ctls_only_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out ", outpath2, "/03_pmb_cis-eQTL_exons_controls_only"))
    print("3. Done with PLINK clumping")    
    
    # 4. Get clumped SNP eQTLs
    # 4.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
    GetClumpedSnpsEqtls(input = "/01_pmb_cis-eQTL_genes.clumped", eqtls_gene, output = "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls_clumped.txt")
    GetClumpedSnpsEqtls(input = "/01_pmb_cis-eQTL_transcripts.clumped", eqtls_transcript, output = "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls_clumped.txt")
    GetClumpedSnpsEqtls(input = "/01_pmb_cis-eQTL_exons.clumped", eqtls_exon, output = "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls_clumped.txt")
    # 4.2 cis-eQTLs of cases only in genes, transcripts and exons
    GetClumpedSnpsEqtls(input = "/02_pmb_cis-eQTL_genes_cases_only.clumped", eqtls_gene_cases, output = "/02_pmb_cis-eQTL_gene_hits_cases_only_clumped.txt")
    GetClumpedSnpsEqtls(input = "/02_pmb_cis-eQTL_transcripts_cases_only.clumped", eqtls_transcript_cases, output = "/02_pmb_cis-eQTL_transcript_hits_cases_only_clumped.txt")
    GetClumpedSnpsEqtls(input = "/02_pmb_cis-eQTL_exons_cases_only.clumped", eqtls_exon_cases, output = "/02_pmb_cis-eQTL_exon_hits_cases_only_clumped.txt")
    # 4.3 
    GetClumpedSnpsEqtls(input = "/03_pmb_cis-eQTL_genes_controls_only.clumped", eqtls_gene_controls, output = "/03_pmb_cis-eQTL_gene_hits_controls_only_clumped.txt")
    GetClumpedSnpsEqtls(input = "/03_pmb_cis-eQTL_transcripts_controls_only.clumped", eqtls_transcript_controls, output = "/03_pmb_cis-eQTL_transcript_hits_controls_only_clumped.txt")
    GetClumpedSnpsEqtls(input = "/03_pmb_cis-eQTL_exons_controls_only.clumped", eqtls_exon_controls, output = "/03_pmb_cis-eQTL_exon_hits_controls_only_clumped.txt")
    print("4. Done generating result tables")
}
