################################
### Generate GRanges objects ###
################################
# This R-Script generates GRanges objects for the eQTL data as well as several public datasets (CMC, GTEx, VEP)
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args frq=<path/maf-file> geno=<path/to/geno-files> eqtls=<path/to/eQTL-files> clump=<path/to/clumped-eQTL-files> cmc=<path/to/CMC-files> gtex=<path/to/gtex-v8> cortex=<path/to/cortex-files> gwas=<path/to/gwas-sumstats> outpath=<path/for/output>" 05_generate_GRanges_objects.R &>> generate_GRanges_objects.log
# or using Rscript: Rscript --slave 04a_generate_result_table.R frq=<path/maf-file> geno=<path/to/geno-files> eqtls=<path/to/eQTL-files> clump=<path/to/clumped-eQTL-files> cmc=<path/to/CMC-files> gtex=<path/to/gtex-v8> cortex=<path/to/cortex-files> gwas=<path/to/gwas-sumstats> outpath=<path/for/output> &>> generate_GRanges_objects.log

# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(IRanges)
library(arules)
library(reshape2)
library(biomaRt) # to get SNP positions

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
    cat("\nScript to generate GRanges objects from the eQTL data as well as public datasets, such as CMC, GTEx, VEP, etc.\n\nYou can invoke the script as follows:\nRscript 05_generate_GRanges_objects.R inpath=\"/home/.../to/input-files\" frq=\"/home/.../to/maf.frq/\" geno=\"path/to/genotype-files\" eqtls=\"path/to/eqtl-files\" clump=\"path/to/clumping-files\" cmc=\"path/to/cmc-files\" gtex=\"path/to/gtes-files\" cortex=\"/home/.../cortex.file\" vep=\"path/to/vep.file\" chromhmm=\"path/to/chromhmm-files\" gwas=\"path/to/gwas-files\" outpath=\"path/to/store/output\" &>> generate_GRanges_objects.log \n\nCommandline parameters:\n\nfrq\t\tfrequency file including the maf values for all SNPs\ngeno\t\tpath and file name without ending of the SNP genotype files in PLINK format (.bim/.bed/.fam)\neqtls\t\tpath to eqtl result files\nclump\t\tpath to clumped eQTL result files\ncmc\t\tpath to CMC (Fromer et al.) eQTL and isoQTL files\ngtex\t\tpath to GTEx v8 eQTL and sQTL files\ncortex\t\tpath to Cortex meQTL files\ngwas\t\tpath to GWAS summary statistic files\noutpath\t\tpath to store all output files into\n\n")
} else {
    for(i in 1:length(args)){
        eval(parse(text = args[[i]]))
    }
    # function to generate an eQTL GRanges object (with MAF bins)
    # @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
    # @param mdata Merging data (e.g. Background SNPs)
    # @param ofile Name for the rds output file of the generated GRanges object 
    GenerateGrangesObjecteQTL <- function(input, mdata, ofile){
        # merge with MAF info
        input <- merge(input, mdata, by = "SNP")
        # generate GRanges object
        input_gr <- GenomicRanges::GRanges(seqnames = input$CHR, 
                                           ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)), 
                                                                     end = as.numeric(as.character(input$POS))),
                                           snp_id = input$SNP, 
                                           bin = input$bin)
        # remove duplicates
        input_gr <- input_gr[!duplicated(input_gr)] 
        GenomeInfoDb::seqlevelsStyle(input_gr) <- "NCBI"
        # save GRanges object to an rds file
        saveRDS(input_gr, file = paste0(outpath, "/", ofile))
        # return finalized GRanges object
        return(input_gr)
    }

    # function to generate a GRanges object (for all datasets with Chromosome, base pair position and rs SNP ID information)
    # @param input Input data of which an GRanges object should be made of (needs CHR, BP, P and SNP as column names)
    # @param ofile Name for the rds output file of the generated GRanges object 
    GenerateGrangesObject <- function(input, ofile){
        # filter for p-value (P) < 0.05
        input <- input[input$P < 0.05, ]
        # remove all rows including NA in any column
        input <- input[complete.cases(input), ]
        # generate GRanges object
        input_gr <- GenomicRanges::GRanges(seqnames = input$CHR,
                                           ranges = IRanges::IRanges(start = as.numeric(as.character(input$BP)),
                                                                     end = as.numeric(as.character(input$BP))),
                                           snp_id = input$SNP,
                                           p_value = input$P)  # generate GRanges object only with chr, position, rs SNP ID and p-value
        # remove duplicates
        input_gr <- input_gr[!duplicated(input_gr)] 
        GenomeInfoDb::seqlevelsStyle(input_gr) <- "NCBI"
        print(length(input_gr)) # print length
        # save GRanges object to an rds file
        saveRDS(input_gr, file = paste0(outpath, "/", ofile))
        rm(input_gr)
    }

    ### Load BiomaRt SNP informations
    snpmart <- useMart(biomart = "ENSEMBL_MART_SNP", 
                       host = "grch37.ensembl.org", 
                       path = "/biomart/martservice", 
                       dataset = "hsapiens_snp")


    #######################
    ### Background SNPs ###
    #######################
    # Read MAF values for each SNP
    mafs <- data.table::fread(frq, data.table = F)

    # generate bins of MAF in 0.05 steps
    boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
    mafs$bin <- arules::discretize(mafs$MAF, method = "fixed", breaks = boundaries, labels = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

    # load background SNPs
    snps <- data.table::fread(paste0(geno, ".bim"), data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos

    # Generate background SNP data as GRanges object
    background_gr <- GenerateGrangesObjecteQTL(input = snps, 
                                               mdata = mafs, 
                                               ofile = "background_SNPs.rds")
    table(seqnames(background_gr))  # no chr X/Y/M (remove in dataset if included there)
    length(background_gr$snp_id) # 6,798,048
    rm(snps)

    # Generate merging file for eQTL data from background SNPs 
    background <- as.data.frame(background_gr)[, c("seqnames", "start", "snp_id", "bin")]
    colnames(background) <- c("CHR", "POS", "SNP", "bin")
    print("Done with background SNPs")


    #################
    ### PMB eQTLs ###
    #################
    ### Gene eSNPs (geSNPs) ###
    # Cases vs. controls
    geSNPs <- fread(paste0(eqtls, "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    geSNPs_gr <- GenerateGrangesObjecteQTL(input = geSNPs,
                                           mdata = background,
                                           ofile = "01_pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
    length(geSNPs_gr$snp_id) 
    rm(geSNPs, geSNPs_gr)
    print("Done with geSNPs cases vs. controls")
    
    # Cases only
    geSNPs_cases <- fread(paste0(eqtls, "/02_pmb_cis-eQTL_gene_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    geSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = geSNPs_cases,
                                                 mdata = background,
                                                 ofile = "02_pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds")
    length(geSNPs_cases_gr$snp_id) 
    rm(geSNPs_cases, geSNPs_cases_gr)
    print("Done with geSNPs cases only")
    
    # Controls only
    geSNPs_controls <- fread(paste0(eqtls, "/03_pmb_cis-eQTL_gene_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    geSNPs_controls_gr <- GenerateGrangesObjecteQTL(input = geSNPs_controls,
                                                    mdata = background,
                                                    ofile = "03_pmb_gene_hits_controls_only_cis_eQTL_SNPs.rds")
    length(geSNPs_controls_gr$snp_id) 
    rm(geSNPs_controls, geSNPs_controls_gr)
    print("Done with geSNPs controls only")
    
    # Clumped cases vs. controls
    geSNPs_clumped <- fread(paste0(clump, "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls_clumped.txt"), header = T, stringsAsFactors = F)
    geSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = geSNPs_clumped,
                                                   mdata = background,
                                                   ofile = "01_pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
    length(geSNPs_clumped_gr$snp_id) 
    rm(geSNPs_clumped, geSNPs_clumped_gr)
    print("Done with geSNPs cases vs. controls clumped")
    
    # Clumped cases only
    geSNPs_cases_clumped <- fread(paste0(clump, "/02_pmb_cis-eQTL_gene_hits_cases_only_clumped.txt"), header = T, stringsAsFactors = F)
    geSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = geSNPs_cases_clumped,
                                                         mdata = background,
                                                         ofile = "02_pmb_gene_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
    length(geSNPs_cases_clumped_gr$snp_id) 
    rm(geSNPs_cases_clumped, geSNPs_cases_clumped_gr)
    print("Done with geSNPs cases only clumed")
    
    # Clumped controls only
    geSNPs_controls_clumped <- fread(paste0(clump, "/03_pmb_cis-eQTL_gene_hits_controls_only_clumped.txt"), header = T, stringsAsFactors = F)
    geSNPs_controls_clumped_gr <- GenerateGrangesObjecteQTL(input = geSNPs_controls_clumped,
                                                            mdata = background,
                                                            ofile = "03_pmb_gene_hits_controls_only_cis_eQTL_SNPs_clumped.rds")
    length(geSNPs_controls_clumped_gr$snp_id) 
    rm(geSNPs_controls_clumped, geSNPs_controls_clumped_gr)
    print("Done with geSNPs controls only clumed")

    ### Transcript eSNPs (teSNPs) ###
    # Cases vs. controls
    teSNPs <- fread(paste0(eqtls, "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    teSNPs_gr <- GenerateGrangesObjecteQTL(input = teSNPs,
                                           mdata = background,
                                           ofile = "01_pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
    length(teSNPs_gr$snp_id) 
    rm(teSNPs, teSNPs_gr)
    print("Done with teSNPs cases vs. controls")
    
    # Cases only
    teSNPs_cases <- fread(paste0(eqtls, "/02_pmb_cis-eQTL_transcript_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    teSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = teSNPs_cases,
                                                 mdata = background,
                                                 ofile = "02_pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds")
    length(teSNPs_cases_gr$snp_id) 
    rm(teSNPs_cases, teSNPs_cases_gr)
    print("Done with teSNPs cases only")
    
    # Controls only
    teSNPs_controls <- fread(paste0(eqtls, "/03_pmb_cis-eQTL_transcript_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    teSNPs_controls_gr <- GenerateGrangesObjecteQTL(input = teSNPs_controls,
                                                 mdata = background,
                                                 ofile = "03_pmb_transcript_hits_controls_only_cis_eQTL_SNPs.rds")
    length(teSNPs_controls_gr$snp_id) 
    rm(teSNPs_controls, teSNPs_controls_gr)
    print("Done with teSNPs controls only")
    
    # Clumped cases vs. controls
    teSNPs_clumped <- fread(paste0(clump, "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls_clumped.txt"), header = T, stringsAsFactors = F)
    teSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = teSNPs_clumped,
                                                   mdata = background,
                                                   ofile = "01_pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
    length(teSNPs_clumped_gr$snp_id) 
    rm(teSNPs_clumped, teSNPs_clumped_gr)
    print("Done with teSNPs cases vs. controls clumped")
    
    # Clumped cases only
    teSNPs_cases_clumped <- fread(paste0(clump, "/02_pmb_cis-eQTL_transcript_hits_cases_only_clumped.txt"), header = T, stringsAsFactors = F)
    teSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = teSNPs_cases_clumped,
                                                         mdata = background,
                                                         ofile = "02_pmb_transcript_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
    length(teSNPs_cases_clumped_gr$snp_id) 
    rm(teSNPs_cases_clumped, teSNPs_cases_clumped_gr)
    print("Done with teSNPs cases only clumped")
    
    # Clumped controls only
    teSNPs_controls_clumped <- fread(paste0(clump, "/03_pmb_cis-eQTL_transcript_hits_controls_only_clumped.txt"), header = T, stringsAsFactors = F)
    teSNPs_controls_clumped_gr <- GenerateGrangesObjecteQTL(input = teSNPs_controls_clumped,
                                                         mdata = background,
                                                         ofile = "03_pmb_transcript_hits_controls_only_cis_eQTL_SNPs_clumped.rds")
    length(teSNPs_controls_clumped_gr$snp_id) 
    rm(teSNPs_controls_clumped, teSNPs_controls_clumped_gr)
    print("Done with teSNPs controls only clumped")
    
    ### Exon eSNPs (eeSNPs) ###
    # Cases vs. controls
    eeSNPs <- fread(paste0(eqtls, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
    eeSNPs_gr <- GenerateGrangesObjecteQTL(input = eeSNPs,
                                       mdata = background,
                                       ofile = "01_pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
    length(eeSNPs_gr$snp_id) 
    rm(eeSNPs, eeSNPs_gr)
    print("Done with eeSNPs cases vs. controls")
    
    # Cases only
    eeSNPs_cases <- fread(paste0(eqtls, "/02_pmb_cis-eQTL_exon_hits_cases_only.txt"), header = T, stringsAsFactors = F)
    eeSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_cases,
                                                 mdata = background,
                                                 ofile = "02_pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds")
    length(eeSNPs_cases_gr$snp_id) 
    rm(eeSNPs_cases, eeSNPs_cases_gr)
    print("Done with eeSNPs cases only")
    
    # Controls only
    eeSNPs_controls <- fread(paste0(eqtls, "/03_pmb_cis-eQTL_exon_hits_controls_only.txt"), header = T, stringsAsFactors = F)
    eeSNPs_controls_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_controls,
                                                 mdata = background,
                                                 ofile = "03_pmb_exon_hits_controls_only_cis_eQTL_SNPs.rds")
    length(eeSNPs_controls_gr$snp_id) 
    rm(eeSNPs_controls, eeSNPs_controls_gr)
    print("Done with eeSNPs controls only")
    
    # Clumped cases vs. controls
    eeSNPs_clumped <- fread(paste0(clump, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls_clumped.txt"), header = T, stringsAsFactors = F)
    eeSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_clumped,
                                                   mdata = background,
                                                   ofile = "01_pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
    length(eeSNPs_clumped_gr$snp_id) 
    rm(eeSNPs_clumped, eeSNPs_clumped_gr)
    print("Done with eeSNPs cases vs. controls clumped")
    
    # Clumped cases only
    eeSNPs_cases_clumped <- fread(paste0(clump, "/02_pmb_cis-eQTL_exon_hits_cases_only_clumped.txt"), header = T, stringsAsFactors = F)
    eeSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_cases_clumped,
                                                         mdata = background,
                                                         ofile = "02_pmb_exon_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
    length(eeSNPs_cases_clumped_gr$snp_id) 
    rm(eeSNPs_cases_clumped, eeSNPs_cases_clumped_gr)
    print("Done with eeSNPs cases only clumped")
    
    # Clumped controls only
    eeSNPs_controls_clumped <- fread(paste0(clump, "/03_pmb_cis-eQTL_exon_hits_controls_only_clumped.txt"), header = T, stringsAsFactors = F)
    eeSNPs_controls_clumped_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_controls_clumped,
                                                            mdata = background,
                                                            ofile = "03_pmb_exon_hits_controls_only_cis_eQTL_SNPs_clumped.rds")
    length(eeSNPs_controls_clumped_gr$snp_id) 
    rm(eeSNPs_controls_clumped, eeSNPs_controls_clumped_gr)
    print("Done with eeSNPs controls only clumped")

    
    ###########
    ### CMC ###
    ###########
    # CMC (Fromer et al.) gene eQTLs: geSNPs
    cmc_gene <- data.table::fread(paste0(cmc, "/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt"))
    cmc_gene <- cmc_gene[cmc_gene$FDR < 0.05, ]  # filter for FDR <= 0.05
    cmc_gene_gr <- GenomicRanges::GRanges(seqnames = cmc_gene$SNP_location, 
                                      snp_id = cmc_gene$SNP)  # generate GRanges object 
    cmc_gene_gr <- cmc_gene_gr[!duplicated(cmc_gene_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
    GenomeInfoDb::seqlevelsStyle(cmc_gene_gr) <- "NCBI"
    seqlevels(cmc_gene_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cmc_gene_gr), c("X", "Y", "M"))  # remove SNPs located on chromosomes X, Y and M
    saveRDS(cmc_gene_gr, paste0(outpath, "/CMC_adjustedSVA_eQTL_SNPs.rds"))
    rm(cmc_gene, cmc_gene_gr)
    print("Done with CMC (Fromer et al.) gene eQTLs")
    
    # CMC (Fromer et al.) isoform eQTLs: ieSNPs
    cmc_iso <- data.table::fread(paste0(cmc, "/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt"))
    cmc_iso <- cmc_iso[cmc_iso$FDR < 0.05, ]  # filter for FDR <= 0.05
    cmc_iso_gr <- GenomicRanges::GRanges(seqnames = cmc_iso$SNP_location,
                                     snp_id = cmc_iso$SNP)  # generate GRanges object 
    cmc_iso_gr <- cmc_iso_gr[!duplicated(cmc_iso_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
    GenomeInfoDb::seqlevelsStyle(cmc_iso_gr) <- "NCBI"
    seqlevels(cmc_iso_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cmc_iso_gr), c("X", "Y", "M"))
    saveRDS(cmc_iso_gr, paste0(outpath, "/CMC_adjustedSVA_isoQTL_SNPs.rds"))
    rm(cmc_iso, cmc_iso_gr)
    print("Done with CMC (Fromer et al.) isoform eQTLs")


    ############
    ### GTEx ###
    ############
    ### eQTL SNPs (eSNPs) ###
    # GTEx v8 is based on the new human genome GRCh38, as we are working on GRCH37/hg19 we need to downlift to this genome
    # Load hg38 to hg19 chain (chain downloaded from UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)
    hg38_to_hg19 <- rtracklayer::import.chain("/storage/groups/ccm01/datasets/2019_GTEx_v8/hg38ToHg19.over.chain")

    # GTEx Whole Blood
    gtex_wb <- data.table::fread(paste0(gtex, "/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt"))
    gtex_wb <- gtex_wb[gtex_wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
    gtex_wb_gr <- GenomicRanges::GRanges(seqnames = gtex_wb$chr,
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_wb$variant_pos)), 
                                                                   end = as.numeric(as.character(gtex_wb$variant_pos))),
                                         snp_id = gtex_wb$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_wb_gr) <- "UCSC"  # UCSC needed for the lift over
    gtex_wb_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_wb_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 12,360 => 12,292
    gtex_wb_gr_hg19 <- gtex_wb_gr_hg19[!duplicated(gtex_wb_gr_hg19)]  # no duplicated rows found: 12,292 => 11,785 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_wb_gr_hg19) <- "NCBI"
    seqlevels(gtex_wb_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_wb_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_wb_gr_hg19, paste0(outpath, "/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
    rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)
    print("Done with GTEx whole blood eQTLs")
    
    # GTEx Brain Prefrontal Cortex BA9
    gtex_ba9 <- data.table::fread(paste0(gtex, "/GTEx_Analysis_v8_eQTL/Brain_Frontal_Cortex_BA9.v8.egenes.txt"))
    gtex_ba9 <- gtex_ba9[gtex_ba9$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
    gtex_ba9_gr <- GenomicRanges::GRanges(seqnames = gtex_ba9$chr,
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_ba9$variant_pos)), 
                                                                   end = as.numeric(as.character(gtex_ba9$variant_pos))),
                                         snp_id = gtex_ba9$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_ba9_gr) <- "UCSC"  # UCSC needed for the lift over
    gtex_ba9_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_ba9_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19)
    gtex_ba9_gr_hg19 <- gtex_ba9_gr_hg19[!duplicated(gtex_ba9_gr_hg19)]  # no duplicated rows found: 12,292 => 11,785 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_ba9_gr_hg19) <- "NCBI"
    seqlevels(gtex_ba9_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_ba9_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_ba9_gr_hg19, paste0(outpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_eQTL_SNPs.rds"))
    rm(gtex_ba9, gtex_ba9_gr, gtex_ba9_gr_hg19)
    print("Done with GTEx brain frontal cortex BA9 eQTLs")
    
    # GTEx Brain (all brain tissues)
    brains <- c("Brain_Amygdala.v8.egenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v8.egenes.txt", "Brain_Caudate_basal_ganglia.v8.egenes.txt", "Brain_Cerebellar_Hemisphere.v8.egenes.txt", "Brain_Cerebellum.v8.egenes.txt", "Brain_Cortex.v8.egenes.txt", "Brain_Frontal_Cortex_BA9.v8.egenes.txt", "Brain_Hippocampus.v8.egenes.txt", "Brain_Hypothalamus.v8.egenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v8.egenes.txt", "Brain_Putamen_basal_ganglia.v8.egenes.txt", "Brain_Spinal_cord_cervical_c-1.v8.egenes.txt", "Brain_Substantia_nigra.v8.egenes.txt")
    gtex_brains <- data.frame()
    for (i in brains) {
        filename <- paste0(gtex, "/GTEx_Analysis_v8_eQTL/", i)
        data <- data.table::fread(filename)
        data <- data[data$qval < 0.05, ]  # filter for significant results
        gtex_brains <- rbind(gtex_brains, data) # merge different brain files together
    }
    gtex_brains_gr <- GenomicRanges::GRanges(seqnames = gtex_brains$chr, 
                                             ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_brains$variant_pos)), 
                                                                       end = as.numeric(as.character(gtex_brains$variant_pos))),
                                             snp_id = gtex_brains$rs_id_dbSNP151_GRCh38p7,
                                             state = gtex_brains)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_brains_gr) <- "UCSC"
    gtex_brains_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_brains_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 89,312 => 88,885
    gtex_brains_gr_hg19 <- gtex_brains_gr_hg19[!duplicated(gtex_brains_gr_hg19)]  # no duplicated rows found: 88,885 => 61,935 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_brains_gr_hg19) <- "NCBI"
    seqlevels(gtex_brains_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_brains_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_brains_gr_hg19, paste0(outpath, "/GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
    rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19)
    print("Done with GTEx brain (all brain tissues) eQTLs")
    
    ### splice QTL (sQTL) SNPs (seSNPs) ###
    # GTEx Whole Blood
    gtex_wb <- data.table::fread(paste0(gtex, "/GTEx_Analysis_v8_sQTL/Whole_Blood.v8.sgenes.txt"))
    gtex_wb <- gtex_wb[gtex_wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
    gtex_wb_gr <- GenomicRanges::GRanges(seqnames = gtex_wb$chr,
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_wb$variant_pos)),
                                                                   end = as.numeric(as.character(gtex_wb$variant_pos))),
                                         snp_id = gtex_wb$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_wb_gr) <- "UCSC"  # UCSC needed for the lift over
    gtex_wb_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_wb_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 3,013 => 3,004
    gtex_wb_gr_hg19 <- gtex_wb_gr_hg19[!duplicated(gtex_wb_gr_hg19)]  # no duplicated rows found: 3,004 => 2,837 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_wb_gr_hg19) <- "NCBI"
    seqlevels(gtex_wb_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_wb_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_wb_gr_hg19, paste0(outpath, "/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))
    rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)
    print("Done with GTEx whole blood sQTLs")
    
    # GTEx Brain Prefrontal Cortex BA9
    gtex_ba9 <- data.table::fread(paste0(gtex, "/GTEx_Analysis_v8_sQTL/Brain_Frontal_Cortex_BA9.v8.sgenes.txt"))
    gtex_ba9 <- gtex_ba9[gtex_ba9$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
    gtex_ba9_gr <- GenomicRanges::GRanges(seqnames = gtex_ba9$chr,
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_ba9$variant_pos)),
                                                                   end = as.numeric(as.character(gtex_ba9$variant_pos))),
                                         snp_id = gtex_ba9$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_ba9_gr) <- "UCSC"  # UCSC needed for the lift over
    gtex_ba9_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_ba9_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 3,013 => 3,004
    gtex_ba9_gr_hg19 <- gtex_ba9_gr_hg19[!duplicated(gtex_ba9_gr_hg19)]  # no duplicated rows found: 3,004 => 2,837 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_ba9_gr_hg19) <- "NCBI"
    seqlevels(gtex_ba9_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_ba9_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_ba9_gr_hg19, paste0(outpath, "/GTEx_v8_Brain_Frontal_Cortex_BA9_hg19_sQTL_SNPs.rds"))
    rm(gtex_ba9, gtex_ba9_gr, gtex_ba9_gr_hg19)
    print("Done with GTEx brain frontal cortex BA9 sQTLs")

    # GTEx Brain (all brain tissues)
    brains <- c("Brain_Amygdala.v8.sgenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v8.sgenes.txt", "Brain_Caudate_basal_ganglia.v8.sgenes.txt", "Brain_Cerebellar_Hemisphere.v8.sgenes.txt", "Brain_Cerebellum.v8.sgenes.txt", "Brain_Cortex.v8.sgenes.txt", "Brain_Frontal_Cortex_BA9.v8.sgenes.txt", "Brain_Hippocampus.v8.sgenes.txt", "Brain_Hypothalamus.v8.sgenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v8.sgenes.txt", "Brain_Putamen_basal_ganglia.v8.sgenes.txt", "Brain_Spinal_cord_cervical_c-1.v8.sgenes.txt", "Brain_Substantia_nigra.v8.sgenes.txt")
    gtex_brains <- data.frame()
    for (i in brains) {
        filename <- paste0(gtex, "/GTEx_Analysis_v8_sQTL/", i)
        data <- data.table::fread(filename)
        data <- data[data$qval < 0.05, ]  # filter for significant results
        gtex_brains <- rbind(gtex_brains, data) # merge different brain files together
    }
    gtex_brains_gr <- GenomicRanges::GRanges(seqnames = gtex_brains$chr,
                                             ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_brains$variant_pos)),
                                                                       end = as.numeric(as.character(gtex_brains$variant_pos))),
                                             snp_id = gtex_brains$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
    GenomeInfoDb::seqlevelsStyle(gtex_brains_gr) <- "UCSC"
    gtex_brains_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_brains_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 20,642 => 20,527
    gtex_brains_gr_hg19 <- gtex_brains_gr_hg19[!duplicated(gtex_brains_gr_hg19)]  # no duplicated rows found: 20,527 => 11,804 SNPs
    GenomeInfoDb::seqlevelsStyle(gtex_brains_gr_hg19) <- "NCBI"
    seqlevels(gtex_brains_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_brains_gr_hg19), c("X", "Y", "M"))
    saveRDS(gtex_brains_gr_hg19, paste0(outpath, "/GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))
    rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19)
    print("Done with GTEx brain (all brain tissues) eQTLs")

    
    ####################
    ### Cortex meQTL ###
    ####################
    cortex_meqtl <- fread(cortex)
    cortex_meqtl <- cortex_meqtl[cortex_meqtl$p < 0.05, ]
    cortex_meqtl_gr <- GRanges(seqnames = cortex_meqtl$snpChr, 
                               ranges = IRanges(start = as.numeric(as.character(cortex_meqtl$snpPos)), 
                                                end = as.numeric(as.character(cortex_meqtl$snpPos))),
                               snp_id = cortex_meqtl$snpRsNum)  # generate GRanges object only with gene ID, name and location
    cortex_meqtl_gr <- cortex_meqtl_gr[!duplicated(cortex_meqtl_gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
    GenomeInfoDb::seqlevelsStyle(cortex_meqtl_gr) <- "NCBI"
    seqlevels(cortex_meqtl_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cortex_meqtl_gr), c("X", "Y", "M"))
    saveRDS(cortex_meqtl_gr, paste0(outpath, "/Cortex_meQTL_SNPs.rds"))
    rm(cortex_meqtl, cortex_meqtl_gr)
    print("Done with cortex meQTLs")

    
    ######################################
    ### Variant Effect Predictor (VEP) ###
    ######################################
    vep_maf01 <- data.table::fread(vep, data.table = F)
    vep_maf01$SNP.Position <- as.numeric(sapply(as.character(vep_maf01$SNP), function(x) unlist(strsplit(x, "_"))[2]))
    vep_maf01$SNP.Chr <- as.character(sapply(as.character(vep_maf01$SNP), function(x) unlist(strsplit(x, "_"))[1]))
    vep_melt <- melt(vep_maf01, id.vars = c("SNP", "SNP.Chr", "SNP.Position"))
    vep_gr <- makeGRangesFromDataFrame(
        data.frame(
            chr = vep_melt$SNP.Chr, 
            start = vep_melt$SNP.Position, 
            end = vep_melt$SNP.Position, 
            code = as.character(vep_melt$variable)), 
        keep.extra.columns = TRUE) 

    # downlift from GRCh38 to GRCh37 (hg19)
    vep_gr_hg19 <- unlist(rtracklayer::liftOver(vep_gr, hg38_to_hg19)) # 10,719,405
    vep_gr_hg19 <- vep_gr_hg19[!duplicated(vep_gr_hg19)]  # remove duplicated rows
    GenomeInfoDb::seqlevelsStyle(vep_gr_hg19) <- "NCBI"
    seqlevels(vep_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(vep_gr_hg19), c("X", "Y", "M"))
    saveRDS(vep_gr_hg19, paste0(outpath, "/VEP_hg19_collapsed_MAF01.rds"))
    rm(vep_maf01, vep_melt, vep_gr, hg38_to_hg19, vep_gr_hg19)
    print("Done with VEP")
    

    ################
    ### ChromHMM ###
    ################
    # read data
    hmm <- readRDS(paste0(chromhmm, "/chromHMM_all_states.Rds"))
    # remove chromosome X, Y and M
    seqlevels(hmm, pruning.mode = "coarse") <- setdiff(seqlevels(hmm), c("chrX", "chrY", "chrM"))
    # load ChromHMM data information
    hmm_info <- read.csv2(paste0(chromhmm, "/epigenomes.tsv"), header = T, stringsAsFactors = F, sep = "\t")
    
    # get ChromHMM blood data
    # Blood has multiple  cell lines (filter for Group == "Blood & T-cell")
    blood <- hmm_info[hmm_info$GROUP == "Blood & T-cell", "EID"]  # get IDs for all blood related ChromHMM data
    hmm_blood <- hmm[hmm$code %in% blood, ]  # get blood data
    seqlevelsStyle(hmm_blood) <- "NCBI"
    saveRDS(hmm_blood, paste0(outpath, "/ChromHMM_Blood_and_T-cells.rds"))
    rm(blood, hmm_blood)
    print("Done with ChromHMM Blood & T-cells")

    # get ChromHMM brain data
    # Brain has multiple  cell lines (filter for Group == "Brain")
    brain <- hmm_info[hmm_info$GROUP == "Brain", "EID"]  # get IDs for all brain related ChromHMM data
    hmm_brain <- hmm[hmm$code %in% brain, ]  # get brain data
    seqlevelsStyle(hmm_brain) <- "NCBI"
    saveRDS(hmm_brain, paste0(outpath, "/ChromHMM_all_Brains.rds"))
    rm(brain, hmm_brain) 
    print("Done with ChromHMM all Brains")
    
    # get ChromHMM Brain DLPFC (Dorsolateral Prefrontal Cortex, EID = E073)
    hmm_dlpfc <- hmm[hmm$code == "E073"]
    seqlevelsStyle(hmm_dlpfc) <- "NCBI"
    saveRDS(hmm_dlpfc, paste0(outpath, "/ChromHMM_Brain_Dorsolateral_Prefrontal_Cortex.rds"))
    rm(hmm, hmm_info, hmm_dlpfc) 
    print("Done with ChromHMM Brain DLPFC")
    
    ###############################
    ### GWAS Summary Statistics ###
    ###############################
    ### AD-2019 ###
    ad <- data.table::fread(paste0(gwas, "/AD-2019/AD_sumstats_Jansenetal_2019sept.txt"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(ad, ofile = "GWAS_sumstats_AD-2019_p0.05.rds")
    rm(ad)

    ### ADHD-2019 ###
    # data from 2017, but Paper published 2019
    adhd <- data.table::fread(paste0(gwas, "/ADHD-2019/adhd_jul2017"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(adhd, ofile = "GWAS_sumstats_ADHD-2019_p0.05.rds")
    rm(adhd)

    ### ALCDEP-2018 ###
    alcdep <- data.table::fread(paste0(gwas, "/ALCDEP-2018/pgc_alcdep.trans_fe_unrel_geno.aug2018_release.txt"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(alcdep, ofile = "GWAS_sumstats_ALCDEP-2018_p0.05.rds")
    rm(alcdep)

    ### AN-2019 ###
    an <- data.table::fread(paste0(gwas, "/AN-2019/pgcAN2.2019-07.vcf.tsv"), header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"), skip = 70)
    colnames(an) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(an, ofile = "GWAS_sumstats_AN-2019_p0.05.rds")
    rm(an)

    ### ASD-2019 ###
    # data from 2017 but Paper published 2019, that is why we call it ASD-2019 not ASD-2017!
    asd <- data.table::fread(paste0(gwas, "/ASD-2019/iPSYCH-PGC_ASD_Nov2017"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(asd, ofile = "GWAS_sumstats_ASD-2019_p0.05.rds")
    rm(asd)

    ### AUDIT-2019 ###
    # data from 2018, but paper published 2019
    # positions were not available, so load the version including bp positions
    audit <- data.table::fread(paste0(gwas, "/AUDIT-2019/AUDIT_UKB_2018_AJP_with_positions.txt"), header = T, stringsAsFactors = F, select = c("chr", "bp", "rsid", "p_T"))
    colnames(audit) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(audit, ofile = "GWAS_sumstats_AUDIT-2019_p0.05.rds")
    rm(audit)

    ### BIP-2019 ###
    # data from 2018, but published 2019
    bip <- data.table::fread(paste0(gwas, "/BIP-2019/daner_PGC_BIP32b_mds7a_0416a"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(bip, ofile = "GWAS_sumstats_BIP-2019_p0.05.rds")
    rm(bip)
    
    ### BIPvsSCZ-2018 ###
    bip <- data.table::fread(paste0(gwas, "/BIP+SCZ-2018/BDvsSCZ.sumstats"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(bip, ofile = "GWAS_sumstats_BIPvsSCZ-2018_p0.05.rds")
    rm(bip)
    
    ### BIPSCZvsCONT-2018 ###
    bip <- data.table::fread(paste0(gwas, "/BIP+SCZ-2018/BDSCZvsCONT.sumstats"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(bip, ofile = "GWAS_sumstats_BIPSCZvsCONT-2018_p0.05.rds")
    rm(bip)

    ### BMI-2015 (ENGAGE) ###
    bmi <- data.table::fread(paste0(gwas, "/BMI-2015_ENGAGE/ENGAGE1000G_BMI.txt"), header = T, stringsAsFactors = F, select = c("chromosome", "position", "SNP", "P"))
    colnames(bmi) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(bmi, ofile = "GWAS_sumstats_BMI-2015_ENGAGE_p0.05.rds")
    rm(bmi)
    
    ### BMI-2015 (GIANT) ###
    # add positions using HapMap II release 22 (done already)
    bmi <- data.table::fread(paste0(gwas, "/BMI-2015_GIANT/SNP_gwas_mc_merge_nogc_with_locations.tbl.uniq"), header = T, stringsAsFactors = F, select = c("Chromosome", "Position", "SNP", "p"))
    colnames(bmi) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(bmi, ofile = "GWAS_sumstats_BMI-2015_GIANT_p0.05.rds")
    rm(bmi)
    
    ### BMI-2018 (GIANT) ###
    bmi <- data.table::fread(paste0(gwas, "/BMI-2018/Meta-analysis_Locke_et_al+UKBiobank_2018_UPDATED.txt"), header = T, stringsAsFactors = F, select = c("CHR", "POS", "SNP", "P"))
    colnames(bmi) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(bmi, ofile = "GWAS_sumstats_BMI-2018_p0.05.rds")
    rm(bmi)
    
    ### BMI-2019 (GIANT) ###
    bmi <- data.table::fread(paste0(gwas, "/BMI-2019/Bmi.giant-ukbb.meta-analysis.combined.23May2018.HapMap2_only.txt"), header = T, stringsAsFactors = F, select = c("CHR", "POS", "SNP", "P"))
    colnames(bmi) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(bmi, ofile = "GWAS_sumstats_BMI-2019_p0.05.rds")
    rm(bmi)

    ### CAD-2017 ###
    # downlift from hg38 to hg19 done already just take the hg19 version
    cad <- data.table::fread(paste0(gwas, "/CAD-2017/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.hg19.OR.txt"), header = T, stringsAsFactors = F, select = c("chr", "bp_hg19", "snptestid", "p-value_gc"))
    colnames(cad) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(cad, ofile = "GWAS_sumstats_CAD-2017_p0.05.rds")
    rm(cad)

    ### Cannabis-2018 ###
    cannabis <- data.table::fread(paste0(gwas, "/Cannabis-2018/Cannabis_ICC_UKB_het.txt"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    cannabis <- cannabis[complete.cases(cannabis)] # remove rows with missing information 
    GenerateGrangesObject(cannabis, ofile = "GWAS_sumstats_Cannabis-2018_p0.05.rds")
    rm(cannabis)

    ### CDG1-2013 ###
    # uplift from hg18 to hg19 done already, just load the hg19 version
    cdg1 <- data.table::fread(paste0(gwas, "/CDG1-2013/pgc.cross.full.2013-03.hg19.txt"), header = T, stringsAsFactors = F, select = c("hg18chr", "bp", "snpid", "pval"))
    colnames(cdg1) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(cdg1, ofile = "GWAS_sumstats_CDG1-2013_p0.05.rds")
    rm(cdg1)
    
    ### CDG2-2019 ###
    cdg2 <- data.table::fread(paste0(gwas, "/CDG2-2019/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt"), header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"))
    colnames(cdg2) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(cdg2, ofile = "GWAS_sumstats_CDG2-2019_p0.05.rds")
    rm(cdg2)

    ### EA-2018 ###
    ea <- data.table::fread(paste0(gwas, "/EA-2018/GWAS_EA_excl23andMe.txt"), header = T, stringsAsFactors = F, select = c("CHR", "POS", "MarkerName", "Pval"))
    colnames(ea) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(ea, ofile = "GWAS_sumstats_EA-2018_p0.05.rds")
    rm(ea)

    ### ED-2017 ###
    ed <- data.table::fread(paste0(gwas, "/ED-2017/pgc.ed.freeze1.summarystatistics.July2017.txt"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(ed, ofile = "GWAS_sumstats_ED-2017_p0.05.rds")
    rm(ed)

    ### Height-2014 ###
    # uplift from hg18 to hg19 done already, just load the hg19 version
    height <- data.table::fread(paste0(gwas, "/Height-2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq_hg19_with_locations.txt"), header = T, stringsAsFactors = F, select = c("Chromosome", "Position", "MarkerName", "p"))
    colnames(height) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(height, ofile = "GWAS_sumstats_HE-2014_p0.05.rds")
    rm(height)

    ### Intelligence-2018 ###
    intelligence <- data.table::fread(paste0(gwas, "/Intelligence-2018/SavageJansen_2018_intelligence_metaanalysis.txt"), header = T, stringsAsFactors = F, select = c("CHR", "POS", "SNP", "P"))
    colnames(intelligence) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(intelligence, ofile = "GWAS_sumstats_IN-2018_p0.05.rds")
    rm(intelligence)

    ### MDD1-2012 ###
    # uplift from hg18 to hg19 done already, just load the hg19 version
    mdd1 <- data.table::fread(paste0(gwas, "/MDD1-2012/pgc.mdd.full.2012-04.hg19.txt"), header = T, stringsAsFactors = F, select = c("hg18chr", "bp", "snpid", "pval"))
    colnames(mdd1) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(mdd1, ofile = "GWAS_sumstats_MDD1-2012_p0.05.rds")
    rm(mdd1)

    ### MDD2-2018 ###
    mdd2 <- data.table::fread(paste0(gwas, "/MDD2-2018/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(mdd2, ofile = "GWAS_sumstats_MDD2-2018_p0.05.rds")
    rm(mdd2)

    ### MDD3-2019 ###
    mdd3 <- data.table::fread(paste0(gwas, "/MDD3-2019/daner_pgc_mdd_meta_w2_no23andMe_rmUKBB"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(mdd3, ofile = "GWAS_sumstats_MDD3-2019_p0.05.rds")
    rm(mdd3)

    ### OCD-2017 ###
    ocd <- data.table::fread(paste0(gwas, "/OCD-2017/ocd_aug2017"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(ocd, ofile = "GWAS_sumstats_OCD-2017_p0.05.rds")
    rm(ocd)

    ### PTSD-2019 ###
    ptsd <- data.table::fread(paste0(gwas, "/PTSD-2019/pts_all_freeze2_overall.results"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(ptsd, ofile = "GWAS_sumstats_PTSD-2019_p0.05.rds")
    rm(ptsd)

    ### SCZ2-2014 ###
    scz2 <- data.table::fread(paste0(gwas, "/SCZ2-2014/scz2.snp.results.txt"), header = T, stringsAsFactors = F, select = c("hg19chrc", "bp", "snpid", "p"))
    colnames(scz2) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(scz2, ofile = "GWAS_sumstats_SCZ2-2014_p0.05.rds")
    rm(scz2)

    ### T2D-2017 ###
    t2d <- data.table::fread(paste0(gwas, "/T2D-2017/METAANALYSIS_DIAGRAM_SE1_BP_POS.txt"), header = T, stringsAsFactors = F, select = c("Chr", "Position", "Chr:Position", "P-value"))
    colnames(t2d) <- c("CHR", "BP", "SNP", "P")
    GenerateGrangesObject(t2d, ofile = "GWAS_sumstats_T2D-2017_p0.05.rds")
    rm(t2d)

    ### TS-2018 ###
    ts <- data.table::fread(paste0(gwas, "/TS-2018/TS_Oct2018"), header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
    GenerateGrangesObject(ts, ofile = "GWAS_sumstats_TS-2018_p0.05.rds")
    rm(ts)

    print(sessionInfo())
}
