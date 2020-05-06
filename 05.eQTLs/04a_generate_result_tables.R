####################################################################
### Generate eQTL result tables from the MatrixEQTL output files ###
####################################################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args inpath=<path/input-file> snps=<path/snp-file> gtf=<path/gtf-file> gpos=<path/genespos-file> tpos=<path/transpos-file> epos=<path/exonspos-file> outpath=<path/for/output> 04a_generate_result_table.R
# or using Rscript: Rscript --slave 04a_generate_result_table.R inpath=<path/input-file> snps=<path/snp-file> gtf=<path/gtf-file> gpos=<path/genespos-file> tpos=<path/transpos-file> epos=<path/exonspos-file> outpath=<path/for/output>

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)
library(rtracklayer)  # to load gtf file
library(tidyr)  # to separate IDs

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly=TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate eQTL result tables from the MatrixEQTL output files.\n\nYou can invoke the script as follows:\nRscript 04a_generate_result_tables.R inpath=\"/home/.../to/input-files\" snps=\"/home/.../SNP.file\" gtf=\"/home/.../gencode.gtf\" gpos=\"/home/.../gene-positions.file\" tpos=\"/home/.../transcript-positions.file\" epos=\"/home/.../exons-positions.file\" outpath=\"/home/.../for/output-files\" \n\nCommandline parameters:\n\ninpath\t\tpath, where the input files are stored\nsnps\t\tpath and file name of the SNP genotype data\ngtf\t\tpath and file name of the gencode .gtf file\ngpos\t\tpath and file name of the positions for all genes\ntpos\t\tpath and file name of the positions for all transcripts\nepos\t\tpath and file name of the positions for all exons\noutpath\t\tpath to store the output files\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  # 1. Load all matrixEQTL output files and filter them for significant hits (FDR < 0.05)
  # 1.1 Gene cis-eQTL
  # 1.1.1 Cases and controls
  data <- fread(paste0(inpath, "/01_pmb_CASE_CTL_matrixEQTL_cis_output_genes.txt"), header = T, stringsAsFactors = F)
  cis_genes_all <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_genes_all))
  # 1.1.2 Cases only
  data <- fread(paste0(inpath, "/02_pmb_CASE_only_matrixEQTL_cis_output_genes.txt"), header = T, stringsAsFactors = F)
  cis_genes_cases <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_genes_cases))
  # 1.1.3 Controls (CTL) only
  data <- fread(paste0(inpath, "/03_pmb_CTL_only_matrixEQTL_cis_output_genes.txt"), header = T, stringsAsFactors = F)
  cis_genes_ctl <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_genes_ctl))
  # 1.2 Transcript cis-eQTL
  # 1.2.1 Cases and controls
  data <- fread(paste0(inpath, "/01_pmb_CASE_CTL_matrixEQTL_cis_output_transcripts.txt"), header = T, stringsAsFactors = F)
  cis_transcripts_all <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_transcripts_all))
  # 1.2.2 Cases only
  data <- fread(paste0(inpath, "/02_pmb_CASE_only_matrixEQTL_cis_output_transcripts.txt"), header = T, stringsAsFactors = F)
  cis_transcripts_cases <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_transcripts_cases))
  # 1.2.3 Controls (CTL) only
  data <- fread(paste0(inpath, "/03_pmb_CTL_only_matrixEQTL_cis_output_transcripts.txt"), header = T, stringsAsFactors = F)
  cis_transcripts_ctl <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
  rm(data)
  print(dim(cis_transcripts_ctl))  
  # 1.3 Exon cis-eQTL
  # 1.3.1 Cases and controls
  cis_exons_all <- data.table("SNP" = character(), "gene" = character(), "beta" = numeric(), "t-stat" = numeric(), "p-value" = numeric(), "FDR" = numeric())
  chunk <- 1:5
  for (value in chunk) {
    # read chunk file
    data <- fread(paste0(inpath, "/01_pmb_CASE_CTL_matrixEQTL_cis_output_exons_chunk", value, ".txt"), header = T, stringsAsFactors = F)
    exon_chunk <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
    rm(data)
    cis_exons_all <- rbind(cis_exons_all, exon_chunk)  # add chunk to eqtl output dataset
  }  
  print(dim(cis_exons_all))
  # 1.3.2 Cases only
  cis_exons_cases <- data.table("SNP" = character(), "gene" = character(), "beta" = numeric(), "t-stat" = numeric(), "p-value" = numeric(), "FDR" = numeric())
  chunk <- 1:5
  for (value in chunk) {
    # read chunk file
    data <- fread(paste0(inpath, "/02_pmb_CASE_only_matrixEQTL_cis_output_exons_chunk", value, ".txt"), header = T, stringsAsFactors = F)
    exon_chunk <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
    rm(data)
    cis_exons_cases <- rbind(cis_exons_cases, exon_chunk)  # add chunk to eqtl output dataset
  }  
  print(dim(cis_exons_cases))
  # 1.3.3 Controls (CTL) only
  cis_exons_ctl <- data.table("SNP" = character(), "gene" = character(), "beta" = numeric(), "t-stat" = numeric(), "p-value" = numeric(), "FDR" = numeric())
  chunk <- 1:5
  for (value in chunk) {
    # read chunk file
    data <- fread(paste0(inpath, "/03_pmb_CTL_only_matrixEQTL_cis_output_exons_chunk", value, ".txt"), header = T, stringsAsFactors = F)
    exon_chunk <- data[data$FDR < 0.05, ]  # filter for only eQTLs with adjusted p-value (FDR) < 0.05
    rm(data)
    cis_exons_ctl <- rbind(cis_exons_ctl, exon_chunk)  # add chunk to eqtl output dataset
  }  
  print(dim(cis_exons_ctl))
  print("1. Done loading eQTL data")
  
  # 2. Load SNP, gene, transcript and exon positions
  # 2.1 Load SNP positions
  snpspos <- fread(snps, header = T, stringsAsFactors = F, data.table = F) # snp chr pos
  colnames(snpspos) <- c("chromosome", "SNP", "distance", "position", "SNP_allele1", "SNP_allele2")
  snpspos$SNP_location <- paste0(snpspos$chromosome, ":", snpspos$position)
  snpspos <- snpspos[, c("SNP_location", "SNP", "SNP_allele1", "SNP_allele2")]
  snpspos <- snpspos[!duplicated(snpspos$SNP), ]
  # 2.2 Load the GENCODE file and get informations for genes, transcripts and exons
  gencode <- import.gff(gtf, format = "gtf")
  gencode <- gencode[gencode$gene_type == "protein_coding"]  # get only protein-coding genes, transcripts and exons
  gencode <- gencode[gencode$source == "HAVANA"]  # get only data from HAVANA source
  # 2.2.1 Gene positions
  genes <- as.data.frame(gencode[gencode$type == "gene"])  # get only genes (no exons, cds, transcripts etc.)
  genes <- genes[, c("gene_id", "gene_name", "strand")]
  genespos <- fread(gpos, header = T, stringsAsFactors = F, data.table = F)
  genespos$ids <- genespos$geneid
  genespos <- separate(genespos, ids, c("firstid", "restids"), sep = "\\|", extra = "merge", fill = "right")
  genespos <- merge(genespos, genes, by.x = "firstid", by.y = "gene_id", all.x = T)
  genespos$Gene_location <- paste0(genespos$chr, ":", genespos$start, "..", genespos$end)
  genespos <- genespos[, c("geneid", "gene_name", "Gene_location", "strand")]
  colnames(genespos) <- c("Gene", "Gene_symbol", "Gene_location", "Gene_strand")
  rm(genes)
  # 2.2.2 Transcript positions
  # For the transcripts we can not use the GENCODE information alone as we merged several transcript IDs and used the min and max position of these merged transcripts. 
  # So, we first need to load the transcript positions file (used for the eQTLs) and then merge it with the gencode information.
  transcripts <- as.data.frame(gencode[gencode$type == "transcript"])  # get only transcripts
  transcripts <- transcripts[, c("transcript_id", "strand", "gene_id", "gene_name")]
  transpos <- fread(tpos, header = T, stringsAsFactors = F, data.table = F)
  transpos$ids <- transpos$transcriptid
  transpos <- separate(transpos, ids, c("firstid", "restids"), sep = "\\|", extra = "merge", fill = "right")
  transpos <- merge(transpos, transcripts, by.x = "firstid", by.y = "transcript_id", all.x = T)
  transpos$Transcript_location <- paste(transpos$chr, ":", transpos$start, "..", transpos$end, sep = "")
  transpos <- transpos[, c("transcriptid", "Transcript_location", "gene_id", "gene_name")]
  colnames(transpos) <- c("Transcript", "Transcript_location", "Gene", "Gene_symbol")
  rm(transcripts)
  # 2.2.3 Exons
  exons <- as.data.frame(gencode[gencode$type == "exon"])  # get only exons
  exons <- exons[, c("exon_id", "strand", "transcript_id", "gene_id", "gene_name")]
  exonspos <- fread(epos, header = T, stringsAsFactors = F, data.table = F)
  exonspos$ids <- exonspos$exonid
  exonspos <- separate(exonspos, ids, c("firstid", "restids"), sep = "\\|", extra = "merge", fill = "right")
  exonspos <- merge(exonspos, exons, by.x = "firstid", by.y = "exon_id", all.x = T)
  merged_transcripts <- aggregate(transcript_id ~ firstid, exonspos, paste, collapse = '|')
  exonspos <- exonspos[, c("firstid", "exonid", "chr", "start", "end", "strand", "gene_id", "gene_name")]
  exonspos <- merge(merged_transcripts, exonspos)
  exonspos <- exonspos[!duplicated(exonspos),]
  exonspos$Exon_location <- paste(exonspos$chr, ":", exonspos$start, "..", exonspos$end, sep = "")
  exonspos <- exonspos[, c("exonid", "Exon_location", "transcript_id", "gene_id", "gene_name")]
  colnames(exonspos) <- c("Exon", "Exon_location", "Transcript", "Gene", "Gene_symbol")
  rm(gencode, exons)
  print("2. Done loading gene, transcript and exon positions")
  
  # 3. Generate result tables
  # We will add some more information to the significant eQTL hits (annotation) and save the data into result tables (.txt files)
  # 3.1 Gene cis-eQTLs
  # 3.1.1 Cases and controls
  eqtls_gene <- merge(cis_genes_all, snpspos, by = "SNP", all.x = T)  # add SNP positions
  eqtls_gene <- merge(eqtls_gene, genespos, by.x = "gene", by.y = "Gene", all.x = T)  # add gene positions
  eqtls_gene <- eqtls_gene[, c("SNP", "gene", "Gene_symbol", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
  eqtls_gene$eQTL_type <- "cis" 
  colnames(eqtls_gene) <- c("SNP ID", "Ensembl Gene ID", "Gene Symbol", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Gene Location",  "Gene Strand", "eQTL Type")
  eqtls_gene <- eqtls_gene[order(eqtls_gene$FDR), ]
  print(paste0("Gene eQTLs cases vs. controls: Rows: ", dim(eqtls_gene)[[1]], "  Columns: ", dim(eqtls_gene)[[2]]))
  write.table(eqtls_gene, file = paste0(outpath, "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.1.2 Cases only
  eqtls_gene_cases <- merge(cis_genes_cases, snpspos, by = "SNP", all.x = T)
  eqtls_gene_cases <- merge(eqtls_gene_cases, genespos, by.x = "gene", by.y = "Gene", all.x = T)
  eqtls_gene_cases[grep("\\|", eqtls_gene_cases$gene), c("Gene_symbol", "Gene_location", "Gene_strand")] <- genespos[genespos$Gene == strsplit(as.character(eqtls_gene_cases[grep("\\|", eqtls_gene_cases$gene), "gene"]), "\\|")[[1]][1], c("Gene_symbol", "Gene_location", "Gene_strand")]
  eqtls_gene_cases <- eqtls_gene_cases[complete.cases(eqtls_gene_cases), ]
  eqtls_gene_cases <- eqtls_gene_cases[, c("SNP", "gene", "Gene_symbol", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
  eqtls_gene_cases$eQTL_type <- "cis" 
  colnames(eqtls_gene_cases) <- c("SNP ID", "Ensembl Gene ID", "Gene Symbol", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Gene Location",  "Gene Strand", "eQTL Type")
  eqtls_gene_cases <- eqtls_gene_cases[order(eqtls_gene_cases$FDR), ]
  print(paste0("Gene eQTLs cases only: Rows: ", dim(eqtls_gene_cases)[[1]], "  Columns: ", dim(eqtls_gene_cases)[[2]]))
  write.table(eqtls_gene_cases, file = paste0(outpath, "/02_pmb_cis-eQTL_gene_hits_cases_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.1.3 Controls only
  eqtls_gene_ctl <- merge(cis_genes_ctl, snpspos, by = "SNP", all.x = T)
  eqtls_gene_ctl <- merge(eqtls_gene_ctl, genespos, by.x = "gene", by.y = "Gene", all.x = T)
  eqtls_gene_ctl[grep("\\|", eqtls_gene_ctl$gene), c("Gene_symbol", "Gene_location", "Gene_strand")] <- genespos[genespos$Gene == strsplit(as.character(eqtls_gene_ctl[grep("\\|", eqtls_gene_ctl$gene), "gene"]), "\\|")[[1]][1], c("Gene_symbol", "Gene_location", "Gene_strand")]
  eqtls_gene_ctl <- eqtls_gene_ctl[complete.cases(eqtls_gene_ctl), ]
  eqtls_gene_ctl <- eqtls_gene_ctl[, c("SNP", "gene", "Gene_symbol", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
  eqtls_gene_ctl$eQTL_type <- "cis" 
  colnames(eqtls_gene_ctl) <- c("SNP ID", "Ensembl Gene ID", "Gene Symbol", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Gene Location",  "Gene Strand", "eQTL Type")
  eqtls_gene_ctl <- eqtls_gene_ctl[order(eqtls_gene_ctl$FDR), ]
  print(paste0("Gene eQTLs controls only: Rows: ", dim(eqtls_gene_ctl)[[1]], "  Columns: ", dim(eqtls_gene_ctl)[[2]]))
  write.table(eqtls_gene_ctl, file = paste0(outpath, "/03_pmb_cis-eQTL_gene_hits_controls_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.2 Transcript cis-eQTLs
  # 3.2.1 Cases and controls
  eqtls_transcript <- merge(cis_transcripts_all, snpspos, by = "SNP", all.x = T)
  eqtls_transcript <- merge(eqtls_transcript, transpos, by.x = "gene", by.y = "Transcript", all.x = T)
  eqtls_transcript <- merge(eqtls_transcript, genespos, by = "Gene", all.x = T)
  eqtls_transcript <- eqtls_transcript[complete.cases(eqtls_transcript), ]
  eqtls_transcript <- eqtls_transcript[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Transcript_location", "Gene",  "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_transcript$eQTL_type <- "cis"  
  colnames(eqtls_transcript) <- c("SNP ID", "Ensembl Transcript ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Transcript Location", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_transcript <- eqtls_transcript[order(eqtls_transcript$FDR), ]
  print(paste0("Transcript eQTLs cases vs. controls: Rows: ", dim(eqtls_transcript)[[1]], "  Columns: ", dim(eqtls_transcript)[[2]]))
  write.table(eqtls_transcript, file = paste0(outpath, "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.2.2 Cases only
  eqtls_transcript_cases <- merge(cis_transcripts_cases, snpspos, by = "SNP", all.x = T)
  eqtls_transcript_cases <- merge(eqtls_transcript_cases, transpos, by.x = "gene", by.y = "Transcript", all.x = T)
  eqtls_transcript_cases <- merge(eqtls_transcript_cases, genespos, by = "Gene", all.x = T)
  eqtls_transcript_cases <- eqtls_transcript_cases[complete.cases(eqtls_transcript_cases), ]
  eqtls_transcript_cases <- eqtls_transcript_cases[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Transcript_location", "Gene", "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_transcript_cases$eQTL_type <- "cis"
  colnames(eqtls_transcript_cases) <- c("SNP ID", "Ensembl Transcript ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Transcript Location", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_transcript_cases <- eqtls_transcript_cases[order(eqtls_transcript_cases$FDR),]
  print(paste0("Transcript eQTLs cases only: Rows: ", dim(eqtls_transcript_cases)[[1]], "  Columns: ", dim(eqtls_transcript_cases)[[2]]))
  write.table(eqtls_transcript_cases, file = paste0(outpath, "/02_pmb_cis-eQTL_transcript_hits_cases_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.2.3 Controls only
  eqtls_transcript_ctl <- merge(cis_transcripts_ctl, snpspos, by = "SNP", all.x = T)
  eqtls_transcript_ctl <- merge(eqtls_transcript_ctl, transpos, by.x = "gene", by.y = "Transcript", all.x = T)
  eqtls_transcript_ctl <- merge(eqtls_transcript_ctl, genespos, by = "Gene", all.x = T)
  eqtls_transcript_ctl <- eqtls_transcript_ctl[complete.cases(eqtls_transcript_ctl), ]
  eqtls_transcript_ctl <- eqtls_transcript_ctl[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Transcript_location", "Gene", "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_transcript_ctl$eQTL_type <- "cis"
  colnames(eqtls_transcript_ctl) <- c("SNP ID", "Ensembl Transcript ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Transcript Location", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_transcript_ctl <- eqtls_transcript_ctl[order(eqtls_transcript_ctl$FDR),]
  print(paste0("Transcript eQTLs controls only: Rows: ", dim(eqtls_transcript_ctl)[[1]], "  Columns: ", dim(eqtls_transcript_ctl)[[2]]))
  write.table(eqtls_transcript_ctl, file = paste0(outpath, "/03_pmb_cis-eQTL_transcript_hits_controls_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.3 Exon cis-eQTLs
  # 3.3.1 Cases and controls
  eqtls_exon <- merge(cis_exons_all, snpspos, by = "SNP", all.x = T)
  eqtls_exon <- merge(eqtls_exon, exonspos, by.x = "gene", by.y = "Exon", all.x = T)
  eqtls_exon <- merge(eqtls_exon, genespos, by = "Gene", all.x = T)
  eqtls_exon <- eqtls_exon[complete.cases(eqtls_exon), ]
  eqtls_exon <- eqtls_exon[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Transcript", "Gene", "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_exon$eQTL_type <- "cis"  
  colnames(eqtls_exon) <- c("SNP ID", "Ensembl Exon ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Exon Location", "Ensembl Transcript ID", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_exon <- eqtls_exon[order(eqtls_exon$FDR), ]
  print(paste0("Exon eQTLs cases vs. controls: Rows: ", dim(eqtls_exon)[[1]], "  Columns: ", dim(eqtls_exon)[[2]]))
  write.table(eqtls_exon, file = paste0(outpath, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.3.2 Cases only
  eqtls_exon_cases <- merge(cis_exons_cases, snpspos, by = "SNP", all.x = T)
  eqtls_exon_cases <- merge(eqtls_exon_cases, exonspos, by.x = "gene", by.y = "Exon", all.x = T)
  eqtls_exon_cases <- merge(eqtls_exon_cases, genespos, by = "Gene", all.x = T)
  eqtls_exon_cases <- eqtls_exon_cases[complete.cases(eqtls_exon_cases), ]
  eqtls_exon_cases <- eqtls_exon_cases[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Transcript", "Gene", "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_exon_cases$eQTL_type <- "cis"   
  colnames(eqtls_exon_cases) <- c("SNP ID", "Ensembl Exon ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Exon Location", "Ensembl Transcript ID", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_exon_cases <- eqtls_exon_cases[order(eqtls_exon_cases$FDR), ]
  print(paste0("Exon eQTLs cases only: Rows: ", dim(eqtls_exon_cases)[[1]], "  Columns: ", dim(eqtls_exon_cases)[[2]]))
  write.table(eqtls_exon_cases, file = paste0(outpath, "/02_pmb_cis-eQTL_exon_hits_cases_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  # 3.3.3 Controls only
  eqtls_exon_ctl <- merge(cis_exons_ctl, snpspos, by = "SNP", all.x = T)
  eqtls_exon_ctl <- merge(eqtls_exon_ctl, exonspos, by.x = "gene", by.y = "Exon", all.x = T)
  eqtls_exon_ctl <- merge(eqtls_exon_ctl, genespos, by = "Gene", all.x = T)
  eqtls_exon_ctl <- eqtls_exon_ctl[complete.cases(eqtls_exon_ctl), ]
  eqtls_exon_ctl <- eqtls_exon_ctl[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Transcript", "Gene", "Gene_symbol.x", "Gene_location", "Gene_strand")]
  eqtls_exon_ctl$eQTL_type <- "cis"   
  colnames(eqtls_exon_ctl) <- c("SNP ID", "Ensembl Exon ID", "Beta coefficient", "T-statistic", "P-value", "FDR estimate", "SNP Location", "SNP Allele1", "SNP Allele2", "Exon Location", "Ensembl Transcript ID", "Ensembl Gene ID", "Gene Symbol", "Gene Location", "Gene Strand", "eQTL Type")
  eqtls_exon_ctl <- eqtls_exon_ctl[order(eqtls_exon_ctl$FDR), ]
  print(paste0("Exon eQTLs controls only: Rows: ", dim(eqtls_exon_ctl)[[1]], "  Columns: ", dim(eqtls_exon_ctl)[[2]]))
  write.table(eqtls_exon_ctl, file = paste0(outpath, "/03_pmb_cis-eQTL_exon_hits_controls_only.txt"), row.names = F, quote = F, dec = ".", sep = "\t")
  print("3. Done generating result tables")
}
