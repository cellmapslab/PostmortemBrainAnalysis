#############################
### Generate overlap file ###
#############################
# Here, we will generate the overlap file for the index-patient rare variant genes, exon eQTL genes (eeGenes) and ADHD, BIP, MDD and SCZ (ABMS) eQTS genes.
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args snvs=<path/to/snv.file> eqts=<path/to/eqts.file> eqtl=<path/to/etl.file> gencode=<path/to/gencode.file> outpath=<path/for/output>" 02_generate_gene_overlap_file.R &>> generate_gene_overlap_file.R
# or with Rscript: Rscript --slave 02_generate_gene_overlap_file.R snvs=<path/to/snv.file> eqts=<path/to/eqts.file> eqtl=<path/to/etl.file> gencode=<path/to/gencode.file> outpath=<path/for/output> &>> generate_gene_overlap_file.log

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)  # to load data with fread (faster when read.table etc.)
library(GenomicRanges)  # to load the GRanges object
library(rtracklayer)  # to read Gencodes .gtf file

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate the overlap file for the index-patient rare variant genes, exon eQTL genes (eeGenes) and ADHD, BIP, MDD and SCZ (ABMS) eQTS genes.\n\nYou can invoke the script as follows:\nRscript 02_generate_gene_overlap_file.R snvs=\"/home/.../to/snv.file\" eqts=\"/home/.../to/eqts.file/\" eqtl=\"/home/.../to/eqtl.file/\" gencode=\"path/to/gencode.gtf\" outpath=\"path/to/store/output\" &>> generate_gene_overlap_file.log \n\nCommandline parameters:\n\nsnvs\t\tpath and file name of rare mutation data files from index patients (IntegraMent, Bonn) \neqts\t\tpath and file name of ADHD, BIP, MDD and SCZ joint exon eQTS file \neqtl\t\tpath and name of exon eqtl result file \noutpath\t\tpath to store all output files into\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  # We want to generate the overlap between exon eQTL Genes (eeGenes), joint exon eQTS ABMS (ADHD + BIP + MDD + SCZ) 
  # genes and rare variant damaged genes from index-patients.
  
  #################
  ### Load Data ###
  #################
  # 1. Load ABMS eQTS
  abms <- fread(eqts, data.table = F, stringsAsFactors = F, header = F)
  print(paste0("# joint eeQTS genes: ", length(unique(abms$V1))))
  
  # 2. Load exon eQTL genes (eeGenes)
  eqtls_exon <- fread(eqtl, header = T, stringsAsFactors = F)
  eqtls_exon$`Ensembl Gene ID` <- as.vector(sapply(unlist(strsplit(eqtls_exon$`Ensembl Gene ID`, "\\|")), function(x) strsplit(x, "\\.")[[1]][1]))
  print(paste0("# eeGenes: ", length(unique(eqtls_exon$`Ensembl Gene ID`))))
  
  # 3. Load Index RData file
  index <- as.data.frame(readRDS(snvs))
  print(paste0("# rare variant genes: ", length(unique(index$Ensembl.Gene.ID))))
  
  #######################################################
  ### Generate gene symbol to Ensembl gene ID mapping ###
  #######################################################
  gencode_all <- import.gff(gencode, format = "gtf")  # Load GENCODE .gtf file for gene symbol to Ensembl Gene ID mapping
  genes <- as.data.frame(gencode_all[gencode_all$type == "gene"])  # get only genes
  symbol_to_id <- genes[, c("gene_id", "gene_name")]  # map gene symbol to Ensembl gene ID
  symbol_to_id <- symbol_to_id[!duplicated(symbol_to_id), ]  # remove duplicated rows
  symbol_to_id$gene_id <- sapply(symbol_to_id$gene_id, function(x) strsplit(x, "\\.")[[1]][1])  # remove .xx endings
  rm(gencode_all, genes)
  
  ######################
  ### Generate Files ###
  ######################
  # 1. ABMS eQTS file
  abms_new <- merge(abms, symbol_to_id, by.x = "V1", by.y = "gene_id")
  abms_new <- abms_new[, c("gene_name", "V1")]
  colnames(abms_new) <- c("Gene Symbol", "Ensembl Gene ID")
  dim(abms_new)
  write.table(abms_new, file = paste0(outpath, "/1_ABMS_eQTS_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  # 2. exon eQTL genes (eeGenes) file
  eeGenes <- eqtls_exon[, c("Gene Symbol", "Ensembl Gene ID")]
  eeGenes <- eeGenes[!duplicated(eeGenes), ]  # remove duplicates
  eeGenes$`Ensembl Gene ID` <- unlist(lapply(eeGenes$`Ensembl Gene ID`, function(x) strsplit(x, "\\.")[[1]][1]))
  dim(eeGenes)
  write.table(eeGenes, file = paste0(outpath, "/2_exon_eQTL_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  # 3. genes of rare mutations from index patients file
  index_new <- index[, c("GeneSymbol", "Ensembl.Gene.ID")]
  colnames(index_new) <- c("Gene Symbol", "Ensembl Gene ID")
  dim(index_new)
  write.table(index_new, file = paste0(outpath, "/3_rare_variant_index_genes.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  # 4. Overlapping genes file
  overlap <- intersect(eeGenes$`Ensembl Gene ID`, intersect(abms_new$`Ensembl Gene ID`, index_new$`Ensembl Gene ID`))
  length(overlap)  # 832
  overlap_new <- index_new[index_new$`Ensembl Gene ID` %in% overlap, ]
  overlap_new <- overlap_new[!duplicated(overlap_new), ]
  dim(overlap_new)  # 832
  write.table(overlap_new, file = paste0(outpath, "/4_overlap_eeGenes_ABMS_rare.txt"), sep = "\t", col.names = T, row.names = F, quote = F)
  
  ####################
  ### Session Info ###
  ####################
  print(sessionInfo())
}
