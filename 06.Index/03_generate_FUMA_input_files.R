###################################
### FUMA - Generate input files ###
###################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log
# or with Rscript: Rscript --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)
library(limma)
library(rtracklayer)
library(GenomicRanges)

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate the overlap file for the index-patient rare variant genes, exon eQTL genes (eeGenes) and ADHD, BIP, MDD and SCZ (ABMS) eQTS genes.\n\nYou can invoke the script as follows:\nRscript 02_generate_gene_overlap_file.R snvs=\"/home/.../to/snv.file\" eqts=\"/home/.../to/eqts.file/\" eqtl=\"/home/.../to/eqtl.file/\" gencode=\"path/to/gencode.gtf\" outpath=\"path/to/store/output\" &>> generate_gene_overlap_file.log \n\nCommandline parameters:\n\nsnvs\t\tpath and file name of rare mutation data files from index patients (IntegraMent, Bonn) \neqts\t\tpath and file name of ADHD, BIP, MDD and SCZ joint exon eQTS file \neqtl\t\tpath and name of exon eqtl result file \noutpath\t\tpath to store all output files into\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  #######################################
  ### Rare variants of index-patients ###
  #######################################
  index <- fread(paste0(inpath, "/3_rare_variant_index_genes.txt"), data.table = F, stringsAsFactors = F, header = T)  # load data
  index <- unique(index$`Ensembl Gene ID`)  # get unique Ensembl gene IDs
  print(paste0("# rare variant genes: ", length(index)))  # check number of genes
  write.table(index, file = paste0(outpath, "/1_index.txt"), row.names = F, col.names = F, quote = F)  # write FUMA input file
  
  #####################
  ### Overlap genes ###
  #####################
  # overlap of exon eQTL genes (eeGenes), ADHD, BIP, MDD and SCZ (ABMS) exon eQTS and rare variants of index patients
  overlap <- fread(paste0(inpath, "/4_overlap_eeGenes_ABMS_rare.txt"), data.table = F, stringsAsFactors = F, header = T)  # load data
  overlap <- unique(overlap$`Ensembl Gene ID`)  # get unique Ensembl gene IDs
  print(paste0("# genes from rare variants, eeQTL and eeQTS overlap: ", length(overlap)))  # check number of genes
  write.table(overlap, file = paste0(outpath, "/2_overlap.txt"), row.names = F, col.names = F, quote = F)  # write FUMA input file
  
  ############################################
  ### Load exon expression background file ###
  ############################################
  background_exon <- fread(ebkgd, data.table = F, stringsAsFactors = F, header = F)
  
  ##################################
  ### Load Index Background File ###
  ##################################
  # We will load the index background file, but here we need to get the gene symbols and map them to Ensembl gene IDs
  agilent <- fread(surselect, data.table = F, stringsAsFactors = F, header = F)  # read surselect background gene set
  agilent <- agilent[!duplicated(agilent), ]  # remove duplicated rows
  agilent$V6 <- unlist(lapply(agilent$V5, function(x) ifelse(length(x[grep("ref", x)]) != 0, strsplit(strsplit(x, ",")[[1]][1], "\\|")[[1]][2], "-")))  # get unique genes symbols
  agilent <- agilent[, c("V2", "V3", "V4", "V6")]
  colnames(agilent) <- c("seqnames", "start", "end", "symbol")
  dim(agilent)
  agilent <- agilent[!is.na(agilent$start), ] # some rows do not have a start or end positon, so we remove them
  dim(agilent)
  agilent_gr <- GenomicRanges::makeGRangesFromDataFrame(agilent, keep.extra.columns = T)
  length(agilent_gr)
  rm(agilent)
  
  # generate symbol to Ensembl gene ID mapping dataset
  gencode_all <- import.gff(gencode, format = "gtf")  # Load GENCODE .gtf file for gene symbol to Ensembl Gene ID mapping
  genes <- as.data.frame(gencode_all[gencode_all$type == "gene"])  # get only genes
  genes <- genes[, c("seqnames","start", "end", "strand", "gene_id", "gene_name")]
  genes_gr <- GenomicRanges::makeGRangesFromDataFrame(genes, keep.extra.columns = T) # generate a GRanges object
  rm(gencode_all, genes)
  
  # map gene symbols to Ensembl gene IDs
  background_genes <- mergeByOverlaps(agilent_gr, genes_gr)
  background_genes$gene_id <- sapply(background_genes$gene_id, function(x) strsplit(x, "\\.")[[1]][1])  # remove .xx endings from Ensembl gene IDs
  background_index <- unique(background_genes$gene_id)
 
  ##############################
  ### Overlap of Backgrounds ###
  ##############################
  background_overlap <- intersect(background_exon$V1, background_index)
  cat(paste0("# background exon eQTLs: ", length(background_exon$V1), "\n# background index: ", length(background_index), "\n# background overlap:", length(background_overlap)))
  
  ############################
  ### Save Background Data ###
  ############################
  write.table(background_index, file = paste0(outpath, "/1_index_bkgd.txt"), row.names = F, col.names = F, quote = F)  # write FUMA input file
  write.table(background_overlap, file = paste0(outpath, "/2_overlap_bkgd.txt"), row.names = F, col.names = F, quote = F)  # write FUMA input file
  
  ####################
  ### Session Info ###
  ####################
  print(sessionInfo()) 
}