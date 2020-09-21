###############################
### Generate GRanges object ###
###############################
# This R-Script generates a GRanges object for the rare mutations of the index-patient dataset.
# Invoke the script in the commandline as follows: R CMD BATCH --slave "--args inpath=<path/to/input-files> gencode=<path/to/gencode.file> outpath=<path/for/output>" 01_generate_Granges_object.R &>> generate_GRanges_object.R
# or using Rscript: Rscript --slave 01_generate_Granges_object.R inpath=<path/to/input-files> gencode=<path/to/gencode.file> outpath=<path/for/output> &>> generate_GRanges_object.log

# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(readxl)  # for reading in excel files
library(tidyr)  # for the separate_rows function
library(limma)  # for the mapping of gene symbol to Ensembl gene ID
library(rtracklayer)  # to read Gencodes .gtf file
library(GenomicRanges)  # to generate the GRanges object

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate a GRanges object for the rare mutations of the index-patient dataset.\n\nYou can invoke the script as follows:\nRscript --slave 01_generate_Granges_object.R inpath=<path/to/input-files> gencode=<path/to/gencode.file> outpath=<path/for/output> &>> generate_GRanges_object.log \n\nCommandline parameters:\n\ninpath\t\tpath to the IntegraMent (Bonn) rare variant dataset files from index-patients\ngencode\t\tpath and file name of Gencode .gtf file with all human GRCh37 (hg19) gene annotations \noutpath\t\tpath to store all output files into\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  
  ##########################
  ### Load rare variants ###
  ##########################
  # We found out, that we have a huge overlap to rare variants that are gene damaging from schizophrenia 
  # index-patients. 
  # Therefore, we first need to load the list from Lisa and from Frederik and filter it for following criterias: 
  # 1. 5 or 2 or INS/DEL/INDEL
  # 2. CADD score > 20
  # 3. ExAC < 0.1% (ExAC < 0.001) or NA

  # 1. load Lisa's data
  lisa <- as.data.frame(read_excel(paste0(inpath, "/2020-08-13_alleVariantenLisa.xlsx"), sheet = "Tabelle1"))
  lisa <- lisa[lisa$Anzahl...18 == 5 | lisa$Anzahl...21 == 2 | lisa$Anzahl...18 == "NA", ]  # filter for 1. 5 or 2
  lisa <- lisa[lisa$CADD_07.10.2016 > 20, ]  # filter for 2. CADD > 20
  lisa <- lisa[lisa$`Frequenz Exac` <  0.001 | lisa$`Frequenz Exac` == "NA", ]  # filter for 3. ExAC < 0.1% (< 0.001)
  # filter for Mutationstyp, Chromosom, Start, Stop, Gene Symbol, Ref and Alt
  lisa <- lisa[, c("Patienten-ID", "Mutationstyp", "Chromosom", "Start_neu", "Stop_neu", "Genname", "Ref_neu", "Allel_neu")]
  colnames(lisa) <- c("PatientID", "Type", "Chromosome", "Start", "End", "GeneSymbol", "Ref", "Alt")

  # 2. Load Fred's data
  fred <- as.data.frame(read_excel(paste0(inpath, "/2020-08-30_Frederick_Neukirch_Liste.xlsx"), sheet = "5 von 5 bzw. 2 von 2"))
  fred <- fred[!is.na(fred$Patient), ] # remove first row with NAs
  fred <- fred[fred$Anzahl...16 == "5" | fred$Anzahl...19 == 2 | is.na(fred$Anzahl...16), ]  # filter for 1. 5 or 2
  fred <- fred[!is.na(fred$`CADD_phred-Score`) & fred$`CADD_phred-Score` > 20, ]  # filter for 2. CADD > 20
  fred <- fred[is.numeric(fred$`Frequenz Exac`) <  0.001 | fred$`Frequenz Exac` == "NA", ]  # filter for 3. ExAC < 0.1% (< 0.001)
  # filter for Mutationstyp, Chromosom, Start, Stop, Gene Symbol, Ref and Alt
  fred <- fred[, c("Patient", "Mutationstyp", "Chromosom", "Start", "Stop", "Genname", "Ref", "Allel")]
  colnames(fred) <- c("PatientID", "Type", "Chromosome", "Start", "End", "GeneSymbol", "Ref", "Alt")

  # 3. Merge Lisa's and Fred's dataset
  index <- rbind(lisa, fred)
  index <- index[!duplicated(index), ]  # remove duplicated rows
  index <- separate_rows(index, "GeneSymbol", sep = ";")  # split multiple gene names
  rm(lisa, fred)
  
  ############################################
  ### Map Gene Symbols to Ensembl Gene IDs ###
  ############################################
  # We need to update some gene symbols as one is changed to a date in excel and two include different symbols
  index <- index[!(is.na(index$PatientID)), ]  # remove NA row
  index[index$GeneSymbol == "39692", "GeneSymbol"] <- "SEPT8"  # excel reads SEPT8 as the date 8th September, so we need to change it back to gene symbol SEPT8
  index[index$GeneSymbol == "DNAH7 (§x)", "GeneSymbol"] <- "DNAH7" # DNAH7 (§x) <- remove the (§x) behind the gene symbol
  index[index$GeneSymbol == "MTBP (2x)", "GeneSymbol"] <- "MTBP" # MTBP (2x) <- remove the (2x) behind the gene symbol
  
  # Then, we map the gene symbols to Ensembl Gene IDs: 
  gencode_all <- import.gff(gencode, format = "gtf")  # Load GENCODE .gtf file for gene symbol to Ensembl Gene ID mapping
  genes <- as.data.frame(gencode_all[gencode_all$type == "gene"])  # get only genes
  symbol_to_id <- genes[, c("gene_id", "gene_name")]  # map gene symbol to Ensembl gene ID
  symbol_to_id <- symbol_to_id[!duplicated(symbol_to_id), ]  # remove duplicated rows
  symbol_to_id$gene_id <- sapply(symbol_to_id$gene_id, function(x) strsplit(x, "\\.")[[1]][1])  # remove .xx endings
  rm(gencode_all, genes)
  
  # map eight "old" to current gene symbols
  index[index$GeneSymbol == "AC113191.1", "GeneSymbol"] <- "RP11-598P20.3"
  index[index$GeneSymbol == "AC226150.2", "GeneSymbol"] <- "ANHX"
  index[index$GeneSymbol == "GFY", "GeneSymbol"] <- "CTD-3148I10.1"
  index[index$GeneSymbol == "KIAA1797/FOCAD", "GeneSymbol"] <- "FOCAD"
  index[index$GeneSymbol == "RP11-123B3.6", "GeneSymbol"] <- "ERCC6-PGBD3"
  index[index$GeneSymbol == "RP11-240B13.2", "GeneSymbol"] <- "OC90"
  index[index$GeneSymbol == "RP11-332P22.1", "GeneSymbol"] <- "GMNC"
  index[index$GeneSymbol == "RP3-402G11.5", "GeneSymbol"] <- "SELO"
  
  # map Index patient data to Ensembl Gene IDs
  index <- merge(index, symbol_to_id, by.x = "GeneSymbol", by.y = "gene_name", all.x = T)
  index[is.na(index$gene_id), "GeneSymbol"] <- unlist(lapply(index[is.na(index$gene_id), "GeneSymbol"], function(x) ifelse(length(alias2Symbol(x)) != 0, alias2Symbol(x), x)))
  index <- merge(index, symbol_to_id, by.x = "GeneSymbol", by.y = "gene_name", all.x = T)
  colnames(index) <- c(colnames(index)[1:length(index)-1], "gene_id")
  index[is.na(index$gene_id), "GeneSymbol"]
  
  # for four gene symbols Gencode has no Ensembl ID given, so we map them manually
  index[index$GeneSymbol == "CFAP46", "gene_id"] <- "ENSG00000171811"
  index[index$GeneSymbol == "NPY4R2", "gene_id"] <- "ENSG00000264717"
  index[index$GeneSymbol == "TASOR", "gene_id"] <- "ENSG00000163946"
  index[index$GeneSymbol == "TEPSIN", "gene_id"] <- "ENSG00000167302"
  index$gene_id.x <- NULL

  ###############################
  ### Generate GRanges object ###
  ###############################
  # generate SNP IDs (chr:start)
  index$SNP <- paste0(index$Chromosome, ":", index$Start)
  colnames(index) <- c("GeneSymbol", "PatientID", "Type", "chr", "start", "end", "Ref", "Alt", "Ensembl Gene ID", "SNP")
  
  # generate GRanges object
  index_gr <- GenomicRanges::makeGRangesFromDataFrame(index, keep.extra.columns = TRUE, ignore.strand = TRUE)   # generate GRanges object only with chr, position, rs SNP ID and p-value
  GenomeInfoDb::seqlevelsStyle(index_gr) <- "NCBI"
  
  # save GRanges object to an rds file
  saveRDS(index_gr, file = paste0(outpath, "/Index-patients_rare_variants_filtered.rds"))
  
  ####################
  ### Session Info ###
  ####################
  print(sessionInfo())
}
