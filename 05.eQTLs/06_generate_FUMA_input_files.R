###################################
### FUMA - Generate input files ###
###################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log
# or with Rscript: Rscript --slave 06_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)
library(oligo)
library(GSEABase)

# First read the arguments listed in the commandline, trailingOnly uses only the parameters after --args
args = (commandArgs(trailingOnly = TRUE))	

# go through the list of arguments
if(length(args) == 0){
  cat("\nScript to generate input files for matrixEQTL with samples keeping in the same order.\n\nYou can invoke the script as follows:\nRscript 02_generate_input_files.R bim=\"/home/.../genotype-file.bim\" snps=\"/home/.../genotype-file.gen\" gene=\"/home/.../gene-expr.RData\" transcript=\"/home/.../transcript-expr.RData\" exon=\"/home/.../exon-expr.RData\" outpath=\"/home/to/store/output/\"\n\nCommandline parameters:\n\nbim\t\tpath and filename of the duplicated removed and merged genotype file with .bim ending\nsnps\t\tpath and filename of the duplicated removed, merged and probabilities changed genotype file with .gen (or .gen.gz) ending\ngene\t\tpath and filename of the .RData file including the preprocessed gene expression data\ntranscript\tpath and filename of the .RData file including the preprocessed transcript expression data\nexon\t\tpath and filename of the .RData file including the preprocessed exon expression data\noutpath\t\tpath to store the output into\n\n")
} else {
  for(i in 1:length(args)){
    eval(parse(text = args[[i]]))
  }
  # 1. Save input path and generate output paths if not existing already
  input_files_path <- paste0(outpath, "/7.1_input_files")
  geGenes_path <- paste0(outpath, "/7.2_gene_eQTLs")
  geGenes_cases_path <- paste0(outpath, "/7.2_gene_eQTLs_cases")
  teGenes_path <- paste0(outpath, "/7.3_transcript_eQTLs")
  teGenes_cases_path <- paste0(outpath, "/7.3_transcript_eQTLs_cases")
  eeGenes_path <- paste0(outpath, "/7.4_exon_eQTLs")
  eeGenes_cases_path <- paste0(outpath, "/7.4_exon_eQTLs_cases")
  system(paste0("mkdir -p ", input_files_path, " ", geGenes_path, " ", geGenes_cases_path, " ", teGenes_path, " ", teGenes_cases_path, " ", eeGenes_path, " ", eeGenes_cases_path))
  
  # 2. Load eQTL result tables
  geGenes <- fread(paste0(inpath, "/01_pmb_cis-eQTL_gene_hits_cases_vs_controls.txt"), select = "Ensembl Gene ID", data.table = F)
  geGenes <- unique(unlist(lapply(geGenes$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(geGenes))))  # 1,115
  geGenes_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_gene_hits_cases_only.txt"), select = "Ensembl Gene ID", data.table = F)
  geGenes_cases <- unique(unlist(lapply(geGenes_cases$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(geGenes_cases))))  # 433
  teGenes <- fread(paste0(inpath, "/01_pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt"), select = "Ensembl Gene ID", data.table = F)
  teGenes <- unique(unlist(lapply(teGenes$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(teGenes))))  # 3,159
  teGenes_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_transcript_hits_cases_only.txt"), select = "Ensembl Gene ID", data.table = F)
  teGenes_cases <- unique(unlist(lapply(teGenes_cases$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(teGenes_cases))))  # 1,364
  eeGenes <- fread(paste0(inpath, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"), select = "Ensembl Gene ID", data.table = F)
  eeGenes <- unique(unlist(lapply(eeGenes$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(eeGenes))))  # 6,401
  eeGenes_cases <- fread(paste0(inpath, "/02_pmb_cis-eQTL_exon_hits_cases_only.txt"), select = "Ensembl Gene ID", data.table = F)
  eeGenes_cases <- unique(unlist(lapply(eeGenes_cases$`Ensembl Gene ID`, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Total gene counts: ", length(unique(eeGenes_cases))))  # 3,674

  # 3. Here, we already have significant eQTL hits. So we only need to get the background gene sets from the expression datasets
  load(paste0(dea, "/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData"))
  gene_bkgd <- unique(unlist(lapply(fData(gene_summaries)$gene_id, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Background genes (gene-level): ", length(gene_bkgd)))  # 17,478
  load(paste0(dea, "/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData"))
  transcript_bkgd <- unique(unlist(lapply(fData(transcript_summaries)$gene_id, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Background genes (transcript-level): ", length(transcript_bkgd)))  # 17,485
  load(paste0(dea, "/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData"))
  exon_bkgd <- unique(unlist(lapply(fData(exon_summaries)$gene_id, function(x) strsplit(x, "\\|"))))  # get all gene IDs and split them if they were merged with a "|"
  print(paste0("Background genes (exon-level): ", length(exon_bkgd)))  # 17,496

  # 5. Remove version numbers from Ensembl Gene IDs, e.g. ENSG00000000419.8 to ENSG00000000419
  geGenes <- sapply(strsplit(as.character(geGenes), "\\."), "[", 1)
  geGenes_cases <- sapply(strsplit(as.character(geGenes_cases), "\\."), "[", 1)
  teGenes <- sapply(strsplit(as.character(teGenes), "\\."), "[", 1)
  teGenes_cases <- sapply(strsplit(as.character(teGenes_cases), "\\."), "[", 1)
  eeGenes <- sapply(strsplit(as.character(eeGenes), "\\."), "[", 1)
  eeGenes_cases <- sapply(strsplit(as.character(eeGenes_cases), "\\."), "[", 1)
  gene_bkgd <- sapply(strsplit(as.character(gene_bkgd), "\\."), "[", 1)
  transcript_bkgd <- sapply(strsplit(as.character(transcript_bkgd), "\\."), "[", 1)
  exon_bkgd <- sapply(strsplit(as.character(exon_bkgd), "\\."), "[", 1)

  # 4. Generate FUMA input files for background genes (based on Ensembl Gene IDs)
  write.table(unique(gene_bkgd), file = paste0(input_files_path, "/1_gene_bkgd.txt"), row.names = F,col.names = F, quote = F)
  write.table(unique(transcript_bkgd), file = paste0(input_files_path, "/2_transcript_bkgd.txt"), row.names = F,col.names = F, quote = F)
  write.table(unique(exon_bkgd), file = paste0(input_files_path, "/3_exon_bkgd.txt"), row.names = F,col.names = F, quote = F)

  # 5. Generate FUMA input files for eGenes/eQTL genes (based on Ensembl Gene IDs)
  write.table(unique(geGenes), file = paste0(input_files_path, "/1_geGenes.txt"), row.names = F, col.names = F, quote = F)
  write.table(unique(geGenes_cases), file = paste0(input_files_path, "/1_geGenes_cases.txt"), row.names = F, col.names = F, quote = F)
  write.table(unique(teGenes), file = paste0(input_files_path, "/2_teGenes.txt"), row.names = F, col.names = F, quote = F)
  write.table(unique(teGenes_cases), file = paste0(input_files_path, "/2_teGenes_cases.txt"), row.names = F, col.names = F, quote = F)
  write.table(unique(eeGenes), file = paste0(input_files_path, "/3_eeGenes.txt"), row.names = F, col.names = F, quote = F)
  write.table(unique(eeGenes_cases), file = paste0(input_files_path, "/3_eeGenes_cases.txt"), row.names = F, col.names = F, quote = F)

  # 6. Generate FUMA input file for KEGG pathways (just load and save as needed)
  load(keggf)
  toGmt(kegg, paste0(input_files_path, "/4_kegg_entrez_no_disease_pathways.txt"))
  
  print(sessionInfo())
}