###################################
### FUMA - Generate input files ###
###################################
# Invoke the script in the commandline as follows: R CMD BATCH --slave 05_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log
# or with Rscript: Rscript --slave 05_generate_FUMA_input_files.R &>> generate_FUMA_input_files.log

# add installed packages from home path to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(GSEABase)

# 1. Save input path and generate output paths if not existing already
kegg_file <- "/storage/groups/ccm01/datasets/2020_KEGG/2020-06-14_KEGG_entrez_no_disease_pathways.rds"
input_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results"
output_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/04_FUMA"
input_files_path <- paste0(output_path, "/4.1_input_files")
tdegenes_path <- paste0(output_path, "/4.2_tDEGenes_result")
edegenes_path <- paste0(output_path, "/4.3_eDEGenes_result")
system(paste0("mkdir -p ", output_path, " ", input_files_path, " ", tdegenes_path, " ", edegenes_path))
  
# 2. Load DEA results
gene <- read.delim(paste0(input_path, "/1.1_diff_expr_genes_for_diagnosis.csv"))
print(paste0("Total gene counts: ", dim(gene)[1]))  # 17,478
transcript <- read.delim(paste0(input_path, "/2.4_diff_expr_transcripts_for_diagnosis_new_pvalue_and_fdr.csv"))
print(paste0("Total transcript counts: ", dim(transcript)[1]))  # 17,507  
exon <- read.delim(paste0(input_path, "/3.4_diff_expr_exons_for_diagnosis_new_pvalue_and_fdr.csv"))
print(paste0("Total exon counts: ", dim(exon)[1]))  # 17,522

# 3. Get only significant hits (as unique Ensembl Gene IDs)
gene_sig <- gene[gene$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("gDEGenes: ", length(unique(gene_sig))))  # 0
transcript_sig  <- transcript[transcript$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("tDEGenes: ", length(unique(transcript_sig))))  # 6  
exon_sig <- exon[exon$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("eDEGenes: ", length(unique(exon_sig))))  # 2,223
print(paste0("Overlapping tDEGenes and eDEGenes: ", length(intersect(transcript_sig, exon_sig))))

# 4. Get background genes (as unique Ensembl Gene IDs)
gene_bkgd <- unique(gene$Ensembl.Gene.ID)
print(paste0("Background genes (gene-level): ", length(gene_bkgd)))  # 17,478
transcript_bkgd <- unique(transcript$Ensembl.Gene.ID)
print(paste0("Background genes (transcript-level): ", length(transcript_bkgd)))  # 17,485
exon_bkgd <- unique(exon$Ensembl.Gene.ID)
print(paste0("Background genes (exon-level): ", length(exon_bkgd)))  # 17,496

# 5. Remove version numbers from Ensembl Gene IDs, e.g. ENSG00000000419.8 to ENSG00000000419
transcript_sig <- sapply(strsplit(as.character(transcript_sig), "\\."), "[", 1)
exon_sig <- sapply(strsplit(as.character(exon_sig), "\\."), "[", 1)
gene_bkgd <- sapply(strsplit(as.character(gene_bkgd), "\\."), "[", 1)
transcript_bkgd <- sapply(strsplit(as.character(transcript_bkgd), "\\."), "[", 1)
exon_bkgd <- sapply(strsplit(as.character(exon_bkgd), "\\."), "[", 1)

# 4. Generate FUMA input files for background genes (based on Ensembl Gene IDs)
write.table(unique(gene_bkgd), file = paste0(input_files_path, "/1_gene_bkgd.txt"), row.names = F,col.names = F, quote = F)
write.table(unique(transcript_bkgd), file = paste0(input_files_path, "/2_transcript_bkgd.txt"), row.names = F,col.names = F, quote = F)
write.table(unique(exon_bkgd), file = paste0(input_files_path, "/3_exon_bkgd.txt"), row.names = F,col.names = F, quote = F)

# 5. Generate FUMA input files for t/eDEGenes (based on Ensembl Gene IDs)
write.table(unique(transcript_sig), file = paste0(input_files_path, "/2_tDEGenes.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(exon_sig), file = paste0(input_files_path, "/3_eDEGenes.txt"), row.names = F, col.names = F, quote = F)

# 6. Generate FUMA input file for KEGG pathways (just load and save as needed)
load(kegg_file)
toGmt(kegg, paste0(input_files_path, "/4_kegg_entrez_no_disease_pathways.txt"))

print(sessionInfo())
