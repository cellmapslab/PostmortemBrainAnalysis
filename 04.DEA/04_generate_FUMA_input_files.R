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
  
# 2. Load DEA and eQTS results
gene <- read.delim(paste0(input_path, "/1.1_diff_expr_genes_for_diagnosis.csv"))
print(paste0("Total gene counts: ", dim(gene)[1]))  # 17,478
transcript <- read.delim(paste0(input_path, "/2.4_diff_expr_transcripts_for_diagnosis_new_pvalue_and_fdr.csv"))
print(paste0("Total transcript counts: ", dim(transcript)[1]))  # 17,507  
exon <- read.delim(paste0(input_path, "/3.4_diff_expr_exons_for_diagnosis_new_pvalue_and_fdr.csv"))
print(paste0("Total exon counts: ", dim(exon)[1]))  # 17,522
asd_teqts <- read.delim(paste0(input_path, "/2.5.2_diff_expr_transcripts_for_ASD-2019_PRS_new_pvalue_and_fdr.csv"))
mdd_teqts <- read.delim(paste0(input_path, "/2.5.13_diff_expr_transcripts_for_MDD3-2019_PRS_new_pvalue_and_fdr.csv"))
scz_teqts <- read.delim(paste0(input_path, "/2.5.14_diff_expr_transcripts_for_SCZ2-2014_PRS_new_pvalue_and_fdr.csv"))
adhd_eeqts <- read.delim(paste0(input_path, "/3.5.1_diff_expr_exons_for_ADHD-2019_PRS_new_pvalue_and_fdr.csv"))
asd_eeqts <- read.delim(paste0(input_path, "/3.5.2_diff_expr_exons_for_ASD-2019_PRS_new_pvalue_and_fdr.csv"))
bip_eeqts <- read.delim(paste0(input_path, "/3.5.3_diff_expr_exons_for_BIP-2019_PRS_new_pvalue_and_fdr.csv"))
cdg_eeqts <- read.delim(paste0(input_path, "/3.5.7_diff_expr_exons_for_CDG2-2019_PRS_new_pvalue_and_fdr.csv"))
mdd_eeqts <- read.delim(paste0(input_path, "/3.5.13_diff_expr_exons_for_MDD3-2019_PRS_new_pvalue_and_fdr.csv"))
scz_eeqts <- read.delim(paste0(input_path, "/3.5.14_diff_expr_exons_for_SCZ2-2014_PRS_new_pvalue_and_fdr.csv"))

# 3. Get only significant hits (as unique Ensembl Gene IDs)
gene_sig <- gene[gene$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("gDEGenes: ", length(unique(gene_sig))))  # 0
transcript_sig  <- transcript[transcript$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("tDEGenes: ", length(unique(transcript_sig))))  # 6  
exon_sig <- exon[exon$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("eDEGenes: ", length(unique(exon_sig))))  # 2,223
print(paste0("Overlapping tDEGenes and eDEGenes: ", length(intersect(transcript_sig, exon_sig))))
asd_teqts <- asd_teqts[asd_teqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("ASD-2019 transcript eQTS: ", length(unique(asd_teqts))))  # 2
mdd_teqts <- mdd_teqts[mdd_teqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("MDD-2019 transcript eQTS: ", length(unique(mdd_teqts))))  # 23
scz_teqts <- scz_teqts[scz_teqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("SCZ-2014 transcript eQTS: ", length(unique(scz_teqts))))  # 4
adhd_eeqts <- adhd_eeqts[adhd_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("ADHD-2019 exon eQTS: ", length(unique(adhd_eeqts))))  # 3
asd_eeqts <- asd_eeqts[asd_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("ASD-2019 exon eQTS: ", length(unique(asd_eeqts))))  # 4693
bip_eeqts <- bip_eeqts[bip_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("BIP-2019-2014 exon eQTS: ", length(unique(bip_eeqts))))  # 5
cdg_eeqts <- cdg_eeqts[cdg_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("CDG-2019 exon eQTS: ", length(unique(cdg_eeqts))))  # 5000
mdd_eeqts <- mdd_eeqts[mdd_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("MDD-2019 exon eQTS: ", length(unique(mdd_eeqts))))  # 8359
scz_eeqts <- scz_eeqts[scz_eeqts$FDR.estimate <= 0.1, "Ensembl.Gene.ID"]
print(paste0("SCZ-2014 exon eQTS: ", length(unique(scz_eeqts))))  # 3582

# combine the five exon eQTS of ADHD-2019, ASD-2019, BIP-2019, MDD-2019 and SCZ-2014
length(c(adhd_eeqts, asd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts))
length(unique(c(adhd_eeqts, asd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts)))
aabms_eeqts <- unique(c(adhd_eeqts, asd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts))

# combine the four exon eQTS of ADHD-2019, BIP-2019, MDD-2019 and SCZ-2014
length(c(adhd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts))
length(unique(c(adhd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts)))
abms_eeqts <- unique(c(adhd_eeqts, bip_eeqts, mdd_eeqts, scz_eeqts))

# combine the three exon eQTS of BIP-2019, MDD-2019 and SCZ-2014
length(c(bip_eeqts, mdd_eeqts, scz_eeqts))
length(unique(c(bip_eeqts, mdd_eeqts, scz_eeqts)))
bms_eeqts <- unique(c(bip_eeqts, mdd_eeqts, scz_eeqts))

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
asd_teqts <- sapply(strsplit(as.character(asd_teqts), "\\."), "[", 1)
mdd_teqts <- sapply(strsplit(as.character(mdd_teqts), "\\."), "[", 1)
scz_teqts <- sapply(strsplit(as.character(scz_teqts), "\\."), "[", 1)
adhd_eeqts <- sapply(strsplit(as.character(adhd_eeqts), "\\."), "[", 1)
asd_eeqts <- sapply(strsplit(as.character(asd_eeqts), "\\."), "[", 1)
bip_eeqts <- sapply(strsplit(as.character(bip_eeqts), "\\."), "[", 1)
cdg_eeqts <- sapply(strsplit(as.character(cdg_eeqts), "\\."), "[", 1)
mdd_eeqts <- sapply(strsplit(as.character(mdd_eeqts), "\\."), "[", 1)
scz_eeqts <- sapply(strsplit(as.character(scz_eeqts), "\\."), "[", 1)
aabms_eeqts <- sapply(strsplit(as.character(aabms_eeqts), "\\."), "[", 1)
abms_eeqts <- sapply(strsplit(as.character(abms_eeqts), "\\."), "[", 1)
bms_eeqts <- sapply(strsplit(as.character(bms_eeqts), "\\."), "[", 1)

# 4. Generate FUMA input files for background genes (based on Ensembl Gene IDs)
write.table(unique(gene_bkgd), file = paste0(input_files_path, "/1_gene_bkgd.txt"), row.names = F,col.names = F, quote = F)
write.table(unique(transcript_bkgd), file = paste0(input_files_path, "/2_transcript_bkgd.txt"), row.names = F,col.names = F, quote = F)
write.table(unique(exon_bkgd), file = paste0(input_files_path, "/3_exon_bkgd.txt"), row.names = F,col.names = F, quote = F)

# 5. Generate FUMA input files for t/eDEGenes (based on Ensembl Gene IDs)
write.table(unique(transcript_sig), file = paste0(input_files_path, "/2_tDEGenes.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(asd_teqts), file = paste0(input_files_path, "/2_transcript_eQTS_ASD-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(mdd_teqts), file = paste0(input_files_path, "/2_transcript_eQTS_MDD-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(scz_teqts), file = paste0(input_files_path, "/2_transcript_eQTS_SCZ-2014.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(adhd_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_ADHD-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(asd_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_ASD-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(bip_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_BIP-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(cdg_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_CDG-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(mdd_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_MDD-2019.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(scz_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_SCZ-2014.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(aabms_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_ADHD_ASD_BIP_MDD_SCZ.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(abms_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_ADHD_BIP_MDD_SCZ.txt"), row.names = F, col.names = F, quote = F)
write.table(unique(bms_eeqts), file = paste0(input_files_path, "/3_exon_eQTS_BIP_MDD_SCZ.txt"), row.names = F, col.names = F, quote = F)

# 6. Generate FUMA input file for KEGG pathways (just load and save as needed)
load(kegg_file)
toGmt(kegg, paste0(input_files_path, "/4_kegg_entrez_no_disease_pathways.txt"))

print(sessionInfo())
