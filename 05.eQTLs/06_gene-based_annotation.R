# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(oligo)
})

### 1. LOAD eQTL results ###
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

### 2. Gene-based anotation ###
### 2.1 Transform eQTL data into gene/transcript/exon-based GRanges objects ###
# Gene eQTLs
eqtls.gene.gr <- GRanges(seqnames = eqtls.gene$Gene_location, strand = eqtls.gene$Gene_strand, gene_id = eqtls.gene$Gene, gene_name = eqtls.gene$Gene_symbol, p_value = eqtls.gene$`p-value`, fdr = eqtls.gene$FDR)
eqtls.gene.gr <- eqtls.gene.gr[!duplicated(eqtls.gene.gr)]  # 864 genes
# Gene cases eQTLs
eqtls.gene.cases.gr <- GRanges(seqnames = eqtls.gene.cases$Gene_location, strand = eqtls.gene.cases$Gene_strand, gene_id = eqtls.gene.cases$Gene, gene_name = eqtls.gene.cases$Gene_symbol, p_value = eqtls.gene.cases$`p-value`, fdr = eqtls.gene.cases$FDR)
eqtls.gene.cases.gr <- eqtls.gene.cases.gr[!duplicated(eqtls.gene.cases.gr)]  # 408 genes
# Transcript eQTLs
eqtls.transcript.gr <- GRanges(seqnames = eqtls.transcript$Transcript_location, strand = eqtls.transcript$Gene_strand, transcript_id = eqtls.transcript$Transcript, p_value = eqtls.transcript$pvalue, fdr = eqtls.transcript$FDR)
eqtls.transcript.gr <- eqtls.transcript.gr[!duplicated(eqtls.transcript.gr)]  # 1,993 transcripts
# Transcript cases eQTLs
eqtls.transcript.cases.gr <- GRanges(seqnames = eqtls.transcript.cases$Transcript_location, strand = eqtls.transcript.cases$Gene_strand, transcript_id = eqtls.transcript.cases$Transcript, p_value = eqtls.transcript.cases$pvalue, fdr = eqtls.transcript.cases$FDR)
eqtls.transcript.cases.gr <- eqtls.transcript.cases.gr[!duplicated(eqtls.transcript.cases.gr)]  # 725 transcripts
# Exon eQTLs
eqtls.exon.gr <- GRanges(seqnames = eqtls.exon$Exon_location, strand = eqtls.exon$Gene_strand, exon_id = eqtls.exon$Exon, p_value = eqtls.exon$pvalue, fdr = eqtls.exon$FDR)
eqtls.exon.gr <- eqtls.exon.gr[!duplicated(eqtls.exon.gr)]  # 6,001 exons
# Exon cases eQTLs
eqtls.exon.cases.gr <- GRanges(seqnames = eqtls.exon.cases$Exon_location, strand = eqtls.exon.cases$Gene_strand, exon_id = eqtls.exon.cases$Exon, p_value = eqtls.exon.cases$pvalue, fdr = eqtls.exon.cases$FDR)
eqtls.exon.cases.gr <- eqtls.exon.cases.gr[!duplicated(eqtls.exon.cases.gr)]  # 2,638 exons
# remove unnecessary data
rm(eqtls.gene, eqtls.gene.cases, eqtls.transcript, eqtls.transcript.cases, eqtls.exon, eqtls.exon.cases)

### 2.2 Public data ###
# load all data sets and generate GRanges based on gene, transcript (isoforms) or exon positions
### 2.2.1 BIOS whole blood gene eQTLs
bios.gene <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/gene_level_eQTLs.txt", sep = "\t")
bios.gene[, length(bios.gene)] <- NULL # remove last column
bios.gene <- bios.gene[bios.gene$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.gene.gr <- GRanges(seqnames = bios.gene$ProbeChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.gene$ProbeCenterChrPos)), 
                                         end = as.numeric(as.character(bios.gene$ProbeCenterChrPos))),
                        gene_id = bios.gene$ProbeName,
                        gene_name = bios.gene$HGNCName)  # generate GRanges object
bios.gene.gr <- bios.gene.gr[!duplicated(bios.gene.gr)]  # remove duplicated rows: 5,526,636 => 23039 genes
seqlevelsStyle(bios.gene.gr) <- "UCSC"
rm(bios.gene)

### 2.2.2 BIOS whole blood exon eQTLs
bios.exon <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/exon_level_eQTLs_independent_effects.txt", sep = "\t")
bios.exon <- bios.exon[bios.exon$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.exon.gr <- GRanges(seqnames = bios.exon$ExonChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.exon$ExonCenterChrPos)), 
                                         end = as.numeric(as.character(bios.exon$ExonCenterChrPos))),
                        exon_id = bios.exon$ExonId,
                        gene_name = bios.exon$GeneName)  # generate GRanges object
bios.exon.gr <- bios.exon.gr[!duplicated(bios.exon.gr)]  # remove duplicated rows: 171,889 => 105,097 exons
#mcols(bios.exon.gr) <- dplyr::select(bios.exon.gr, -c(ExonChr, ExonCenterChrPos, ExonId, GeneName))  # add all other columns to the dataset
seqlevelsStyle(bios.exon.gr) <- "UCSC"
rm(bios.exon)

### 2.2.3 CMC (Fromer et al.) gene eQTLs
cmc.gene <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc.gene <- cmc.gene[cmc.gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.gene.gr <- GRanges(seqnames = cmc.gene$Gene_location,
                       strand = cmc.gene$Gene_strand,
                       gene_id = cmc.gene$Gene,
                       gene_name = cmc.gene$Gene_symbol)  # generate GRanges object 
cmc.gene.gr <- cmc.gene.gr[!duplicated(cmc.gene.gr)]  # remove duplicated rows: 2,199,772 => 13,291 genes
rm(cmc.gene)

### 2.2.4 CMC (Fromer et al.) isoform eQTLs
cmc.iso <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt")
cmc.iso <- cmc.iso[cmc.iso$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.iso.gr <- GRanges(seqnames = cmc.iso$Transcript_location,
                      transcript_id = cmc.iso$Transcript,
                      gene_name = cmc.iso$Gene_symbol)  # generate GRanges object 
cmc.iso.gr <- cmc.iso.gr[!duplicated(cmc.iso.gr)]  # remove duplicated rows: 3,394,514 => 27,759 transcripts/isoforms
rm(cmc.iso)

### 2.2.5 GTEx Whole Blood
gtex.wb <- fread("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/Whole_Blood.v7.egenes.txt")
gtex.wb <- gtex.wb[gtex.wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
gtex.wb.gr <- GRanges(seqnames = gtex.wb$gene_chr, 
                      ranges = IRanges(start = as.numeric(as.character(gtex.wb$gene_start)), 
                                       end = as.numeric(as.character(gtex.wb$gene_end))),
                      strand = gtex.wb$strand,
                      gene_id = gtex.wb$gene_id,
                      gene_name = gtex.wb$gene_name)  # generate GRanges object only with gene ID, name and location
gtex.wb.gr <- gtex.wb.gr[!duplicated(gtex.wb.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
seqlevelsStyle(gtex.wb.gr) <- "UCSC"
rm(gtex.wb)

### 2.2.6 GTEx Brain (all brain tissues)
brains <- c("Brain_Amygdala.v7.egenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt", "Brain_Caudate_basal_ganglia.v7.egenes.txt", "Brain_Cerebellar_Hemisphere.v7.egenes.txt", "Brain_Cerebellum.v7.egenes.txt", "Brain_Cortex.v7.egenes.txt", "Brain_Frontal_Cortex_BA9.v7.egenes.txt", "Brain_Hippocampus.v7.egenes.txt", "Brain_Hypothalamus.v7.egenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v7.egenes.txt", "Brain_Putamen_basal_ganglia.v7.egenes.txt", "Brain_Spinal_cord_cervical_c-1.v7.egenes.txt", "Brain_Substantia_nigra.v7.egenes.txt")
gtex.brains <- data.frame()
for (i in brains) {
  filename <- paste0("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/", i)
  data <- fread(filename)
  data <- data[data$qval < 0.05, ]  # filter for significant hits
  gtex.brains <- rbind(gtex.brains, data) # merge different brain files together
}
gtex.brains.gr <- GRanges(seqnames = gtex.brains$gene_chr, 
                          ranges = IRanges(start = as.numeric(as.character(gtex.brains$gene_start)), 
                                           end = as.numeric(as.character(gtex.brains$gene_end))),
                          strand = gtex.brains$strand,  
                          gene_id = gtex.brains$gene_id,
                          gene_name = gtex.brains$gene_name)  # generate GRanges object only with gene ID, name and location
gtex.brains.gr <- gtex.brains.gr[!duplicated(gtex.brains.gr)]  # remove duplicated rows: 58,589 => 17,373 genes
seqlevelsStyle(gtex.brains.gr) <- "UCSC"
rm(brains, filename, data, i, gtex.brains)

### 2.3 Load background data ###
### 2.3.1 Load all genes (for which we calculated eQTLs)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
fData(geneSummaries) <- separate_rows(fData(geneSummaries), c(gene_id, location, gene_name), sep = "\\|")
genes <- GRanges(seqnames = unlist(lapply(fData(geneSummaries)$location, function(x) strsplit(x, ":")[[1]][1])),
                 ranges = unlist(lapply(fData(geneSummaries)$location, function(x) strsplit(x, ":")[[1]][2])),
                 gene_id = fData(geneSummaries)$gene_id,
                 gene_name = fData(geneSummaries)$gene_name)
genes <- genes[!duplicated(genes)]  # 17,447 genes
rm(geneSummaries)

# 2.3.2 Load all transcripts (for which we calculated eQTLs)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
fData(transcriptSummaries) <- separate_rows(fData(transcriptSummaries), c(transcript_id, location, gene_id, gene_name), sep = "\\|")
transcripts <- GRanges(seqnames = unlist(lapply(fData(transcriptSummaries)$location, function(x) strsplit(x, ":")[[1]][1])),
                 ranges = unlist(lapply(fData(transcriptSummaries)$location, function(x) strsplit(x, ":")[[1]][2])),
                 transcript_id = fData(transcriptSummaries)$transcript_id,
                 gene_name = fData(transcriptSummaries)$gene_name)
transcripts <- transcripts[!duplicated(transcripts)]  # 100,750 transcripts
rm(transcriptSummaries)

# 2.3.3 Load all exons (for which we calculated eQTLs)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
fData(exonSummaries) <- separate_rows(fData(exonSummaries), c(exon_id, location, gene_id, gene_name), sep = "\\|")
exons <- GRanges(seqnames = unlist(lapply(fData(exonSummaries)$location, function(x) strsplit(x, ":")[[1]][1])),
                 ranges = unlist(lapply(fData(exonSummaries)$location, function(x) strsplit(x, ":")[[1]][2])),
                 exon_id = fData(exonSummaries)$exon_id,
                 gene_name = fData(exonSummaries)$gene_name)
exons <- exons[!duplicated(exons)]  # 242,443 exons
rm(exonSummaries)

### 2.4 Annotation ###
# function, to calculate the gene annotation to a public dataset
# own: own cis-eQTL data set as GRanges object
# public: a public cis-eQTL data set as GRanges object
# nperm: number of permutations (usually a 1000 times)
# background: all background genes coming from the gene expression (here: 14,048 unqiue genes)
annotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  # resample nperm times
  resampling <- lapply(1:nperm, function(x) sum(overlapsAny(sample(background, length(own)), public)))
  # IF 0 => echt kleiner als 1e^(-5/6)
  p.value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p.value), as.character(or))) # return the resampling p-value
}

results <- data.frame()
perm = 10000

### 2.4.1 annotation to BIOS whole blood gene cis-eQTLs
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, bios.gene.gr, perm, genes)
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, bios.gene.gr, perm, genes)
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, bios.gene.gr, perm, transcripts)
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, bios.gene.gr, perm, transcripts)
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, bios.gene.gr, perm, exons)
results[c("BIOS_whole_blood_gene_cis-eQTL_overlap", "BIOS_whole_blood_gene_cis-eQTL_pVal", "BIOS_whole_blood_gene_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, bios.gene.gr, perm, exons)

### 2.4.2 annotation to BIOS whole blood exon cis-eQTLs
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, bios.exon.gr, perm, genes)
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, bios.exon.gr, perm, genes)
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, bios.exon.gr, perm, transcripts)
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, bios.exon.gr, perm, transcripts)
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, bios.exon.gr, perm, exons)
results[c("BIOS_whole_blood_exon_cis-eQTL_overlap", "BIOS_whole_blood_exon_cis-eQTL_pVal", "BIOS_whole_blood_exon_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, bios.exon.gr, perm, exons)

### 2.4.3 annotation to CMC (Fromer et al.) DPC gene cis-eQTLs
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, cmc.gene.gr, perm, genes)
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, cmc.gene.gr, perm, genes)
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, cmc.gene.gr, perm, transcripts)
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, cmc.gene.gr, perm, transcripts)
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, cmc.gene.gr, perm, exons)
results[c("CMC_gene_cis-eQTL_overlap", "CMC_gene_cis-eQTL_pVal", "CMC_gene_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, cmc.gene.gr, perm, exons)

### 2.4.4 annotation to CMC (Fromer et al.) DPC isoform cis-eQTLs
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, cmc.iso.gr, perm, genes)
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, cmc.iso.gr, perm, genes)
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, cmc.iso.gr, perm, transcripts)
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, cmc.iso.gr, perm, transcripts)
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, cmc.iso.gr, perm, exons)
results[c("CMC_isoform_cis-eQTL_overlap", "CMC_isoform_cis-eQTL_pVal", "CMC_isoform_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, cmc.iso.gr, perm, exons)

### 2.4.5 annotation to GTEx whole blood gene cis-eQTLs
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, gtex.wb.gr, perm, genes)
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, gtex.wb.gr, perm, genes)
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, gtex.wb.gr, perm, transcripts)
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, gtex.wb.gr, perm, transcripts)
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, gtex.wb.gr, perm, exons)
results[c("GTEx_whole_blood_gene_cis-eQTL_overlap", "GTEx_whole_blood_gene_cis-eQTL_pVal", "GTEx_whole_blood_gene_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, gtex.wb.gr, perm, exons)

### 2.4.6 annotation to GTEx brain gene cis-eQTLs
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, gtex.brains.gr, perm, genes)
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, gtex.brains.gr, perm, genes)
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, gtex.brains.gr, perm, transcripts)
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, gtex.brains.gr, perm, transcripts)
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, gtex.brains.gr, perm, exons)
results[c("GTEx_brain_gene_cis-eQTL_overlap", "GTEx_brain_gene_cis-eQTL_pVal", "GTEx_brain_gene_cis-eQTL_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, gtex.brains.gr, perm, exons)

### 2.4.6 write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
