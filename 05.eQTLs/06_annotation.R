# Load locally installed packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

# load packages
library(viridis)
library(UpSetR)
library(ggplot2)
library(data.table)
library(rtracklayer)
library(reshape2)
library(VennDiagram)
library(oligo)
library(dplyr)

setwd("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/figures/")

### 1. LOAD own eQTLs ###
# cases and controls cis-eQTLs genes
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_SCZ_CTL_matrixEQTL_cis_results_genes.txt", header = T, stringsAsFactors = F)
cis.genes.all <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with FDR < 0.05
# cases and controls cis-eQTLs exons
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/01_pmb_SCZ_CTL_matrixEQTL_cis_results_exons.txt", header = T, stringsAsFactors = F)
cis.exons.all <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with FDR < 0.05
# only cases cis-eQTLs genes
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/02_pmb_SCZ_matrixEQTL_cis_results_genes.txt", header = T, stringsAsFactors = F)
cis.genes.cases <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with FDR < 0.05
# only cases cis-eQTLs exons
data <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/02_pmb_SCZ_matrixEQTL_cis_results_exons.txt", header = T, stringsAsFactors = F)
cis.exons.cases <- data[data$FDR < 0.05,  ]  # filter for only eQTLs with FDR < 0.05
rm(data)

### 2. Generate Output Files for significant eQTLs ###
# Load snp positions
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", header = F, stringsAsFactors = FALSE) # snp chr pos
colnames(snpspos) <- c("chromosome", "SNP", "distance", "position", "SNP_allele1", "SNP_allele2")
snpspos$SNP_location <- paste(snpspos$chromosome, ":", snpspos$position, sep = "")
snpspos <- snpspos[, c("SNP_location", "SNP", "SNP_allele1", "SNP_allele2")]

# Load gene postions
genes <- import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf")
genes <- genes[genes$gene_type == "protein_coding"]  # get only protein-coding genes
genes <- genes[genes$type == "gene"]  # get only exons from the genes
genes <- as.data.frame(genes)
genes <- genes[, c("gene_id", "gene_name", "seqnames", "start", "end", "strand")]
genes$Gene_location <- paste(genes$seqnames, ":", genes$start, "..", genes$end, sep = "")
genespos <- genes[, c("gene_id", "gene_name", "Gene_location", "strand")]
colnames(genespos) <- c("Gene", "Gene_symbol", "Gene_location", "Gene_strand")

# add info to gene eqtl data (cases vs. controls) and save it into a file
eqtls.gene <- merge(cis.genes.all, snpspos, by = "SNP", all.x = T)
eqtls.gene <- merge(eqtls.gene, genespos, by.x = "gene", by.y = "Gene_symbol", all.x = T)
eqtls.gene <- eqtls.gene[complete.cases(eqtls.gene), ]
eqtls.gene <- eqtls.gene[, c("SNP", "Gene", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
eqtls.gene$eQTL_type <- "cis" 
colnames(eqtls.gene) <- c("SNP", "Gene", "Gene_symbol", "Beta", "Test_statistic", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location",  "Gene_strand", "eQTL_type")
eqtls.gene <- eqtls.gene[order(eqtls.gene$FDR),]
write.table(eqtls.gene, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", row.names = F, quote = F, dec = ".", sep = " ")

# add info to gene eqtl data (only cases) and save it into a file
eqtls.gene.cases <- merge(cis.genes.cases, snpspos, by = "SNP", all.x = T)
eqtls.gene.cases <- merge(eqtls.gene.cases, genespos, by.x = "gene", by.y = "Gene_symbol", all.x = T)
eqtls.gene.cases <- eqtls.gene.cases[complete.cases(eqtls.gene.cases), ]
eqtls.gene.cases <- eqtls.gene.cases[, c("SNP", "Gene", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location", "Gene_strand")]
eqtls.gene.cases$eQTL_type <- "cis" 
colnames(eqtls.gene.cases) <- c("SNP", "Gene", "Gene_symbol", "Beta", "Test_statistic", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Gene_location",  "Gene_strand", "eQTL_type")
eqtls.gene.cases <- eqtls.gene.cases[order(eqtls.gene.cases$FDR),]
write.table(eqtls.gene.cases, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_only_cases.txt", row.names = F, quote = F, dec = ".", sep = " ")

# Load exon positions
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_exon_summaries_annotated_combat_batch_corrected.RData")
exonspos <- exonSummaries@featureData@data
exonspos$location <- gsub("-", "..", exonspos$location)

# add info to exon eqtl data (cases vs. controls) and save it into a file
eqtls.exon <- merge(cis.exons.all, snpspos, by = "SNP", all.x = T)
eqtls.exon <- merge(eqtls.exon, exonspos, by.x = "gene", by.y = "exon_id", all.x = T)
eqtls.exon <- merge(eqtls.exon, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.exon <- eqtls.exon[complete.cases(eqtls.exon), ]
eqtls.exon <- eqtls.exon[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id", "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.exon) <- c("SNP", "Exon", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.exon$eQTL_type <- "cis" 
eqtls.exon <- eqtls.exon[order(eqtls.exon$FDR),]
write.table(eqtls.exon, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", row.names = F, quote = F, dec = ".", sep = " ")

# add info to exon eqtl data (only cases) and save it into a file
eqtls.exon.cases <- merge(cis.exons.cases, snpspos, by = "SNP", all.x = T)
eqtls.exon.cases <- merge(eqtls.exon.cases, exonspos, by.x = "gene", by.y = "exon_id", all.x = T)
eqtls.exon.cases <- merge(eqtls.exon.cases, genespos, by.x = "gene_id", by.y = "Gene", all.x = T)
eqtls.exon.cases <- eqtls.exon.cases[complete.cases(eqtls.exon.cases), ]
eqtls.exon.cases <- eqtls.exon.cases[, c("SNP", "gene", "beta", "t-stat", "p-value", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "location", "gene_id", "gene_name", "Gene_location", "Gene_strand")]
colnames(eqtls.exon.cases) <- c("SNP", "Exon", "Beta", "Test_statistic", "pvalue", "FDR", "SNP_location", "SNP_allele1", "SNP_allele2", "Exon_location", "Gene", "Gene_symbol", "Gene_location", "Gene_strand")
eqtls.exon.cases$eQTL_type <- "cis" 
eqtls.exon.cases <- eqtls.exon.cases[order(eqtls.exon.cases$FDR),]
write.table(eqtls.exon.cases, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_only_cases.txt", row.names = F, quote = F, dec = ".", sep = " ")

# get Exon ID to Gene ID mapping
#load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_exon_summaries_annotated_combat_batch_corrected.RData")
#fData <- fData(exonSummaries)[, c("exon_id", "gene_id")]
#fData <- fData[!duplicated(fData), ]
#eqtls.exon <- merge(eqtls.exon, fData, by.x = "gene", by.y = "exon_id", all.x = T)

### 3. Upset Plots and VENN diagrams ###
### 3.1 SNP level
### 3.1.1 Upset plot
eqtls.snps <- list(genes_all = eqtls.gene$SNP, genes_cases_only = eqtls.gene.cases$SNP, exons_all = eqtls.exon$SNP, genes_cases_only = eqtls.exon.cases$SNP)
pdf("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/figures/01_upset_plot_cis_eqtls_snps.pdf", width = 14, height = 7, onefile = F)
upset(fromList(eqtls.snps), mainbar.y.label = "# SNPs", text.scale = 2)
dev.off()

### 3.1.2 VENN diagram
pdf("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/figures/01_venn_cis_eqtls_snps.pdf", width = 12, height = 11.5, onefile = F)
grid.newpage()
draw.pairwise.venn(area1 = length(unique(eqtls.exon$SNP)), area2 = length(unique(eqtls.gene$SNP)), 
                   cross.area = length(intersect(unique(eqtls.exon$SNP), unique(eqtls.gene$SNP))), 
                   category = c("Exon eQTLs", "Gene eQTLs"), lty = rep("blank", 2), fill = viridis(2),
                   alpha = rep(0.5, 2), cat.pos = c(-4, 6), cat.dist = rep(0.03, 1), cex = 4, cat.cex = 4) 
dev.off()

### 3.2 Gene/Exon level
### 3.2.1 Upset plot
eqtls.genes <- list(genes_all = eqtls.gene$Gene_symbol, genes_cases_only = eqtls.gene.cases$Gene_symbol, exons_all = eqtls.exon$Gene_symbol, genes_cases_only = eqtls.exon.cases$Gene_symbol)
pdf("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/figures/02_upset_plot_cis_eqtls_genes.pdf", width = 14, height = 7, onefile = F)
upset(fromList(eqtls.genes), mainbar.y.label = "# Genes", text.scale = 2)
dev.off()

### 3.2.2 VENN diagram
pdf("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/figures/02_venn_cis_eqtls_genes.pdf", width = 12, height = 11.5, onefile = F)
grid.newpage()
draw.pairwise.venn(area1 = length(unique(eqtls.exon$Gene_symbol)), area2 = length(unique(eqtls.gene$Gene_symbol)), 
                   cross.area = length(intersect(unique(eqtls.exon$Gene_symbol), unique(eqtls.gene$Gene_symbol))), 
                   category = c("Exon eQTLs", "Gene eQTLs"), lty = rep("blank", 2), fill = viridis(2),
                   alpha = rep(0.5, 2), cat.pos = c(-4, 6), cat.dist = rep(0.03, 1), cex = 4, cat.cex = 4) 
dev.off()


### 4. Boxplots top hits ###
### 4.1 Find top exon hit without gene hit => nr. 1 top hit with SNP rs662145 and exon ENSE00001139439.2 of gene USP24
# unoverlapping genes
no.overlap = setdiff(eqtls.exon$Gene_symbol, eqtls.gene$Gene_symbol)
length(no.overlap) # 5779
# overlapping genes
overlap = intersect(eqtls.exon$Gene_symbol, eqtls.gene$Gene_symbol)
length(overlap) # 65.9% (n = 663) of 1006 unique eQTL genes have an exon hit, while 5779 eQTL genes are only found with exons
# plot 2 Bsp overlap vs. non
exon.only = eqtls.exon[eqtls.exon$Gene_symbol %in% no.overlap, ]
# get top hit
exon.only[which.min(exon.only$FDR), ]

# load necessary mapping and expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_core_gene_summaries_annotated_combat_batch_corrected.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_exon_summaries_annotated_combat_batch_corrected.RData")
huex <- fread("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/HuEx-1_0-st-v2/HuEx-1_0-st-v2.na36.hg19.probeset.csv", data.table = F, skip = 22)
huex <- huex[huex$start != "---", ]  # remove ids with no start and stop information
gene.code <- as.data.frame(import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf"))
transcript.ids <- unique(huex[grep("USP24", huex$gene_assignment), "transcript_cluster_id"])
# get only the 68 exons of transcript 1 (USP24-201)
#exons <- c("ENSE00001370828.6", "ENSE00001382948.1", "ENSE00001613119.1", "ENSE00001659875.1", "ENSE00001706381.1", "ENSE00001756617.1", "ENSE00001747335.1", "ENSE00001601173.1", "ENSE00001623489.1", "ENSE00001687516.1", "ENSE00001744360.1", "ENSE00001780813.1", "ENSE00001805718.1", "ENSE00001454788.1", "ENSE00001701698.1", "ENSE00001454786.1", "ENSE00001629856.1", "ENSE00001728065.1", "ENSE00001639761.1", "ENSE00001803592.1", "ENSE00001454781.1", "ENSE00001454780.1", "ENSE00001454779.1", "ENSE00001454778.1", "ENSE00001454777.1", "ENSE00001454776.1", "ENSE00001454775.1", "ENSE00001454774.1", "ENSE00001454773.1", "ENSE00001454772.1", "ENSE00001454771.1", "ENSE00001454770.1", "ENSE00001454769.1", "ENSE00001454768.1", "ENSE00001454767.1", "ENSE00001454766.1", "ENSE00001454765.1", "ENSE00001454764.1", "ENSE00001434462.2", "ENSE00001433950.1", "ENSE00001415442.1", "ENSE00001156674.3", "ENSE00001139599.1", "ENSE00001139590.1", "ENSE00001139585.1", "ENSE00001139580.1", "ENSE00001139572.1", "ENSE00001139565.1", "ENSE00001139555.1", "ENSE00001139544.1", "ENSE00001139535.1", "ENSE00001064828.1", "ENSE00001139522.1", "ENSE00001139516.1", "ENSE00001139510.1", "ENSE00001139504.1", "ENSE00003648271.1", "ENSE00003571309.1", "ENSE00001139476.1", "ENSE00001139465.1", "ENSE00001139457.1", "ENSE00001064819.1", "ENSE00001139439.2", "ENSE00003681127.1", "ENSE00003481075.1", "ENSE00003575373.1", "ENSE00003529998.1", "ENSE00001847641.1")
exons <- unique(gene.code[gene.code$gene_name == "USP24" & gene.code$transcript_name == "USP24-001", "exon_id"])
exon.ids <- unique(gene.code[gene.code$gene_id == "ENSG00000162402.8", ]$exon_id)
exon.ids <- exon.ids[exon.ids %in% exons[complete.cases(exons)]]
transcript.exprs = t(exprs(geneSummaries[featureNames(geneSummaries) %in% transcript.ids, ])) # get expression of transcripts
exon.exprs <- t(exprs(exonSummaries[featureNames(exonSummaries) %in% exon.ids[-1], ])) # get expression of exons
both <- merge(transcript.exprs, exon.exprs, by = "row.names")
save(both, file = "USP24_gene_exon_expression.rda") 

#snps
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --snp rs662145 --recode --tab --out USP24_SNPs"))
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --snp rs473498 --recode --tab --out NCALD_SNPs"))
snp = read.delim("USP24_SNPs.ped", head = T)
pheno = pData(geneSummaries)[, c("IID","Dx")]

exp.p = merge(both, pheno, by.x = "Row.names", by.y = "row.names")
exp.p.snp = merge(exp.p, snp, by.x = "IID", by.y = "X1.1")
save(exp.p.snp, file = "USP24_gene_exon_expression_snp.rda")

plot.data <- melt(exp.p.snp, id.vars = "C.T", measure.vars = c("2413943", grep("ENS", colnames(exp.p.snp), value = T)))
plot.data <- plot.data[plot.data$C.T != "0 0",]
#plot.data[plot.data$rs662145_C == "0" & !is.na(plot.data$rs662145_C), "rs662145_C"] <- "CC"  # change 0 to major allele (CC)
#plot.data[plot.data$rs662145_C == "1" & !is.na(plot.data$rs662145_C), "rs662145_C"] <- "CT"  # change 1 to heterozygote allele (CT)
#plot.data[plot.data$rs662145_C == "2" & !is.na(plot.data$rs662145_C), "rs662145_C"] <- "TT"  # change 2 to minor allele (TT)
pdf("03_boxplots_USP24_gene_exon_expression_snp.pdf")
ggplot(na.omit(plot.data), aes(C.T, value)) +
  geom_boxplot() + facet_wrap(~variable) +
  theme(strip.text = element_text(size = 4.1), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  xlab("rs662145") + ylab("gene expression")
dev.off()

plot.data = melt(exp.p.snp, id.vars = "rs662145_C", measure.vars = c("2413943",grep("ENS", colnames(exp.p.snp), value = T)))
pdf("USP24_gene_exon_expression_snp.pdf")
ggplot(na.omit(plot.data), aes(as.factor(rs662145_C),value)) +
  geom_boxplot() + facet_wrap(~variable ) +
  theme(strip.text = element_text(size = 3)) # The new stuff
dev.off()
# ENSE00001811443.1 is also sig
exon.only[ exon.only$Exon == "ENSE00001811443.1" & exon.only$SNP == "rs662145", ]

### ANNOTATION ###
suppressPackageStartupMessages({
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(ChIPseeker)
library(org.Hs.eg.db)
})
#snps <- unlist(lapply(gene.eqtl.snps, function(x){if(!startsWith(x, "rs")) paste("chr", strsplit(x, ":")[[1]][1], ":", strsplit(x, ":")[[1]][2], "-", strsplit(x, ":")[[1]][2], sep = "") else x}))
# 1. Annotation on gene level
# 1.1 Generate GRanges object for all cis-eQTL SNPs
eqtls.gene$chr <- unlist(lapply(eqtls.gene$SNP_location, function(x) strsplit(x, ":")[[1]][1]))
eqtls.gene$start <- unlist(lapply(eqtls.gene$SNP_location, function(x) strsplit(x, ":")[[1]][2]))
eqtls.gene$SNP <- as.factor(eqtls.gene$SNP)
eqtls.gene.gr <- GRanges(seqnames = paste0("chr", eqtls.gene$chr), ranges = IRanges(start = as.numeric(eqtls.gene$start), width = 1))
mcols(eqtls.gene.gr) <- dplyr::select(eqtls.gene, -c(chr, start))
# 1.2 Annotate data
as.gene <- annotatePeak(unique(eqtls.gene.gr), 
                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                        annoDb = "org.Hs.eg.db", verbose = F)
# 2. Annotation on exon level
# 2.1 Generate GRanges object for all cis-eQTL SNPs
eqtls.exon$chr <- unlist(lapply(eqtls.exon$SNP_location, function(x) strsplit(x, ":")[[1]][1]))
eqtls.exon$start <- unlist(lapply(eqtls.exon$SNP_location, function(x) strsplit(x, ":")[[1]][2]))
eqtls.exon$SNP <- as.factor(eqtls.exon$SNP)
eqtls.exon.gr <- GRanges(seqnames = paste0("chr", eqtls.exon$chr), ranges = IRanges(start = as.numeric(eqtls.exon$start), width = 1))
mcols(eqtls.exon.gr) <- dplyr::select(eqtls.exon, -c(chr, start))
# 2.2 Annotate data
as.exon <- annotatePeak(unique(eqtls.exon.gr), 
                        TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, 
                        annoDb = "org.Hs.eg.db", verbose = F)
# 3. Plot the annotations
# 3.1 Combine the two annotations
genomic.annotation <- bind_rows(Gene = as.gene@annoStat,
                                Exon = as.exon@annoStat,
                                .id = 'Type')

# 3.2 Plot Annotations as stacked barplot
ggplot(genomic.annotation, aes(Type, Frequency, fill = Feature)) + 
  geom_bar(stat = 'identity', position = 'fill') + 
  scale_fill_brewer(palette = 'Paired', direction = -1) +
  coord_flip() + 
  labs(title = 'UCSC Genomic Features', y = 'Frequency') +
  theme_minimal() +
  theme(text = element_text(size = 22))
ggsave("/storage/scratch/users/karolina.worf/MPIP/05.eQTLs/Barplot_stacked_eqtl_annotations.pdf")

### LOAD Fromer et al. (CMC) eQTLs
fromer.eqtls.sva <- fread("/storage/groups/ccm01/datasets/2019_Fromer_et_al_Supplements/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt", sep = " ", stringsAsFactors = F, select = c("SNP", "Gene_symbol", "pvalue", "FDR"))
# get only hits (adj.p-val < 0.05)
fromer.eqtls.sva <- fromer.eqtls.sva[fromer.eqtls.sva$FDR < 0.05, ]
# remove all rows there gene symbol is not available symbolized as a dot "." or any rows with NA's in the data
fromer.eqtls.sva <- fromer.eqtls.sva[fromer.eqtls.sva$Gene_symbol != ".", ]
fromer.eqtls.sva <- fromer.eqtls.sva[complete.cases(fromer.eqtls.sva), ]

cmc.snps <- unique(sort(fromer.eqtls.sva$SNP))
gene.eqtl.snps <- unique(sort(eqtls.gene$SNP))
exon.eqtl.snps <- unique(sort(eqtls.exon$SNP))

# percentage of overlap
length(intersect(cmc.snps, gene.eqtl.snps))/(length(gene.eqtl.snps)/100)
length(intersect(cmc.snps, gene.eqtl.snps))/(length(unique(union(cmc.snps, gene.eqtl.snps))))

### LOAD gtex brain: frontal cortex eQTLs ###
gtex <- fread("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/Brain_Frontal_Cortex_BA9.v7.egenes.txt", header = T, sep = "\t", stringsAsFactors = F, select = c("rs_id_dbSNP147_GRCh37p13", "qval"))
gtex <- gtex[gtex$qval < 0.05, ]
gtex.snps <- unique(sort(gtex$rs_id_dbSNP147_GRCh37p13))
length(gtex.snps)

# percentage of overlap
length(unique(intersect(gtex.snps, gene.eqtl.snps)))/(length(gene.eqtl.snps)/100)
length(unique(intersect(gtex.snps, gene.eqtl.snps)))/length(unique(union(gtex.snps, gene.eqtl.snps)))

#library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#library(ChIPseeker)
#library(org.Hs.eg.db)


# SCATTERPLOT: eQTLs vs. own GWAS
# Load eQTLs
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results_old/01_pmb_SCZ_CTL_matrixEQTL.RData")
me.all.gene <- me.all
eqtls.gene <- me.all.gene$cis$eqtls[me.all.gene$cis$eqtls$FDR < 0.05, c("snps", "pvalue")] # filter for significant eQTLs
rm(me.all, me.all.gene)
# Load GWAS
gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc.logistic", header = T, stringsAsFactors = F, data.table = F, select =  c("SNP", "P"))
gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values
gwas.eqtl.snps <- merge(eqtls.gene, gwasResults, by.x = "snps", by.y = "SNP") # merge the two datasets (not overlapping SNPs will be removed from the dataset)
gwas.eqtl.snps$log_pvalue <- -log10(gwas.eqtl.snps$pvalue)
gwas.eqtl.snps$log_P <- -log10(gwas.eqtl.snps$P)
ggplot(gwas.eqtl.snps, aes(x = log_pvalue, y = log_P)) + 
  geom_point(alpha = 0.3) +
  ylab("-log10(P) of GWAS") +
  xlab("-log10(P) of eQTLs")

# SCATTERPLOT: eqTLs vs. SCZ GWAS
# Load GWAS
gwasResults <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics//SCZ2/scz2.snp.results.txt", header = T, stringsAsFactors = F, data.table = F, select = c("snpid", "p"))
gwas.eqtl.snps <- merge(eqtls.gene, gwasResults, by.x = "snps", by.y = "snpid") # merge the two datasets (not overlapping SNPs will be removed from the dataset)
gwas.eqtl.snps$log_pvalue <- -log10(gwas.eqtl.snps$pvalue)
gwas.eqtl.snps$log_p <- -log10(gwas.eqtl.snps$p)
ggplot(gwas.eqtl.snps, aes(x = log_pvalue, y = log_p)) + 
  geom_point(alpha = 0.3) +
  ylab("-log10(P) of GWAS") +
  xlab("-log10(P) of eQTLs")



