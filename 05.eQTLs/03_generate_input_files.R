# add installed packages from home path, to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

library(data.table)
library(huex10sttranscriptcluster.db) # transcript/gene infos for HuEx-1_0-st-v2
library(AnnotationDbi)
library(biomaRt)
library(rtracklayer)

################################
### 1. Gene eQTL Input Files ###
################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_gene_summaries_filtered_annotated_batch_corrected.RData")
# As we have some missing values in the phenotype data, 33 mds (C1-C10) and 2 Type_of_Death values, 
# we first need to remove them from the microarrays as well as the phenotype dataset.
# remove all samples with C1 = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$C1))
geneSummaries <- geneSummaries[, -c(index)]
# remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$Type_of_Death))
geneSummaries <- geneSummaries[, -c(index)]
# check dim of geneSummaries (35 rows should be removed => 169 samples left)
dim(geneSummaries)
pData <- pData(geneSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

#########################
### GET SNP POSITIONS ###
#########################
# get SNP ID, chr and pos and save it
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", header = F, stringsAsFactors = FALSE) # snp chr pos
#snpspos[, 2] <- paste0(snpspos[, 2], ":", snpspos[, 6], ":", snpspos[, 5])
head(snpspos)
snpspos <- snpspos[, c("V2", "V1", "V4")]
colnames(snpspos) <- c("snp", "chr", "pos")
snpspos <- snpspos[!duplicated(snpspos$snp),]
head(snpspos)
write.table(snpspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/SNP_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

##########################
### GET GENE POSITIONS ###
##########################
# Map gene symbols to chr, start and end position and save them into a file
# Get gene symbols from our geneSummaries (microarray gene expression dataset)
symbols <- as.data.frame(geneSummaries@featureData@data$SYMBOL)
colnames(symbols) <- "symbols"
genes <- import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf")
#table(genes$type)
genes <- genes[genes$gene_type == "protein_coding"]  # get only protein-coding genes
genes <- genes[genes$type == "gene"]  # get only exons from the genes
genes <- as.data.frame(genes)
head(genes)
genes <- genes[, c("gene_name", "seqnames", "start", "end")]
geneinfo <- merge(symbols, genes, by.x = "symbols", by.y = "gene_name", all.x = T)
geneinfo <- geneinfo[!duplicated(geneinfo), ]
# But we have some duplicates (190 genes) now:
dups <- geneinfo[duplicated(geneinfo$symbols), ]$symbols
# remove duplicated genes with smaller length (most times one of the genes lies within the other one)
dup.to.remove <- NULL
for(gene in dups){
  length = 0
  info <- geneinfo[geneinfo$symbols == gene, ]
  info$length <- info$end - info$start
  dup.to.remove <- c(dup.to.remove, as.numeric(row.names(info[info$length == min(info$length), ])))
}
geneinfo <- geneinfo[!(rownames(geneinfo) %in% dup.to.remove), ]
geneinfo$chr <- sapply(strsplit(as.character(geneinfo$seqnames),'r'), "[", 2)
geneinfo <- geneinfo[, c("symbols", "chr", "start", "end")]
colnames(geneinfo) <- c("geneid", "chr", "start", "end")
# check, if we now really have the same genes with locations
identical(geneinfo$geneid, unique(sort(symbols$symbols)))
# save data into a file
write.table(geneinfo, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/gene_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# 1. Generate SNP Dataset
# Get SNP data
snps <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_probs_changed/pmb_info06_all_chr_maf0.05_merged_matrixEQTL_transformed.gen.gz" , data.table = F)

# remove duplicated SNPs and change rownames to SNP IDs
snps <- snps[!duplicated(snps$V1), ]
row.names(snps) <- snps$V1  # get row names <- DOES NOT WORK WITH DUPLICATED ROW.NAMES!!! 
snps$V1 <- NULL  # remove column 1 (including rownames)

# Map Sample IID's to short sample names
ids <- pData[, c("IID", "ShortName")]
ids <- ids[order(ids$IID), ]
cols <- as.data.frame(colnames(snps))
snps <- snps[, colnames(snps) %in% ids$IID]  # filter samples

# check if colnames of snps are identical to the IIDs given in pData
identical(as.character(ids$IID), colnames(snps))
# change snp column names (IID's) to short sample names
colnames(snps) <- ids$ShortName
snps <- snps[, order(colnames(snps))]
identical(sort(colnames(exprs(geneSummaries))), colnames(snps))  # NOT TRUE => Problem!
head(snps)
# Save SNP data into a new file
fwrite(snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_snps_SCZ_CTL.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# 2. Generate Gene Expression Dataset
# Get Gene expression data
expr <- exprs(geneSummaries)
rownames(expr) <- geneSummaries@featureData@data$SYMBOL
expr[1:6, 1:6]

# Get only genes also present in the geneinfo dataset and order them the same way
expr <- expr[rownames(expr) %in% geneinfo$geneid,]
expr <- expr[order(rownames(expr)), ]
# check if gene name order is the same
identical(unique(rownames(expr)), as.character(geneinfo$geneid))# Gene expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save gene expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_gene_expr_SCZ_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# 3. Generate Covariate Dataset
# Covariates: get Dx, Sex, Age, pH, PMI, RIN, RIN², Suicide, Type_of_Death, C1, C2, C3, C4 (C1-4: Ancestry information)
rownames(pData) <- pData$ShortName
pData["RIN2"] <- pData$RIN^2
pData$Diagnosis <- pData$Dx
head(pData)

# chance BPD and MDD diagnosis to SCZ
pData$Dx <- as.character(pData$Dx)
pData[pData$Dx == "BPD" | pData$Dx == "MDD", "Dx"] <- "SCZ"
pData$Dx <- as.factor(pData$Dx)

# relevel disease (Dx) to control (CTL)
pData$Dx <- factor(pData$Dx, levels = c("CTL", "SCZ"))
cov <- model.matrix(~ Dx + Sex + Age + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
cov <- t(cov[,-1])
#cov <- apply(cov, 2, function(x) round(x, digits = 4))
cov[1:14, 1:6]
cov <- cov[, order(colnames(cov))]
dim(cov)

# Covariates: check if samples overlap to snps and gene expression
identical(colnames(cov[,-1]), colnames(snps))
cov <- data.table(cov, keep.rownames = T)

# Save covariates into a new file
fwrite(cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_cov_SCZ_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

##########################
### ONLY CASE SAMPLES ###
#########################
# get all cases (107 samples)
cov2 <- as.data.frame(cov)  # change data.table to data.frame
rownames(cov2) <- cov2$rn # get rownames
cov2$rn <- NULL # remove rn column
cov2 <- as.data.frame(t(cov2)) # transpose data.frame
scz.samples <- row.names(cov2[cov2$DxSCZ == 1, ]) # get scz samples
scz.samples

# get all covariates with these samples
scz.cov.samples <- c("rn", scz.samples)
scz.cov <- subset(cov, , scz.cov.samples)
# remove both Dx (AFF and SCZ) rows
scz.cov <- scz.cov[scz.cov$rn != "DxAFF" & scz.cov$rn != "DxSCZ"]
head(scz.cov)
# save covariates into a file
fwrite(scz.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_cov_SCZ.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all snps with these samples
scz.snps <- subset(snps, , scz.samples)
head(scz.snps)
# save snps into a file
fwrite(scz.snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_snps_SCZ.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# get gene expression data for these samples
scz.samples <- c("rn", scz.samples)
scz.expr <- subset(expr, , scz.samples)
head(scz.expr)
# save gene expression into a file
fwrite(scz.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_gene_expr_SCZ.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get all CTL cases (62 samples)
cov3 <- as.data.frame(cov)  # change data.table to data.frame
rownames(cov3) <- cov3$rn # get rownames
cov3$rn <- NULL # remove rn column
cov3 <- as.data.frame(t(cov3)) # transpose data.frame
ctl.samples <- row.names(cov3[cov3$DxSCZ == 0, ]) # get ctl samples

# get all covariates with these samples
ctl.cov.samples <- c("rn", ctl.samples)
ctl.cov <- subset(cov, , ctl.cov.samples)
# remove both Dx (AFF and SCZ) rows and Suicide (as no CTL made suicide & matrixEQTL does not like cov without any changes!!)
ctl.cov <- ctl.cov[ctl.cov$rn != "DxAFF" & ctl.cov$rn != "DxSCZ" & ctl.cov$rn != "SuicideY"]
head(ctl.cov)
# save covariates into a file
fwrite(ctl.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_cov_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all snps with these samples
ctl.snps <- subset(snps, , ctl.samples)
head(ctl.snps)
# save snps into a file
fwrite(ctl.snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_snps_CTL.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# get gene expression data for these samples
ctl.samples <- c("rn", ctl.samples)
ctl.expr <- subset(expr, , ctl.samples)
head(ctl.expr)
# save gene expression into a file
fwrite(ctl.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_gene_expr_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)



################################
### 2. Exon eQTL Input Files ###
################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")

# As we have some missing values in the phenotype data, 33 mds (C1-C10) and 2 Type_of_Death values, 
# we first need to remove them from the microarrays as well as the phenotype dataset.
# remove all samples with C1 = NA from phenotype data and gene expression data
index <- which(is.na(pData(exonSummaries)$C1))
exonSummaries <- exonSummaries[, -c(index)]
# remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(exonSummaries)$Type_of_Death))
exonSummaries <- exonSummaries[, -c(index)]
# check dim of exonSummaries (35 rows should be removed => 169 samples left)
dim(exonSummaries)
pData <- pData(exonSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

##########################
### GET EXON POSITIONS ###
##########################
# Get exon symbols, chr, start and end pos
exonspos <- exonSummaries@featureData@data[, c("exon_id", "location")]
exonspos$chr <- sapply(strsplit(exonspos$location, ':'), "[", 1)
exonspos$chr <- sapply(strsplit(exonspos$chr, 'r'), "[", 2)
exonspos$loc <- sapply(strsplit(exonspos$location, ':'), "[", 2)
exonspos$start <- sapply(strsplit(exonspos$loc, '-'), "[", 1)
exonspos$end <- sapply(strsplit(exonspos$loc, '-'), "[", 2)
exonspos <- exonspos[, c("exon_id", "chr", "start", "end")]
colnames(exonspos) <- c("geneid", "chr", "start", "end")
head(exonspos)
# check, if we now really have the same exons with locations
identical(exonspos$geneid, rownames(exprs(exonSummaries)))
# save data into a file
write.table(exonspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/exon_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# 2. Generate Exon Expression Dataset
# Get exon expression data
expr <- exprs(exonSummaries)
expr[1:6, 1:6]
# check if exon ID order is the same as in exonspos
identical(unique(rownames(expr)), as.character(exonspos$geneid))# Gene expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for exon expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save gene expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_exon_expr_SCZ_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get exon expression data for SCZ (CASE) samples
scz.exon.expr <- subset(expr, , scz.samples)
head(scz.exon.expr)

# save gene expression into a file
fwrite(scz.exon.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_exon_expr_SCZ.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get gene expression data for CTL samples
ctl.exon.expr <- subset(expr, , ctl.samples)
head(ctl.exon.expr)
# save gene expression into a file
fwrite(ctl.exon.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_exon_expr_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)


######################################
### 3. Transcript eQTL Input Files ###
######################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")

# As we have some missing values in the phenotype data, 33 mds (C1-C10) and 2 Type_of_Death values, 
# we first need to remove them from the microarrays as well as the phenotype dataset.
# remove all samples with C1 = NA from phenotype data and gene expression data
index <- which(is.na(pData(transcriptSummaries)$C1))
transcriptSummaries <- transcriptSummaries[, -c(index)]
# remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(transcriptSummaries)$Type_of_Death))
transcriptSummaries <- transcriptSummaries[, -c(index)]
# check dim of transcriptSummaries (35 rows should be removed => 169 samples left)
dim(transcriptSummaries)
pData <- pData(transcriptSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

################################
### GET Transcript POSITIONS ###
################################
# Get exon symbols, chr, start and end pos
transpos <- transcriptSummaries@featureData@data[, c("transcript_id", "location")]
transpos$chr <- sapply(strsplit(transpos$location, ':'), "[", 1)
transpos$chr <- sapply(strsplit(transpos$chr, 'r'), "[", 2)
transpos$loc <- sapply(strsplit(transpos$location, ':'), "[", 2)
transpos$start <- sapply(strsplit(transpos$loc, '-'), "[", 1)
transpos$end <- sapply(strsplit(transpos$loc, '-'), "[", 2)
transpos <- transpos[, c("transcript_id", "chr", "start", "end")]
colnames(transpos) <- c("geneid", "chr", "start", "end")
head(transpos)
# check, if we now really have the same exons with locations
identical(transpos$geneid, rownames(exprs(transcriptSummaries)))
# save data into a file
write.table(transpos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/transcript_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# 2. Generate transcript Expression Dataset
# Get transcript expression data
expr <- exprs(transcriptSummaries)
expr[1:6, 1:6]
# check if exon ID order is the same as in transpos
identical(unique(rownames(expr)), as.character(transpos$geneid))# Transcript expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save transcript expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_transcript_expr_SCZ_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get transcript expression data for SCZ (CASE) samples
scz.transcript.expr <- subset(expr, , scz.samples)
head(scz.transcript.expr)
# save gene expression into a file
fwrite(scz.transcript.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/02_pmb_transcript_expr_SCZ.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get gene expression data for CTL samples
ctl.transcript.expr <- subset(expr, , ctl.samples)
head(ctl.transcript.expr)
# save gene expression into a file
fwrite(ctl.transcript.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/03_pmb_transcript_expr_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)
