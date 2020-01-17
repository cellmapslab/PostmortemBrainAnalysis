# add installed packages from home path, to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)
library(huex10sttranscriptcluster.db) # transcript/gene infos for HuEx-1_0-st-v2
library(AnnotationDbi)
library(biomaRt)
library(rtracklayer)

############################################################
### 1. Gene (based on core transcripts) eQTL Input Files ###
############################################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
gene.coreSummaries <- geneSummaries
# check dim of gene.coreSummaries (14944 features and 169 samples)
dim(gene.coreSummaries)
pData <- pData(gene.coreSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)
# remove uneccesary data
rm(geneSummaries)

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
write.table(snpspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/SNP_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

##########################
### GET GENE POSITIONS ###
##########################
# Map gene symbols to chr, minimum start and maximum end position and save them into a file
genes.fdata <- gene.coreSummaries@featureData@data[, c("gene_id", "location")]
genespos <- data.frame()
for (i in 1:length(genes.fdata[, 1])) {
  name <- genes.fdata[i, "gene_id"]
  chr <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", genes.fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(genes.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(genes.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  genespos <- rbind(genespos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
colnames(genespos) <- c("geneid", "chr", "start", "end")
head(genespos)

# check, if we now really have the same genes with locations
identical(as.vector(genespos$geneid), fData(gene.coreSummaries)$gene_id)
# save data into a file
write.table(genespos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/gene_core_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)
rm(genes.fdata, row, chr, i, max, min, name)

###################
### ALL SAMPLES ###
###################
# 1. Generate genotype/SNP dataset
# Get SNP data (throws a warning, just ignore it!!)
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
identical(sort(colnames(exprs(gene.coreSummaries))), colnames(snps))  # NOT TRUE => Problem!
head(snps)
# Save SNP data into a new file
fwrite(snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_snps_CASE_CTL.txt", row.names = T, col.names = T, sep = "\t", quote = F)


# 2. Generate gene expression dataset
# Get Gene expression data
expr <- exprs(gene.coreSummaries)
rownames(expr) <- gene.coreSummaries@featureData@data$gene_id
expr[1:6, 1:6]

# check if gene ID order is the same as in genespos
identical(rownames(expr), as.vector(genespos$geneid))# Gene expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for exon expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save gene expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_gene_core_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)


# 3. Generate Covariate Dataset
# Covariates: get Dx, Sex, Age, pH, PMI, RIN, RIN², Suicide, Type_of_Death, C1, C2, C3, C4 (C1-4: Ancestry information)
rownames(pData) <- pData$ShortName
pData["RIN2"] <- pData$RIN^2

# generate model
cov <- model.matrix(~ Dx + Sex + Age + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
cov <- t(cov[,-1])
#cov <- apply(cov, 2, function(x) round(x, digits = 4))
cov[1:14, 1:6]
cov <- cov[, order(colnames(cov))]
dim(cov)

# Covariates: check if samples overlap to snps and gene expression
identical(colnames(cov), colnames(snps))
cov <- data.table(cov, keep.rownames = T)

# Save covariates into a new file
fwrite(cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_core_cov_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

##########################
### ONLY CASE SAMPLES ###
#########################
# get all cases (107 samples)
cov2 <- as.data.frame(cov)  # change data.table to data.frame
rownames(cov2) <- cov2$rn # get rownames
cov2$rn <- NULL # remove rn column
cov2 <- as.data.frame(t(cov2)) # transpose data.frame
case.samples <- row.names(cov2[cov2$DxCASE == 1, ]) # get case samples
case.samples

# get all covariates with these samples
case.cov.samples <- c("rn", case.samples)
case.cov <- subset(cov, , case.cov.samples)
# remove Case rows
case.cov <- case.cov[case.cov$rn != "DxCASE"]
head(case.cov)
# save covariates into a file
fwrite(case.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_core_cov_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all snps with these samples
case.snps <- subset(snps, , case.samples)
head(case.snps)
identical(colnames(case.snps), colnames(case.cov[, -1]))
# save snps into a file
fwrite(case.snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_snps_CASE_only.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# get gene expression data for these samples
case.samples <- c("rn", case.samples)
case.expr <- subset(expr, , case.samples)
head(case.expr)
identical(colnames(case.expr[, -1]), colnames(case.cov[, -1]))
identical(colnames(case.expr[, -1]), colnames(case.snps))
# save gene expression into a file
fwrite(case.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_gene_core_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get all CTL cases (62 samples)
cov3 <- as.data.frame(cov)  # change data.table to data.frame
rownames(cov3) <- cov3$rn # get rownames
cov3$rn <- NULL # remove rn column
cov3 <- as.data.frame(t(cov3)) # transpose data.frame
ctl.samples <- row.names(cov3[cov3$DxCASE == 0, ]) # get ctl samples

# get all covariates with these samples
ctl.cov.samples <- c("rn", ctl.samples)
ctl.cov <- subset(cov, , ctl.cov.samples)
# remove both Dx (Cases) rows and Suicide (as no CTL made suicide & matrixEQTL does not like cov without any changes!!)
ctl.cov <- ctl.cov[ctl.cov$rn != "DxCASE" & ctl.cov$rn != "SuicideY"]
head(ctl.cov)
# save covariates into a file
fwrite(ctl.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_core_cov_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all snps with these samples
ctl.snps <- subset(snps, , ctl.samples)
head(ctl.snps)
identical(colnames(ctl.snps), colnames(ctl.cov[, -1]))
# save snps into a file
fwrite(ctl.snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_snps_CTL_only.txt", row.names = T, col.names = T, sep = "\t", quote = F)

# get gene expression data for these samples
ctl.samples <- c("rn", ctl.samples)
ctl.expr <- subset(expr, , ctl.samples)
head(ctl.expr)
identical(colnames(ctl.expr[, -1]), colnames(ctl.cov[, -1]))
identical(colnames(ctl.expr[, -1]), colnames(ctl.snps))
# save gene expression into a file
fwrite(ctl.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_gene_core_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)


#####################################################
### 2. Gene (based on probesets) eQTL Input Files ###
#####################################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
gene.probesetSummaries <- geneSummaries

# check dim of probesetSummaries (17447 features and 169 samples)
dim(gene.probesetSummaries)

# get new pData (including the SV1 variable)
pData <- pData(gene.probesetSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

##########################
### GET GENE POSITIONS ###
##########################
# Map gene symbols to chr, minimum start and maximum end position and save them into a file
genes.fdata <- gene.probesetSummaries@featureData@data[, c("gene_id", "location")]
genespos <- data.frame()
for (i in 1:length(genes.fdata[, 1])) {
  name <- genes.fdata[i, "gene_id"]
  chr <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", genes.fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(genes.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(genes.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(genes.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  genespos <- rbind(genespos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
colnames(genespos) <- c("geneid", "chr", "start", "end")
head(genespos)

# check, if we now really have the same genes with locations
identical(as.vector(genespos$geneid), fData(gene.probesetSummaries)$gene_id)
# save data into a file
write.table(genespos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/gene_probeset_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# 2. Generate Gene Expression Dataset
# Get Gene expression data
expr <- exprs(gene.probesetSummaries)
expr[1:6, 1:6]

# check if exon ID order is the same as in genespos
identical(as.character(rownames(expr)), as.character(genespos$geneid))# Transcript expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save gene probeset expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_gene_probeset_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# 3. Generate Covariate Dataset
# Covariates: get Dx, Sex, Age, pH, PMI, RIN, RIN², Suicide, Type_of_Death, C1, C2, C3, C4 (C1-4: Ancestry information)
rownames(pData) <- pData$ShortName

# generate model
cov <- model.matrix(~ Dx + Sex + Age + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
cov <- t(cov[,-1])
#cov <- apply(cov, 2, function(x) round(x, digits = 4))
cov[1:15, 1:6]
cov <- cov[, order(colnames(cov))]
dim(cov)

# Covariates: check if samples overlap to snps and gene expression
identical(colnames(cov), colnames(snps))
cov <- data.table(cov, keep.rownames = T)

# Save covariates into a new file
fwrite(cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_probeset_cov_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)


#########################
### ONLY CASE SAMPLES ###
#########################
# get  gene probeset expression data for SCZ (CASE) samples
case.gene.probeset.expr <- subset(expr, , case.samples)
head(case.gene.probeset.expr)
identical(colnames(case.gene.probeset.expr), colnames(case.expr))
# savegene probeset expression into a file
fwrite(case.gene.probeset.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_gene_probeset_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all covariates with these samples
case.cov.samples <- c("rn", case.samples)
case.cov <- subset(cov, , case.cov.samples)
head(case.cov)
# remove Case rows
case.cov <- case.cov[case.cov$rn != "DxCASE"]
head(case.cov)
# save covariates into a file
fwrite(case.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_probeset_cov_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get gene probeset expression data for CTL samples
ctl.gene.probeset.expr <- subset(expr, , ctl.samples)
head(ctl.gene.probeset.expr)
identical(colnames(ctl.gene.probeset.expr), colnames(ctl.expr))
# save gene probeset expression into a file
fwrite(ctl.gene.probeset.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_gene_probeset_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

# get all covariates with these samples
ctl.cov.samples <- c("rn", ctl.samples)
ctl.cov <- subset(cov, , ctl.cov.samples)
head(ctl.cov)
# remove both Dx (Cases) rows and Suicide (as no CTL made suicide & matrixEQTL does not like cov without any changes!!)
ctl.cov <- ctl.cov[ctl.cov$rn != "DxCASE" & ctl.cov$rn != "SuicideY"]
head(ctl.cov)
# save covariates into a file
fwrite(ctl.cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_probeset_cov_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

######################################
### 3. Transcript eQTL Input Files ###
######################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")

# check dim of transcriptSummaries (35 rows should be removed => 169 samples left)
dim(transcriptSummaries)
pData <- pData(transcriptSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

################################
### GET Transcript POSITIONS ###
################################
# Map gene symbols to chr, minimum start and maximum end position and save them into a file
transcript.fdata <- transcriptSummaries@featureData@data[, c("transcript_id", "location")]
transpos <- data.frame()
for (i in 1:length(transcript.fdata[, 1])) {
  name <- rownames(transcript.fdata[i, ])
  chr <- strsplit(strsplit(as.character(transcript.fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", transcript.fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(transcript.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(transcript.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(transcript.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(transcript.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  transpos <- rbind(transpos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}

transpos <- as.data.frame(transpos)
colnames(transpos) <- c("transcriptid", "chr", "start", "end")
head(transpos)

# check, if we now really have the same genes with locations
identical(as.vector(transpos$transcriptid), as.vector(fData(transcriptSummaries)$transcript_id))
# save data into a file
write.table(transpos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/transcript_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# 2. Generate transcript Expression Dataset
# Get transcript expression data
expr <- exprs(transcriptSummaries)
expr[1:6, 1:6]
# check if exon ID order is the same as in transpos
identical(unique(rownames(expr)), as.character(transpos$transcriptid))# Transcript expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save transcript expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_transcript_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get transcript expression data for SCZ (CASE) samples
case.transcript.expr <- subset(expr, , case.samples)
head(case.transcript.expr)
identical(colnames(case.transcript.expr), colnames(case.expr))
# save gene expression into a file
fwrite(case.transcript.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_transcript_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get gene expression data for CTL samples
ctl.transcript.expr <- subset(expr, , ctl.samples)
head(ctl.transcript.expr)
# save gene expression into a file
fwrite(ctl.transcript.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_transcript_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)


################################
### 3. Exon eQTL Input Files ###
################################
# Load exon expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")

# check dim of exonSummaries (35 rows should be removed => 169 samples left)
dim(exonSummaries)
pData <- pData(exonSummaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

##########################
### GET EXON POSITIONS ###
##########################
# Get exon symbols, chr, start and end pos
exons.fdata <- exonSummaries@featureData@data[, c("exon_id", "location")]
exonspos <- data.frame()
for (i in 1:length(exons.fdata[, 1])) {
  name <- rownames(exons.fdata[i, ])
  chr <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", exons.fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(exons.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(exons.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  exonspos <- rbind(exonspos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
colnames(exonspos) <- c("exonid", "chr", "start", "end")
head(exonspos)



exons.fdata <- exonSummaries@featureData@data[, c("exon_id", "location")]
#exonspos <- data.frame()
# parallel version
require(foreach)
require(doSNOW)
# set up cluster
threads <- parallel::detectCores()
cl <- makeCluster(threads, type = "SOCK")
registerDoSNOW(cl)
exonspos <- foreach(i = seq.int(exons.fdata[, 1]), .combine = rbind) %dopar% {
  name <- rownames(exons.fdata[i, ])
  chr <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", exons.fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(exons.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(exons.fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(exons.fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  # transpos <- rbind(transpos, row)
  if ((i %% 1000) == 0) {
    print(i)
  }
  return(row)
}

stopCluster(cl)

exonspos <- as.data.frame(exonspos)
colnames(exonspos) <- c("exonid", "chr", "start", "end")
head(exonspos)

# check, if we now really have the same exons with locations
identical(as.vector(exonspos$exonid), rownames(exprs(exonSummaries)))
# save data into a file
write.table(exonspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/exon_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)

###################
### ALL SAMPLES ###
###################
# Get exon expression data
expr <- exprs(exonSummaries)
expr[1:6, 1:6]
# check if exon ID order is the same as in exonspos
identical(unique(expr$rn), as.character(exonspos$exonid))# Gene expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for exon expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
head(expr)
identical(colnames(expr)[-1], cols)

# Save gene expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/01_pmb_exon_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get exon expression data for SCZ (CASE) samples
case.exon.expr <- subset(expr, , case.samples)
head(case.exon.expr)
identical(colnames(case.exon.expr), colnames(case.expr))

# save gene expression into a file
fwrite(case.exon.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/02_pmb_exon_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)


########################
### ONLY CTL SAMPLES ###
########################
# get gene expression data for CTL samples
ctl.exon.expr <- subset(expr, , ctl.samples)
head(ctl.exon.expr)
identical(colnames(ctl.exon.expr), colnames(ctl.expr))

# save gene expression into a file
fwrite(ctl.exon.expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/03_pmb_exon_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)
