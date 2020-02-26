########################################
### Generate Matrix eQTL input files ###
########################################
# expression data is always coming from the probeset-based summarization (NO core summarization for the genes)

# add installed packages from home path, to the library search path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)
library(huex10sttranscriptcluster.db) # transcript/gene infos for HuEx-1_0-st-v2
library(AnnotationDbi)
library(biomaRt)
library(rtracklayer)

#####################
### SNP LOCATIONS ###
######################
# get SNP ID, chromosome (chr) and bp position (pos) and save it into a file
snpspos <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.4_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup.bim", header = FALSE, stringsAsFactors = FALSE) # snp chr pos
snpspos <- snpspos[, c("V2", "V1", "V4")] # get only relevant columns
colnames(snpspos) <- c("snp", "chr", "pos") # change colnames
snpspos <- snpspos[!duplicated(snpspos$snp),] # remove duplicates
snpspos$chr <- sapply(snpspos$chr, function(x) paste0("chr", x))  # add chr to chromosome number
dim(snpspos)  # 6,802,325 SNPs and 3 columns
snpspos[1:3, ]
write.table(snpspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/SNP_locations.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

########################################
### 1. Gene eQTL (geQTL) Input Files ###
########################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
dim(gene_summaries) # 17,447 features and 169 samples

# get phenotype data (including the SV1 variable)
pData <- pData(gene_summaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

######################
### GENE LOCATIONS ###
######################
# Map gene Ensembl IDs to chr, minimum start and maximum end position and save them into a file
genes_fdata <- gene_summaries@featureData@data[, c("gene_id", "location")]
genespos <- data.frame()
for (i in 1:length(genes_fdata[, 1])) {
  name <- genes_fdata[i, "gene_id"]
  chr <- strsplit(strsplit(as.character(genes_fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", genes_fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(genes_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(genes_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(genes_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(genes_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  genespos <- rbind(genespos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
colnames(genespos) <- c("geneid", "chr", "start", "end")
genespos$chr <- sapply(genespos$chr, function(x) paste0("chr", x)) # add chr to chromosome number
dim(genespos)  # 17,447 genes and four columns
genespos[1:3, ]

# check, if we now really have the same genes with locations
identical(as.vector(genespos$geneid), fData(gene_summaries)$gene_id)
# save data into a file
write.table(genespos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/gene_locations.txt", col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
rm(genes_fdata, i, name, chr, min, max, row)

###################
### ALL SAMPLES ###
###################
# 1. Generate genotype (SNP) dataset
# Get SNP data (throws a warning, just ignore it!!)
snps <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.5_probs_changed/pmb_info06_all_chr_maf0.05_merged_noDup_matrixEQTL_transformed.gen" , data.table = FALSE)

# remove duplicated SNPs and change rownames to SNP IDs
snps <- snps[!duplicated(snps$V1), ]
row.names(snps) <- snps$V1  # get row names <- DOES NOT WORK WITH DUPLICATED ROW.NAMES!!! 
snps$V1 <- NULL  # remove column 1 (saved as rownames now)

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
identical(sort(colnames(exprs(gene_summaries))), colnames(snps))  # NOT TRUE => Problem!
head(snps)
# check if the SNPs are in the same order as the SNP location (snpspos) file
identical(snpspos$V2, row.names(snps))
snps$id <- row.names(snps)
snps <- snps[, c(length(snps), 1:length(snps)-1)]
# Save SNP data into a new file
fwrite(snps, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_snps_CASE_CTL.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# 2. Generate Gene Expression Dataset
# Get Gene expression data
expr <- exprs(gene_summaries)
expr[1:3, 1:6]
# check if gene ID order is the same as in genespos
identical(as.character(rownames(expr)), as.character(genespos$geneid))

# Transcript expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
colnames(expr)[1] <- "id"
identical(colnames(expr), colnames(snps))

# Save gene (probeset based) expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_gene_expr_CASE_CTL.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# 3. Generate Covariate Dataset
# Covariates: get Dx, Sex, Age, pH, PMI, RIN, RIN², Suicide, Type_of_Death, C1, C2, C3, C4 (C1-4: Ancestry information)
rownames(pData) <- pData$ShortName

# generate model
cov <- model.matrix(~ Dx + Sex + Age + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
cov <- t(cov[,-1])
cov[1:2, 1:5]
cov <- as.data.frame(cov[, order(colnames(cov))])
cov$id <- c("Diagnosis", "Sex", "Age", "pH", "PMI", "RIN", "RIN2", "Suicide", "Type_of_Death_Non-violent", "Type_of_Death_Violent", "C1", "C2", "C3", "C4", "SV1")
cov <- cov[, c(length(cov), 1:length(cov) - 1)]
dim(cov)

# Covariates: check if samples overlap to snps and gene expression
identical(colnames(cov), colnames(snps))
cov <- data.table(cov)

# Save covariates into a new file
fwrite(cov, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_cov_CASE_CTL.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
rm(ids, cols)

#########################
### CASE SAMPLES ONLY ###
#########################
# get all cases (107 samples) including the diagnosis for SCZ, BPD and MDD (cross-disorder)
cov_case_only <- as.data.frame(cov) # get covariates
rownames(cov_case_only) <- cov_case_only$id # get rownames
cov_case_only$id <- NULL # remove rn column
cov_case_only <- as.data.frame(t(cov_case_only)) # transpose data.frame
case_samples <- row.names(cov_case_only[cov_case_only$Diagnosis == 1, ]) # get case samples
length(case_samples)
case_samples <- c("id", case_samples)
case_samples

# get gene (probeset based) expression data for CASE samples
expr_case_only <- subset(expr, select = case_samples)
head(expr_case_only)
identical(colnames(expr_case_only), case_samples)
# savegene probeset expression into a file
fwrite(expr_case_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_gene_expr_CASE_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# get all covariates of CASE samples and remove diagnosis (all are the same)
cov_case_only <- subset(cov, select = case_samples)
head(cov_case_only)
# remove Diagnosis row
cov_case_only <- cov_case_only[cov_case_only$id != "Diagnosis", ]
head(cov_case_only)
# save covariates into a file
fwrite(cov_case_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_cov_CASE_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# get all snps with these samples
snps_case_only <- subset(snps, select = case_samples)
head(snps_case_only)
identical(colnames(snps_case_only), colnames(cov_case_only))
# save snps into a file
fwrite(snps_case_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_snps_CASE_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
rm(cov_case_only, expr_case_only, snps_case_only)

########################
### CTL SAMPLES ONLY ###
########################
# get gene probeset expression data for control (CTL) samples
# get all controls (62 samples) samples
cov_ctl_only <- as.data.frame(cov) # get covariates
rownames(cov_ctl_only) <- cov_ctl_only$id # get rownames
cov_ctl_only$id <- NULL # remove rn column
cov_ctl_only <- as.data.frame(t(cov_ctl_only)) # transpose data.frame
ctl_samples <- row.names(cov_ctl_only[cov_ctl_only$Diagnosis == 0, ]) # get control (CTL) samples
length(ctl_samples)
ctl_samples <- c("id", ctl_samples)
ctl_samples

# get gene (probeset based) expression data for CTL samples
expr_ctl_only <- subset(expr, select = ctl_samples)
head(expr_ctl_only)
identical(colnames(expr_ctl_only), ctl_samples)
# savegene probeset expression into a file
fwrite(expr_ctl_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_gene_expr_CTL_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# get all covariates of CTL samples and remove diagnosis (all are the same)
cov_ctl_only <- subset(cov, select = ctl_samples)
head(cov_ctl_only)
# remove Diagnosis and Suicide row
cov_ctl_only <- cov_ctl_only[cov_ctl_only$id != "Diagnosis", ]
cov_ctl_only <- cov_ctl_only[cov_ctl_only$id != "Suicide", ]
head(cov_ctl_only)
# save covariates into a file
fwrite(cov_ctl_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_cov_CTL_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

# get all snps with these samples
snps_ctl_only <- subset(snps, select = ctl_samples)
head(snps_ctl_only)
identical(colnames(snps_ctl_only), colnames(cov_ctl_only))
# save snps into a file
fwrite(snps_ctl_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_snps_CTL_only.txt", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
rm(cov_ctl_only, expr_ctl_only, snps_ctl_only, gene_summaries, pData)


##############################################
### 2. Transcript eQTL (teQTL) Input Files ###
##############################################
# Load gene expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
dim(transcript_summaries) # 100,750 features and 169 samples

# get phenotype data (including the SV1 variable)
pData <- pData(transcript_summaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

############################
### Transcript LOCATIONS ###
############################
# Map transcript Ensembl IDs to chr, minimum start and maximum end position and save them into a file
transcript_fdata <- transcript_summaries@featureData@data[, c("transcript_id", "location")]
transpos <- data.frame()
for (i in 1:length(transcript_fdata[, 1])) {
  name <- rownames(transcript_fdata[i, ])
  chr <- strsplit(strsplit(as.character(transcript_fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", transcript_fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(transcript_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(transcript_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(transcript_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(transcript_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  transpos <- rbind(transpos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
transpos <- as.data.frame(transpos)
colnames(transpos) <- c("transcriptid", "chr", "start", "end")
transpos$chr <- sapply(transpos$chr, function(x) paste0("chr", x)) # add chr to chromosome number
dim(transpos)  # 100,750 transcripts and four columns
transpos[1:3, ]

# check, if we now really have the same genes with locations
identical(as.vector(transpos$transcriptid), row.names(exprs(transcript_summaries)))
# save data into a file
write.table(transpos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/transcript_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)
rm(transcript_fdata, i, name, chr, min, max, row)

###################
### ALL SAMPLES ###
###################
# As we already generated the covariate and SNP file above, we only need 
# to generate the transcript expression file here
expr <- exprs(transcript_summaries)
expr[1:3, 1:6]
# check if transcript ID order is the same as in transpos
identical(unique(rownames(expr)), as.character(transpos$transcriptid))# Transcript expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for gene expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
colnames(expr)[1] <- "id"
identical(colnames(expr), colnames(snps))

# Save transcript expression data into a new file
fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_transcript_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get transcript expression data for SCZ (CASE) samples
transcript_expr_case_only <- subset(expr, select = case_samples)
transcript_expr_case_only[1:3, 1:6]
identical(colnames(transcript_expr_case_only), case_samples)
# save transcript expression into a file
fwrite(transcript_expr_case_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_transcript_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get transcript expression data for CTL samples
transcript_expr_ctl_only <- subset(expr, select = ctl_samples)
transcript_expr_ctl_only[1:3, 1:6]
identical(colnames(transcript_expr_ctl_only), ctl_samples)
# save transcript expression into a file
fwrite(transcript_expr_ctl_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_transcript_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)
rm(transcript_summaries, transpos, cols, expr, transcript_expr_case_only, transcript_expr_ctl_only)

################################
### 3. Exon eQTL Input Files ###
################################
# Load exon expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
dim(exon_summaries) # 242,443 features and 169 samples

# get phenotype data (including the SV1 variable)
pData <- pData(exon_summaries)
pData$ShortName <- factor(pData$ShortName, levels = pData$ShortName)

##########################
### GET EXON POSITIONS ###
##########################
# Get exon symbols, chr, start and end pos
exons_fdata <- exon_summaries@featureData@data[, c("exon_id", "location")]
exonspos <- data.frame()
for (i in 1:length(exons_fdata[, 1])) {
  name <- rownames(exons_fdata[i, ])
  chr <- strsplit(strsplit(as.character(exons_fdata[i, "location"]), ":")[[1]][1], "r")[[1]][2]
  min <- 0
  max <- 0
  # if several locations are known, get only the minium start and maximum end position of all locations
  if (grepl("\\|", exons_fdata[i, "location"])) {
    min <- min(sapply(strsplit(as.character(exons_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])))
    max <- max(sapply(strsplit(as.character(exons_fdata[i, "location"]), "\\|"), function(x) as.numeric(strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2])))
  } else {
    min <- strsplit(strsplit(as.character(exons_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][1]
    max <- strsplit(strsplit(as.character(exons_fdata[i, "location"]), ":")[[1]][2], "-")[[1]][2]
  }
  row <- cbind(name, chr, min, max)
  exonspos <- rbind(exonspos, row)
  if ((i %% 10000) == 0) {
    print(i)
  }
}
exonspos <- as.data.frame(exonspos)
colnames(exonspos) <- c("exonid", "chr", "start", "end")
exonspos$chr <- sapply(exonspos$chr, function(x) paste0("chr", x)) # add chr before chromosome number 
dim(exonspos)  # 242,443 and four columns
exonspos[1:3, ]

# check, if we now really have the same exons with locations
identical(as.vector(exonspos$exonid), rownames(exprs(exon_summaries)))
# save data into a file
write.table(exonspos, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/exon_locations.txt", col.names = T, row.names = F, sep = "\t", quote = F)
rm(exons_fdata, i, name, chr, min, max, row)

###################
### ALL SAMPLES ###
###################
# Get exon expression data
expr <- exprs(exon_summaries)
expr[1:3, 1:6]
# check if exon ID order is the same as in exonspos
identical(unique(rownames(expr)), as.character(exonspos$exonid)) # Exon expression: get overlapping samples to snps data
cols <- colnames(snps)  # get new colnames of snps
expr <- data.table(expr[, colnames(expr) %in% cols], keep.rownames = TRUE)  # get same samples for exon expression dataset
setcolorder(expr, c("rn", sort(colnames(expr)[-1])))
colnames(expr)[1] <- "id"
identical(colnames(expr), colnames(snps))

# Generate exon expression input files CASE & CTL
count <- 1
for (i in seq(1, nrow(expr), 50000)) {
  if (i < 200000) {
    chunk <- expr[i : (i + 50000 - 1), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_exon_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- expr[i : nrow(expr), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_exon_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}
# Save gene expression data into a new file
#fwrite(expr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/01_pmb_exon_expr_CASE_CTL.txt", row.names = F, col.names = T, sep = "\t", quote = F)

#########################
### ONLY CASE SAMPLES ###
#########################
# get exon expression data for SCZ (CASE) samples
exon_expr_case_only <- subset(expr, select = case_samples)
exon_expr_case_only[1:3, 1:6]
identical(colnames(exon_expr_case_only), case_samples)
# Generate exon expression input files CASE only
count <- 1
for (i in seq(1, nrow(exon_expr_case_only), 50000)) {
  if (i < 200000) {
    chunk <- exon_expr_case_only[i : (i + 50000 - 1), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_exon_expr_CASE_only_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- exon_expr_case_only[i : nrow(exon_expr_case_only), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_exon_expr_CASE_only_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}
#fwrite(exon_expr_case_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/02_pmb_exon_expr_CASE_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)

########################
### ONLY CTL SAMPLES ###
########################
# get exon expression data for CTL samples
exon_expr_ctl_only <- subset(expr, select = ctl_samples)
exon_expr_ctl_only[1:3, 1:6]
identical(colnames(exon_expr_ctl_only), ctl_samples)
# Generate exon expression input files CTL only
count <- 1
for (i in seq(1, nrow(exon_expr_ctl_only), 50000)) {
  if (i < 200000) {
    chunk <- exon_expr_ctl_only[i : (i + 50000 - 1), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_exon_expr_CTL_only_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- exon_expr_ctl_only[i : nrow(exon_expr_ctl_only), ]
    # Save exon expression data into a new chunk file
    oname <- paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_exon_expr_CTL_only_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}
#fwrite(exon_expr_ctl_only, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/02_matrixeqtl_input/03_pmb_exon_expr_CTL_only.txt", row.names = F, col.names = T, sep = "\t", quote = F)
