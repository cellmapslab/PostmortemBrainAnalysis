library(oligo)
library(rtracklayer)
library(ggplot2)
library(viridis)
library(data.table)
library(limma)
library(gridExtra)

# different colors for genes, transcripts and exons as well as control datasets
colors <- c("#FDAE6B", "#4292C6", "#A6D854")

### 1. Generate FKBP5 result table ###
# rows = patient/controls ID
# columns = FKBP5 gene expression + all FKBP5 transcript and exon expressions + all calculated PRS

# 1.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_gene_summaries_annotated_batch_corrected.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")
# pData
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_pData_with_surrogate_variables.RData")

# 1.2 Remove samples with missing values
# 1.2.1 remove all samples with C1 = NA from phenotype data and gene/transcript/exon expression data
index <- which(is.na(pData(geneSummaries)$C1))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 1.2.2 remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$Type_of_Death))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 1.2.3 check dim of geneSummaries (25 rows should be removed = 169 samples left)
dim(geneSummaries)
dim(transcriptSummaries)
dim(exonSummaries)
identical(rownames(pData(geneSummaries)), rownames(pData))

# 1.3 Get residuals for gene, transcript and exon expression
# 1.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
fit.res <- lmFit(geneSummaries, design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(geneSummaries))  # get residuals
#exprs(geneSummaries) <- res.gene
# 1.3.2 Get residuals for transcript expression
fit.res <- lmFit(transcriptSummaries, design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcriptSummaries))  # get residuals
#exprs(transcriptSummaries) <- res.transcript
# 1.3.3 Get residuals for exon expression (Remember: include SV1!!)
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exonSummaries, design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exonSummaries))  # get residuals
#exprs(exonSummaries) <- res.exon

# 1.4 Get FKBP5 gene, transcript and exon expressions
# 1.4.1 Get FKBP5 gene expression
fkbp5.id <- fData(geneSummaries)[grep("FKBP5", fData(geneSummaries)$gene_name), "gene_id"]
fkbp5.gene.expr <- as.data.frame(exprs(geneSummaries)[fkbp5.id, ])
colnames(fkbp5.gene.expr) <- fkbp5.id

# 1.4.2 Get FKBP5 transcript expression
fkbp5.id <- fData(transcriptSummaries)[grep("FKBP5", fData(transcriptSummaries)$gene_name), "transcript_id"]
fkbp5.transcript.expr <- t(as.data.frame(exprs(transcriptSummaries)[fkbp5.id, ]))

# 1.4.3 Get FKBP5 exon expression
fkbp5.id <- fData(exonSummaries)[grep("FKBP5", fData(exonSummaries)$gene_name), "exon_id"]
fkbp5.exon.expr <- t(as.data.frame(exprs(exonSummaries)[fkbp5.id, ]))

# 1.5 Load PRS
# 1.5.1 load SCZ2 PRS
prs.scz2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/SCZ2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_scz2_filtered_PRS_LDpred-inf.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# 1.5.2 load MDD2 PRS
prs.mdd2 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/MDD2/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_mdd2_filtered_PRS_LDpred-inf.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# 1.5.3 load BIP-2018 PRS
prs.bip <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/BIP-2018/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_bip_filtered_PRS_LDpred-inf.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))
# 1.5.4 load Cross-disorder PRS
prs.cross.disorder <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/LDpred/Cross_Disorder/03_validate/03_updateMap_all_chr_maf0.01_hwe_IDs_changed_noDup_cd_filtered_PRS_LDpred-inf.txt", sep = ",", header = T, stringsAsFactors = F, select = c("IID", "PRS"))

# 1.6 Merge PRS to pData (to get the patient/control IDs)
pData <- pData(geneSummaries)[, c("IID", "ShortName")]
pData$PRS.SCZ <- unlist(prs.scz2[prs.scz2$IID %in% pData$IID, "PRS"])
pData$PRS.MDD <- unlist(prs.mdd2[prs.mdd2$IID %in% pData$IID, "PRS"])
pData$PRS.BIP <- unlist(prs.bip[prs.bip$IID %in% pData$IID, "PRS"])
pData$PRS.CD <- unlist(prs.cross.disorder[prs.cross.disorder$IID %in% pData$IID, "PRS"])
pData$ShortName <- NULL

# 1.7 Merge datasets
fkbp5 <- fkbp5.gene.expr
fkbp5 <- merge(fkbp5, fkbp5.transcript.expr, by = 0)
fkbp5 <- merge(fkbp5, fkbp5.exon.expr, by.x = "Row.names", by.y = "row.names")
fkbp5 <- merge(fkbp5, pData, by.x = "Row.names", by.y = "row.names")
colnames(fkbp5)[1] <- "SampleID"
dim(fkbp5)  # should be 169 rows (samples)

# 1.8 Save dataset
write.table(fkbp5, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/results/fkbp5_expressions_and_PRS.tsv", row.names = FALSE, sep = "\t")
rm(prs.scz2, prs.mdd2, prs.bip, prs.cross.disorder, fkbp5.id, index)

###############################
### 2. Generate FKBP5 plots ###
###############################
# Plot of gene/transcript/exon expression
# 2.1 remove PRS from the data and form data into right order
fkbp5 <- fkbp5[, -c((length(fkbp5)-4):length(fkbp5))]
# exon 8: ENSE00000747303 is filtered out because the crosshyb_type equals 3!
fkbp5 <- fkbp5[, c("SampleID", "ENSG00000096060.10", "ENST00000357266.4|ENST00000539068.1", "ENST00000536438.1", "ENST00000542713.1", "ENST00000540787.1", "ENSE00001519071.2|ENSE00001354888.3|ENSE00002312621.1|ENSE00002286243.1", "ENSE00001370159.1", "ENSE00000747346.1", "ENSE00001391213.1", "ENSE00000747342.1", "ENSE00003659116.1|ENSE00003674427.1", "ENSE00000747304.1", "ENSE00000747302.1", "ENSE00000747301.1", "ENSE00001193215.7|ENSE00001354861.2")]
head(fkbp5)

# 2.2 Change fkbp5 data into ggplot long format
colnames(fkbp5) <- c("SampleID", "FKBP5", "FKBP5-001/FKBP5-202", "FKBP5-201", "FKBP5-203", "FKBP5-204", "ENSE00001354888", "ENSE00001370159", "ENSE00000747346", "ENSE00001391213", "ENSE00000747342", "ENSE00003659116", "ENSE00000747304", "ENSE00000747302", "ENSE00000747301", "ENSE00001354861")
fkbp5.plot.data <- melt(fkbp5, id.vars = "SampleID", value.name = "Expression", variable.name = "EnsemblID" )
fkbp5.plot.data$Type <- "Exons"
fkbp5.plot.data[as.character(fkbp5.plot.data$EnsemblID) == "FKBP5", "Type"] <- "Gene"
fkbp5.plot.data[startsWith(as.character(fkbp5.plot.data$EnsemblID), "FKBP5-"), "Type"] <- "Transcripts"
fkbp5.plot.data$Type <- as.factor(fkbp5.plot.data$Type)
# relevel to gene, transcript, exon order
fkbp5.plot.data$Type <- factor(fkbp5.plot.data$Type, levels = c("Gene", "Transcripts", "Exons"))
# add case/control information
fkbp5.plot.data$SampleType <- "Cases"
controls <-pData(geneSummaries)[pData(geneSummaries)$Dx == "CTL", "ShortName"]
fkbp5.plot.data[fkbp5.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(fkbp5.plot.data)

# 2.3 Generate FKBP5 plot
cols <- c(colors[1], "gray")
#pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/figures/boxplot_FKBP5_gene_vs_transcripts_vs_exons.pdf", width = 15)
sp <- ggplot(fkbp5.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) + #, labels = c("gene", "controls", "transcript 1/3", "transcript 2", "transcript 4", "transcript 5", "exon 1", "exon 2", "exon 3", "exon 4", "exon 5", "exon 6", "exon 7", "exon 9", "exon 10", "exon 11"), name = "") + 
  #scale_x_discrete(labels = c("FKBP5", "FKBP5-001/FKBP5-202", "FKBP5-201", "FKBP5-203", "FKBP5-204", "ENSE00001354888", "ENSE00001370159", "ENSE00000747346", "ENSE00001391213", "ENSE00000747342", "ENSE00003659116", "ENSE00000747304", "ENSE00000747303", "ENSE00000747302", "ENSE00000747301", "ENSE00001354861")) + 
  geom_hline(yintercept = median(fkbp5.plot.data[fkbp5.plot.data$EnsemblID == "FKBP5", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
#dev.off()

# boxplot of cases vs. controls
# remove all samples with C1 = NA from phenotype data and gene/exon expression data
index <- which(is.na(pData(geneSummaries)$C1))
geneSummaries <- geneSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]
# remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$Type_of_Death))
geneSummaries <- geneSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]
# check dim of geneSummaries (25 rows should be removed)
dim(geneSummaries)
pData <- pData(geneSummaries)

# add original Dx as Diagnosis column to pData
pData["Diagnosis"] <- pData$Dx
# change BPD and MDD diagnosis (Dx) to SCZ
pData$Dx <- as.character(pData$Dx) # change factor column to character
pData[pData$Dx == "BPD" | pData$Dx == "MDD", "Dx"] <- "SCZ" # change BPD and MDD to SCZ
pData[, c("Sex", "Dx", "Suicide")] <- lapply(pData[, c("Sex", "Dx", "Suicide")], factor)
cases <- pData[pData$Dx == "SCZ", "ShortName"]
ctrl <- pData[pData$Dx == "CTL", "ShortName"]
fkbp.cases <- as.data.frame(fkbp5.transcript.expr[which(names(fkbp5.transcript.expr) %in% cases)])
fkbp.ctrls <- as.data.frame(fkbp5.transcript.expr[which(names(fkbp5.transcript.expr) %in% ctrl)])
fkbp.cases <- cbind(fkbp.cases, rep("Cases", 107))
colnames(fkbp.cases) <- c("Expression", "Dx")
fkbp.ctrls <- cbind(fkbp.ctrls, rep("Controls", 62))
colnames(fkbp.ctrls) <- c("Expression", "Dx")
fkbp <- rbind(fkbp.cases, fkbp.ctrls)
ggplot(fkbp, aes(x = Dx, y = Expression, fill = Dx)) + 
  geom_boxplot() +
  theme_grey(base_size = 18)


###############################
### 3. Generate EPM2A plots ###
###############################
# 3.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_gene_summaries_annotated_batch_corrected.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")
# pData
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_pData_with_surrogate_variables.RData")

# 3.2 Remove samples with missing values
# 3.2.1 remove all samples with C1 = NA from phenotype data and gene/transcript/exon expression data
index <- which(is.na(pData(geneSummaries)$C1))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 3.2.2 remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$Type_of_Death))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 3.2.3 check dim of geneSummaries (25 rows should be removed = 169 samples left)
dim(geneSummaries)
dim(transcriptSummaries)
dim(exonSummaries)
identical(rownames(pData(geneSummaries)), rownames(pData))

# 3.3 Get residuals for gene, transcript and exon expression
# 3.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
fit.res <- lmFit(geneSummaries, design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(geneSummaries))  # get residuals
exprs(geneSummaries) <- res.gene
# 3.3.2 Get residuals for transcript expression
fit.res <- lmFit(transcriptSummaries, design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcriptSummaries))  # get residuals
exprs(transcriptSummaries) <- res.transcript
# 3.3.3 Get residuals for exon expression (Remember: include SV1!!)
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exonSummaries, design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exonSummaries))  # get residuals
exprs(exonSummaries) <- res.exon

# 3.4 Get EPM2A gene, transcript and exon expressions
# 3.4.1 Get EPM2A gene expression
epm2a.id <- fData(geneSummaries)[fData(geneSummaries)$gene_name == "EPM2A", "gene_id"]
epm2a.gene.expr <- as.data.frame(exprs(geneSummaries)[epm2a.id, ])
colnames(epm2a.gene.expr) <- epm2a.id

# 3.4.2 Get EPM2A transcript expression
epm2a.id <- fData(transcriptSummaries)[fData(transcriptSummaries)$gene_name == "EPM2A", "transcript_id"]
epm2a.transcript.expr <- t(as.data.frame(exprs(transcriptSummaries)[epm2a.id, ]))

# 3.4.3 Get EPM2A exon expression
epm2a.id <- fData(exonSummaries)[fData(exonSummaries)$gene_name == "EPM2A", "exon_id"]
epm2a.exon.expr <- t(as.data.frame(exprs(exonSummaries)[epm2a.id, ]))

# 3.5 Merge datasets
epm2a <- epm2a.gene.expr
epm2a <- merge(epm2a, epm2a.transcript.expr, by = 0)
epm2a <- merge(epm2a, epm2a.exon.expr, by.x = "Row.names", by.y = "row.names")
epm2a <- merge(epm2a, pData, by.x = "Row.names", by.y = "row.names")
colnames(epm2a)[1] <- "SampleID"
dim(epm2a)  # should be 169 rows (samples)

# 3.6 Form data into right order
epm2a <- epm2a[, c("SampleID", "ENSG00000112425.9", "ENST00000450221.1", "ENST00000367519.3", "ENST00000461700.1", "ENST00000435470.1", "ENST00000496228.1", "ENST00000489412.1", "ENSE00001898858.1", "ENSE00001952999.1", "ENSE00001795677.1", "ENSE00001886960.1", "ENSE00001511364.2", "ENSE00001825781.1", "ENSE00001955317.1", "ENSE00001849769.1", "ENSE00001723502.1", "ENSE00001658660.1", "ENSE00001946407.1", "ENSE00001867523.1", "ENSE00001906493.1", "ENSE00001782653.1")]
colnames(epm2a) <- c("SampleID", "EPM2A", "EPM2A-001", "EPM2A-002", "EPM2A-003", "EPM2A-004", "EPM2A-005", "EPM2A-006", "ENSE00001898858.1", "ENSE00001952999.1", "ENSE00001795677.1", "ENSE00001886960.1", "ENSE00001511364.2", "ENSE00001825781.1", "ENSE00001955317.1", "ENSE00001849769.1", "ENSE00001723502.1", "ENSE00001658660.1", "ENSE00001946407.1", "ENSE00001867523.1", "ENSE00001906493.1", "ENSE00001782653.1")
head(epm2a)

# 3.7 Change EPM2A data into ggplot long format
epm2a.plot.data <- melt(epm2a, id.vars = "SampleID", value.name = "Expression", variable.name = "EnsemblID" )
epm2a.plot.data$Type <- "Exons"
epm2a.plot.data[as.character(epm2a.plot.data$EnsemblID) == "EPM2A", "Type"] <- "Gene"
epm2a.plot.data[startsWith(as.character(epm2a.plot.data$EnsemblID), "EPM2A-"), "Type"] <- "Transcripts"
epm2a.plot.data$Type <- as.factor(epm2a.plot.data$Type)
# relevel to gene, transcript, exon order
epm2a.plot.data$Type <- factor(epm2a.plot.data$Type, levels = c("Gene", "Transcripts", "Exons"))
# add case/control information
epm2a.plot.data$SampleType <- "Cases"
controls <-pData(geneSummaries)[pData(geneSummaries)$Dx == "CTL", "ShortName"]
epm2a.plot.data[epm2a.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(epm2a.plot.data)

# 3.8 Generate EPM2A plot
cols <- c(colors[1], "gray")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/figures/boxplot_EPM2A_gene_vs_transcripts_vs_exons_residuals.pdf", width = 15)
sp <- ggplot(epm2a.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) + #, labels = c("gene", "transcript 1", "transcript 2", "transcript 3", "transcript 4", "transcript 5", "transcript 6", "exon 1", "exon 2", "exon 3", "exon 4", "exon 5", "exon 6", "exon 7", "exon 8", "exon 9", "exon 10", "exon 11", "exon 12", "exon 13", "exon 14"), name = "") + 
  #scale_x_discrete(labels = c("EPM2A", "EPM2A-001", "EPM2A-002", "EPM2A-003", "EPM2A-004", "EPM2A-005", "EPM2A-006", "ENSG00000112425.9", "ENST00000450221.1", "ENST00000367519.3", "ENST00000461700.1", "ENST00000435470.1", "ENST00000496228.1", "ENST00000489412.1", "ENSE00001898858.1", "ENSE00001952999.1", "ENSE00001795677.1", "ENSE00001886960.1", "ENSE00001511364.2", "ENSE00001825781.1", "ENSE00001955317.1", "ENSE00001849769.1", "ENSE00001723502.1", "ENSE00001658660.1", "ENSE00001946407.1", "ENSE00001867523.1", "ENSE00001906493.1", "ENSE00001782653.1")) + 
  geom_hline(yintercept = median(epm2a.plot.data[epm2a.plot.data$EnsemblID == "EPM2A", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
dev.off()

#################################
### 3. Generate SERTAD4 plots ###
#################################
# 3.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_gene_summaries_annotated_batch_corrected.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")
# pData
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_pData_with_surrogate_variables.RData")

# 3.2 Remove samples with missing values
# 3.2.1 remove all samples with C1 = NA from phenotype data and gene/transcript/exon expression data
index <- which(is.na(pData(geneSummaries)$C1))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 3.2.2 remove all samples with Type_of_Death = NA from phenotype data and gene expression data
index <- which(is.na(pData(geneSummaries)$Type_of_Death))
geneSummaries <- geneSummaries[, -c(index)]
transcriptSummaries <- transcriptSummaries[, -c(index)]
exonSummaries <- exonSummaries[, -c(index)]

# 3.2.3 check dim of geneSummaries (25 rows should be removed = 169 samples left)
dim(geneSummaries)
dim(transcriptSummaries)
dim(exonSummaries)
identical(rownames(pData(geneSummaries)), rownames(pData))

# 3.3 Get residuals for gene, transcript and exon expression
# 3.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
fit.res <- lmFit(geneSummaries, design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(geneSummaries))  # get residuals
exprs(geneSummaries) <- res.gene
# 3.3.2 Get residuals for transcript expression
fit.res <- lmFit(transcriptSummaries, design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcriptSummaries))  # get residuals
exprs(transcriptSummaries) <- res.transcript
# 3.3.3 Get residuals for exon expression (Remember: include SV1!!)
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exonSummaries, design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exonSummaries))  # get residuals
exprs(exonSummaries) <- res.exon

# 3.4 Get SERTAD4 gene, transcript and exon expressions
# 3.4.1 Get SERTAD4 gene expression
sertad4.id <- fData(geneSummaries)[fData(geneSummaries)$gene_name == "SERTAD4", "gene_id"]
sertad4.gene.expr <- as.data.frame(exprs(geneSummaries)[sertad4.id, ])
colnames(sertad4.gene.expr) <- sertad4.id

# 3.4.2 Get SERTAD4 transcript expression
sertad4.id <- fData(transcriptSummaries)[fData(transcriptSummaries)$gene_name == "SERTAD4", "transcript_id"]
sertad4.transcript.expr <- t(as.data.frame(exprs(transcriptSummaries)[sertad4.id, ]))

# 3.4.3 Get SERTAD4 exon expression
sertad4.id <- fData(exonSummaries)[fData(exonSummaries)$gene_name == "SERTAD4", "exon_id"]
sertad4.exon.expr <- t(as.data.frame(exprs(exonSummaries)[sertad4.id, ]))

# 3.5 Merge datasets
sertad4 <- sertad4.gene.expr
sertad4 <- merge(sertad4, sertad4.transcript.expr, by = 0)
sertad4 <- merge(sertad4, sertad4.exon.expr, by.x = "Row.names", by.y = "row.names")
sertad4 <- merge(sertad4, pData, by.x = "Row.names", by.y = "row.names")
colnames(sertad4)[1] <- "SampleID"
dim(sertad4)  # should be 169 rows (samples)

# 3.6 Form data into right order
sertad4 <- sertad4[, c("SampleID", "ENSG00000082497.7", "ENST00000490620.1", "ENST00000367012.3", "ENSE00001443240.3", "ENSE00001443238.3", "ENSE00001935710.1", "ENSE00001917877.1")]
colnames(sertad4) <- c("SampleID", "SERTAD4", "SERTAD4-001", "SERTAD4-003", "ENSE00001443240.3", "ENSE00001443238.3", "ENSE00001935710.1", "ENSE00001917877.1")
head(sertad4)

# 3.7 Change SERTAD4 data into ggplot long format
sertad4.plot.data <- melt(sertad4, id.vars = "SampleID", value.name = "Expression", variable.name = "EnsemblID" )
sertad4.plot.data$Type <- "Exons"
sertad4.plot.data[as.character(sertad4.plot.data$EnsemblID) == "SERTAD4", "Type"] <- "Gene"
sertad4.plot.data[startsWith(as.character(sertad4.plot.data$EnsemblID), "SERTAD4-"), "Type"] <- "Transcripts"
sertad4.plot.data$Type <- as.factor(sertad4.plot.data$Type)
# relevel to gene, transcript, exon order
sertad4.plot.data$Type <- factor(sertad4.plot.data$Type, levels = c("Gene", "Transcripts", "Exons"))
# add case/control information
sertad4.plot.data$SampleType <- "Cases"
controls <-pData(geneSummaries)[pData(geneSummaries)$Dx == "CTL", "ShortName"]
sertad4.plot.data[sertad4.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(sertad4.plot.data)

# 3.8 Generate SERTAD4 plot
cols <- c(colors[1], "gray")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/figures/boxplot_SERTAD4_gene_vs_transcripts_vs_exons_residuals.pdf", width = 15)
sp <- ggplot(sertad4.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) + #, labels = c("gene", "transcript 1", "transcript 3", "exon 1", "exon 2", "exon 3", "exon 4"), name = "") + 
  #scale_x_discrete(labels = c("SERTAD4", "SERTAD4-001", "SERTAD4-003", "ENSE00001443240.3", "ENSE00001443238.3", "ENSE00001935710.1", "ENSE00001917877.1")) + 
  geom_hline(yintercept = median(sertad4.plot.data[sertad4.plot.data$EnsemblID == "SERTAD4", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
dev.off()
