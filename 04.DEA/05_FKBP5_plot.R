library(oligo)
library(rtracklayer)
library(ggplot2)
library(viridis)
library(data.table)
library(limma)
library(griDEAtra)

# different colors for genes, transcripts and exons as well as control datasets
colors <- c("#FDAE6B", "#4292C6", "#A6D854")

### 1. Generate FKBP5 result table ###
# rows = patient/controls ID
# columns = FKBP5 gene expression + all FKBP5 transcript and exon expressions + all calculated PRS

# 1.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_gene_summaries_annotated_batch_corrected.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")
# pData
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_pData_with_surrogate_variables.RData")

# 1.2 Remove samples with missing values
# 1.2.1 remove all samples with C1 = NA from phenotype data and gene/transcript/exon expression data
inDEA <- which(is.na(pData(gene_summaries)$C1))
gene_summaries <- gene_summaries[, -c(inDEA)]
transcript_summaries <- transcript_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]

# 1.2.2 remove all samples with Type_of_Death = NA from phenotype data and gene expression data
inDEA <- which(is.na(pData(gene_summaries)$Type_of_Death))
gene_summaries <- gene_summaries[, -c(inDEA)]
transcript_summaries <- transcript_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]

# 1.2.3 check dim of gene_summaries (25 rows should be removed = 169 samples left)
dim(gene_summaries)
dim(transcript_summaries)
dim(exon_summaries)
identical(rownames(pData(gene_summaries)), rownames(pData))

# 1.3 Get residuals for gene, transcript and exon expression
# 1.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
fit.res <- lmFit(gene_summaries, design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(gene_summaries))  # get residuals
#exprs(gene_summaries) <- res.gene
# 1.3.2 Get residuals for transcript expression
fit.res <- lmFit(transcript_summaries, design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcript_summaries))  # get residuals
#exprs(transcript_summaries) <- res.transcript
# 1.3.3 Get residuals for exon expression (Remember: include SV1!!)
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exon_summaries, design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exon_summaries))  # get residuals
#exprs(exon_summaries) <- res.exon

# 1.4 Get FKBP5 gene, transcript and exon expressions
# 1.4.1 Get FKBP5 gene expression
fkbp5.id <- fData(gene_summaries)[grep("FKBP5", fData(gene_summaries)$gene_name), "gene_id"]
fkbp5.gene.expr <- as.data.frame(exprs(gene_summaries)[fkbp5.id, ])
colnames(fkbp5.gene.expr) <- fkbp5.id

# 1.4.2 Get FKBP5 transcript expression
fkbp5.id <- fData(transcript_summaries)[grep("FKBP5", fData(transcript_summaries)$gene_name), "transcript_id"]
fkbp5.transcript.expr <- t(as.data.frame(exprs(transcript_summaries)[fkbp5.id, ]))

# 1.4.3 Get FKBP5 exon expression
fkbp5.id <- fData(exon_summaries)[grep("FKBP5", fData(exon_summaries)$gene_name), "exon_id"]
fkbp5.exon.expr <- t(as.data.frame(exprs(exon_summaries)[fkbp5.id, ]))

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
pData <- pData(gene_summaries)[, c("IID", "ShortName")]
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
write.table(fkbp5, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/results/fkbp5_expressions_and_PRS.tsv", row.names = FALSE, sep = "\t")
rm(prs.scz2, prs.mdd2, prs.bip, prs.cross.disorder, fkbp5.id, inDEA)

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
controls <-pData(gene_summaries)[pData(gene_summaries)$Dx == "CTL", "ShortName"]
fkbp5.plot.data[fkbp5.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(fkbp5.plot.data)

# 2.3 Generate FKBP5 plot
cols <- c(colors[1], "gray")
#pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/figures/boxplot_FKBP5_gene_vs_transcripts_vs_exons.pdf", width = 15)
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
inDEA <- which(is.na(pData(gene_summaries)$C1))
gene_summaries <- gene_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]
# remove all samples with Type_of_Death = NA from phenotype data and gene expression data
inDEA <- which(is.na(pData(gene_summaries)$Type_of_Death))
gene_summaries <- gene_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]
# check dim of gene_summaries (25 rows should be removed)
dim(gene_summaries)
pData <- pData(gene_summaries)

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


#############################
### 3. Generate DPP4 plot ###
#############################
# DPP4 is also a hit in the 108 loci (CMC) Nature paper!!!
# 3.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
# pData
pData <- pData(gene_summaries)

# check dim of all Summaries (169 samples each)
dim(gene_summaries)
dim(transcript_summaries)
dim(exon_summaries)

# 3.3 Get residuals for gene, transcript and exon expression
# 3.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(~ Dx + Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exprs(transcript_summaries), design.res)  # fit model
res.table <- topTable(eBayes(fit.res), coef = "DxCASE", adjust = "BH", number = Inf)

design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exprs(gene_summaries), design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(gene_summaries))  # get residuals
exprs(gene_summaries) <- res.gene
# 3.3.2 Get residuals for transcript expression
fit.res <- lmFit(exprs(transcript_summaries), design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcript_summaries))  # get residuals
exprs(transcript_summaries) <- res.transcript
# 3.3.3 Get residuals for exon expression (Remember: include SV1!!)
fit.res <- lmFit(exprs(exon_summaries), design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exon_summaries))  # get residuals
exprs(exon_summaries) <- res.exon

# 3.4 Get DPP4 gene, transcript and exon expressions
# 3.4.1 Get DPP4 gene expression
dpp4.id <- fData(gene_summaries)[fData(gene_summaries)$gene_name == "DPP4", "gene_id"]
dpp4.gene.expr <- as.data.frame(exprs(gene_summaries)[dpp4.id, ])
colnames(dpp4.gene.expr) <- "DPP4"

# 3.4.2 Get DPP4 transcript expression
dpp4.id <- fData(transcript_summaries)[grep("DPP4", fData(transcript_summaries)$gene_name), "transcript_id"]
dpp4.transcript.expr <- t(as.data.frame(exprs(transcript_summaries)[dpp4.id, ]))
dpp4.transcript.expr <- dpp4.transcript.expr[, c("ENST00000360534.3", "ENST00000434918.2", "ENST00000497461.1", "ENST00000416189.1", "ENST00000461836.1", "ENST00000490286.1", "ENST00000413651.1", "ENST00000468903.1", "ENST00000491591.1", "ENST00000494507.1")]
colnames(dpp4.transcript.expr) <- c("DPP4-001", "DPP4-002", "DPP4-003", "DPP4-004", "DPP4-005", "DPP4-006", "DPP4-007", "DPP4-008", "DPP4-010", "DPP4-011")

# 3.4.3 Get DPP4 exon expression
dpp4.id <- fData(exon_summaries)[grep("DPP4", fData(exon_summaries)$gene_name), "exon_id"]
dpp4.exon.expr <- t(as.data.frame(exprs(exon_summaries)[dpp4.id, ]))
# all 7 exons in best scoring transcript of DPP4 (DPP4-008): 
dpp4_08 <- c("ENSE00001845011", "ENSE00003482771", "ENSE00003513529", "ENSE00003668522", "ENSE00003517002", "ENSE00003590666", "ENSE00001879445")
exons <- sapply(dpp4_08, function(x) dpp4.exon.expr[, grep(x, colnames(dpp4.exon.expr))])
dpp4_08_exons <- rbind(exons[[2]], exons[[3]], exons[[5]], exons[[7]])
rownames(dpp4_08_exons) <- names(exons)[c(2, 3, 5, 7)]
dpp4.exon.expr <- t(dpp4_08_exons)

# 3.5 Merge datasets
dpp4 <- dpp4.gene.expr
dpp4 <- merge(dpp4, dpp4.transcript.expr, by = 0)
dpp4 <- merge(dpp4, dpp4.exon.expr, by.x = "Row.names", by.y = "row.names")
#dpp4 <- merge(dpp4, pData, by.x = "Row.names", by.y = "row.names")
colnames(dpp4)[1] <- "SampleID"
dim(dpp4)  # should be 169 rows (samples)
head(dpp4)

# 3.7 Change dpp4 data into ggplot long format
dpp4.plot.data <- reshape2::melt(dpp4, id.vars = "SampleID", value.name = "Expression", variable.name = "EnsemblID" )
dpp4.plot.data$Type <- "Exons"
dpp4.plot.data[as.character(dpp4.plot.data$EnsemblID) == "DPP4", "Type"] <- "Gene"
dpp4.plot.data[startsWith(as.character(dpp4.plot.data$EnsemblID), "DPP4-"), "Type"] <- "Transcripts"
dpp4.plot.data$Type <- as.factor(dpp4.plot.data$Type)
# relevel to gene, transcript, exon order
dpp4.plot.data$Type <- factor(dpp4.plot.data$Type, levels = c("Gene", "Transcripts", "Exons"))
# add case/control information
dpp4.plot.data$SampleType <- "Cases"
controls <-pData(gene_summaries)[pData(gene_summaries)$Dx == "CTL", "ShortName"]
dpp4.plot.data[dpp4.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(dpp4.plot.data)

# 3.8 Generate dpp4 plot
cols <- c(colors[1], "gray", colors[2], "gray", colors[3], "gray")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/figures/boxplot_dpp4_gene_vs_transcripts_vs_exons_residuals.pdf", width = 15)
sp <- ggplot(dpp4.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = median(dpp4.plot.data[dpp4.plot.data$EnsemblID == "DPP4", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
dev.off()

############################
### 3. Generate EYS plot ###
############################
# EYS is a GWAS hit in two different schizophrenia paper!!
# 3.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_transcript_summaries_annotated.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
# pData
pData <- pData(gene_summaries)

# check dim of all Summaries (169 samples each)
dim(gene_summaries)
dim(transcript_summaries)
dim(exon_summaries)

# 3.3 Get residuals for gene, transcript and exon expression
# 3.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exprs(gene_summaries), design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(gene_summaries))  # get residuals
exprs(gene_summaries) <- res.gene
# 3.3.2 Get residuals for transcript expression
fit.res <- lmFit(exprs(transcript_summaries), design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcript_summaries))  # get residuals
exprs(transcript_summaries) <- res.transcript
# 3.3.3 Get residuals for exon expression (Remember: include SV1!!)
fit.res <- lmFit(exprs(exon_summaries), design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exon_summaries))  # get residuals
exprs(exon_summaries) <- res.exon

# 3.4 Get FNDC3A gene, transcript and exon expressions
# 3.4.1 Get FNDC3A gene expression
fndc3a.id <- fData(gene_summaries)[fData(gene_summaries)$gene_name == "FNDC3A", "gene_id"]
fndc3a.gene.expr <- as.data.frame(exprs(gene_summaries)[fndc3a.id, ])
colnames(fndc3a.gene.expr) <- fndc3a.id

# 3.4.2 Get FNDC3A transcript expression
fndc3a.id <- fData(transcript_summaries)[grep("FNDC3A", fData(transcript_summaries)$gene_name), "transcript_id"]
fndc3a.transcript.expr <- t(as.data.frame(exprs(transcript_summaries)[fndc3a.id, ]))

# 3.4.3 Get FNDC3A exon expression
fndc3a.id <- fData(exon_summaries)[grep("FNDC3A", fData(exon_summaries)$gene_name), "exon_id"]
fndc3a.exon.expr <- t(as.data.frame(exprs(exon_summaries)[fndc3a.id, ]))

# 3.5 Merge datasets
fndc3a <- fndc3a.gene.expr
fndc3a <- merge(fndc3a, fndc3a.transcript.expr, by = 0)
fndc3a <- merge(fndc3a, fndc3a.exon.expr, by.x = "Row.names", by.y = "row.names")
#fndc3a <- merge(fndc3a, pData, by.x = "Row.names", by.y = "row.names")
colnames(fndc3a)[1] <- "SampleID"
dim(fndc3a)  # should be 169 rows (samples)

# 3.6 Form data into right order
fndc3a <- fndc3a[, c("SampleID", "ENSG00000102531.12", "ENST00000378383.1", "ENST00000398316.3", "ENST00000492622.2", "ENST00000497644.1", "ENSE00000682482.1", "ENSE00001477281.1", "ENSE00001342704.1", "ENSE00001532665.1", "ENSE00001129140.1", "ENSE00001335285.2", "ENSE00000682481.1", "ENSE00001477293.1", "ENSE00000939645.1", "ENSE00001129152.1", "ENSE00001422888.3", "ENSE00001922094.2", "ENSE00001129131.1", "ENSE00001129122.1", "ENSE00001342700.1")]
colnames(fndc3a) <- c("SampleID", "FNDC3A", "FNDC3A-001", "FNDC3A-002", "FNDC3A-003", "FNDC3A-006", "ENSE00000682482.1", "ENSE00001477281.1", "ENSE00001342704.1", "ENSE00001532665.1", "ENSE00001129140.1", "ENSE00001335285.2", "ENSE00000682481.1", "ENSE00001477293.1", "ENSE00000939645.1", "ENSE00001129152.1", "ENSE00001422888.3", "ENSE00001922094.2", "ENSE00001129131.1", "ENSE00001129122.1", "ENSE00001342700.1")
head(fndc3a)

# 3.7 Change FNDC3A data into ggplot long format
fndc3a.plot.data <- reshape2::melt(fndc3a, id.vars = "SampleID", value.name = "Expression", variable.name = "EnsemblID" )
fndc3a.plot.data$Type <- "Exons"
fndc3a.plot.data[as.character(fndc3a.plot.data$EnsemblID) == "FNDC3A", "Type"] <- "Gene"
fndc3a.plot.data[startsWith(as.character(fndc3a.plot.data$EnsemblID), "FNDC3A-"), "Type"] <- "Transcripts"
fndc3a.plot.data$Type <- as.factor(fndc3a.plot.data$Type)
# relevel to gene, transcript, exon order
fndc3a.plot.data$Type <- factor(fndc3a.plot.data$Type, levels = c("Gene", "Transcripts", "Exons"))
# add case/control information
fndc3a.plot.data$SampleType <- "Cases"
controls <-pData(gene_summaries)[pData(gene_summaries)$Dx == "CTL", "ShortName"]
fndc3a.plot.data[fndc3a.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(fndc3a.plot.data)

# 3.8 Generate FNDC3A plot
cols <- c(colors[1], "gray", colors[2], "gray", colors[3], "gray")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/figures/1.1_boxplot_FNDC3A_gene_vs_transcripts_vs_exons_residuals.pdf", width = 15)
sp <- ggplot(fndc3a.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) +
  geom_hline(yintercept = median(fndc3a.plot.data[fndc3a.plot.data$EnsemblID == "FNDC3A", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
dev.off()


#################################
### 3. Generate SERTAD4 plots ###
#################################
# 3.1 Load data
# Gene expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_gene_summaries_annotated_batch_corrected.RData")
# Transcript expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_transcript_summaries_annotated_batch_corrected.RData")
# Exon expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_exon_summaries_annotated_batch_corrected.RData")
# pData
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/RData_Files/7_pData_with_surrogate_variables.RData")

# 3.2 Remove samples with missing values
# 3.2.1 remove all samples with C1 = NA from phenotype data and gene/transcript/exon expression data
inDEA <- which(is.na(pData(gene_summaries)$C1))
gene_summaries <- gene_summaries[, -c(inDEA)]
transcript_summaries <- transcript_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]

# 3.2.2 remove all samples with Type_of_Death = NA from phenotype data and gene expression data
inDEA <- which(is.na(pData(gene_summaries)$Type_of_Death))
gene_summaries <- gene_summaries[, -c(inDEA)]
transcript_summaries <- transcript_summaries[, -c(inDEA)]
exon_summaries <- exon_summaries[, -c(inDEA)]

# 3.2.3 check dim of gene_summaries (25 rows should be removed = 169 samples left)
dim(gene_summaries)
dim(transcript_summaries)
dim(exon_summaries)
identical(rownames(pData(gene_summaries)), rownames(pData))

# 3.3 Get residuals for gene, transcript and exon expression
# 3.3.1 Get residuals for gene expression
# generate design matrix
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4, pData)
fit.res <- lmFit(gene_summaries, design.res)  # fit model
res.gene <- residuals.MArrayLM(fit.res, exprs(gene_summaries))  # get residuals
exprs(gene_summaries) <- res.gene
# 3.3.2 Get residuals for transcript expression
fit.res <- lmFit(transcript_summaries, design.res)  # fit model
res.transcript <- residuals.MArrayLM(fit.res, exprs(transcript_summaries))  # get residuals
exprs(transcript_summaries) <- res.transcript
# 3.3.3 Get residuals for exon expression (Remember: include SV1!!)
design.res <- model.matrix(Dx ~ Age + Sex + pH + PMI + RIN + RIN2 + Suicide + Type_of_Death + C1 + C2 + C3 + C4 + SV1, pData)
fit.res <- lmFit(exon_summaries, design.res)  # fit model
res.exon <- residuals.MArrayLM(fit.res, exprs(exon_summaries))  # get residuals
exprs(exon_summaries) <- res.exon

# 3.4 Get SERTAD4 gene, transcript and exon expressions
# 3.4.1 Get SERTAD4 gene expression
sertad4.id <- fData(gene_summaries)[fData(gene_summaries)$gene_name == "SERTAD4", "gene_id"]
sertad4.gene.expr <- as.data.frame(exprs(gene_summaries)[sertad4.id, ])
colnames(sertad4.gene.expr) <- sertad4.id

# 3.4.2 Get SERTAD4 transcript expression
sertad4.id <- fData(transcript_summaries)[fData(transcript_summaries)$gene_name == "SERTAD4", "transcript_id"]
sertad4.transcript.expr <- t(as.data.frame(exprs(transcript_summaries)[sertad4.id, ]))

# 3.4.3 Get SERTAD4 exon expression
sertad4.id <- fData(exon_summaries)[fData(exon_summaries)$gene_name == "SERTAD4", "exon_id"]
sertad4.exon.expr <- t(as.data.frame(exprs(exon_summaries)[sertad4.id, ]))

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
controls <-pData(gene_summaries)[pData(gene_summaries)$Dx == "CTL", "ShortName"]
sertad4.plot.data[sertad4.plot.data$SampleID %in% controls, "SampleType"] <- "Controls"
head(sertad4.plot.data)

# 3.8 Generate SERTAD4 plot
cols <- c(colors[1], "gray")
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/figures/boxplot_SERTAD4_gene_vs_transcripts_vs_exons_residuals.pdf", width = 15)
sp <- ggplot(sertad4.plot.data, aes(x = EnsemblID, y = Expression, fill = SampleType)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = cols) + 
  geom_hline(yintercept = median(sertad4.plot.data[sertad4.plot.data$EnsemblID == "SERTAD4", "Expression"]), linetype = "dashed") 
sp + facet_wrap(. ~ Type, scales = "free")
dev.off()


### Overlap DEA hits for p < 0.01
genes <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/results/1.1_diff_expr_genes_for_diagnosis.csv", select = c("p-value", "Gene Symbol"), data.table = F)
genes <- genes[genes$`p-value` < 0.01, "Gene Symbol"]
transcripts <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/results/1.3_diff_expr_transcripts_for_diagnosis_new_pvalue_and_fdr.csv", select = c("p-value", "Gene Symbol"), data.table = F)
transcripts <- transcripts[transcripts$`p-value` < 0.01, "Gene Symbol"]
exons <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/results/1.3_diff_expr_exons_for_diagnosis_new_pvalue_and_fdr.csv", select = c("p-value", "Gene Symbol"), data.table = F)
exons <- exons[exons$`p-value` < 0.01, "Gene Symbol"]

length(genes)
length(transcripts)
length(exons)
length(intersect(genes, transcripts))
length(intersect(genes, exons))
length(intersect(transcripts, exons))

library(VennDiagram)
draw.triple.venn(length(genes), length(transcripts), length(exons), 
                 length(intersect(genes, transcripts)), length(intersect(transcripts, exons)),
                 length(intersect(genes, exons)), length(intersect(exons, (intersect(genes, transcripts)))),
                 category = c("Gene", "Transcript", "Exon"), scaled = F)


