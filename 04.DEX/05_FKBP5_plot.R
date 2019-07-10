library(oligo)
library(rtracklayer)
library(ggplot2)
library(viridis)

load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_exon_summaries_annotated_combat_batch_corrected.RData") # load exons (after batch correction)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/2_affymetrix_exon_microarray_data_after_rma_core_gene_summaries.RData") # load genes before annotation and remove batch afterwards
exprs_rma <- exprs(geneSummaries)
batch <- pData(geneSummaries)$Batch
modcombat <- model.matrix(~1, data = pData(geneSummaries))
combat_edata <- ComBat(dat = exprs_rma, batch = batch, mod = modcombat, par.prior = TRUE, prior.plots = FALSE)
identical(rownames(exprs(geneSummaries)), rownames(combat_edata)) # check if the two matrices are in the same order
exprs(geneSummaries) <- combat_edata
dim(geneSummaries)

# Load Affymetrix Huex1.0stv2 probeset annotation file (downloaded from the Affymetrix website)
huex <- read.csv("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/HuEx-1_0-st-v2/HuEx-1_0-st-v2.na36.hg19.probeset.csv", comment.char = "#")
huex <- huex[huex$start != "---", ]  # remove ids with no start and stop information

### Generate Boxplot for FKBP5 Gene ###
transcript.ids <- huex[grep("FKBP5", huex$gene_assignment), "transcript_cluster_id"]
#exon.ids <- fData(exonSummaries)[fData(exonSummaries)$gene_name == "FKBP5", "exon_id"]
#exon.ids <- c("ENSE00002227265.1", "ENSE00001424352.1", "ENSE00001370159.1", "ENSE00000747346.1", "ENSE00001391213.1", "ENSE00000747342.1", "ENSE00003659116.1", "ENSE00000747304.1", "ENSE00000747303.1", "ENSE00000747302.1", "ENSE00000747301.1", "ENSE00001193215.7") # FKBP5-201
exon.ids <- c("ENSE00002312621.1", "ENSE00001370159.1", "ENSE00000747346.1", "ENSE00001391213.1", "ENSE00000747342.1", "ENSE00003659116.1", "ENSE00000747304.1", "ENSE00000747303.1", "ENSE00000747302.1", "ENSE00000747301.1", "ENSE00001193215.7" )
fkbp.transcript.exprs <- exprs(geneSummaries)[rownames(exprs(geneSummaries)) %in% transcript.ids, ]
fkbp.exon.exprs <- exprs(exonSummaries)[rownames(exprs(exonSummaries)) %in% exon.ids, ]
fkbp <- rbind(fkbp.transcript.exprs, fkbp.exon.exprs)
rownames(fkbp) <- c("FKBP5", rownames(fkbp)[-1])
#par(mar = c(10.1, 4.3, 2.1, 2.1))
#boxplot(t(fkbp), ylab = "Gene expression", cex.lab = 1.5, las =  2)

# load human gencode file
genes <- import.gff("/storage/groups/ccm01/datasets/2018_MPIP-Dean_postmortem_brain/annotation_files/gencode.v19.annotation.gtf", format = "gtf")
genes <- genes[genes$gene_type == "protein_coding"]  # get only protein-coding genes
genes <- genes[genes$type == "exon"]  # get only exons from the genes
genes <- as.data.frame(genes)
exon.mapping <- genes[genes$exon_id %in% exon.ids, c("exon_id", "exon_number", "transcript_name")]
exon.mapping <- exon.mapping[exon.mapping$transcript_name == "FKBP5-202",]
exon.mapping <- exon.mapping[!duplicated(exon.mapping),]
exon.mapping <- rbind(exon.mapping, c("FKBP5", 0, "FKBP5-202"))
exon.mapping$exon_number <- as.numeric(exon.mapping$exon_number)
exon.mapping <- exon.mapping[order(exon.mapping$exon_number),]
fkbp5 <- merge(exon.mapping, fkbp, by.x = "exon_id", by.y = "row.names", all.x = T)
fkbp5.plot.data <- fkbp5[rep(seq_len(nrow(fkbp5)), 204), c("exon_id", "exon_number")]
fkbp5.plot.data$expression <- unlist(fkbp5[4:207], use.names = FALSE)
fkbp5.plot.data$exon_number <- as.factor(fkbp5.plot.data$exon_number)
fkbp5.plot.data$exon_id <- as.factor(fkbp5.plot.data$exon_id)
#fkbp5.plot.data$new_id <- paste(fkbp5.plot.data$exon_id, "_", fkbp5.plot.data$exon_number, sep = "")
#fkbp5.plot.data$new_id <- as.factor(fkbp5.plot.data$new_id)
fkbp5.plot.data$exon_id <- factor(fkbp5.plot.data$exon_id, levels = c("FKBP5", "ENSE00002312621.1", "ENSE00001370159.1", "ENSE00000747346.1", "ENSE00001391213.1", "ENSE00000747342.1", "ENSE00003659116.1", "ENSE00000747304.1", "ENSE00000747303.1", "ENSE00000747302.1", "ENSE00000747301.1", "ENSE00001193215.7"))

pdf("/storage/scratch/users/karolina.worf/MPIP/04.DEX/figures/boxplot_FKBP5_gene_vs_exons.pdf", width = 15)
ggplot(fkbp5.plot.data, aes(x = exon_id, y = expression, fill = exon_number)) +
  theme_classic() +
  theme(text = element_text(size = 22), axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_boxplot(alpha = 0.7) +
  xlab("") +
  ylab("Gene expression") +
  scale_fill_manual(values = c("darkblue", viridis(11)), labels = c("gene", "exon 1", "exon 2", "exon 3", "exon 4", "exon 5", "exon 6", "exon 7", "exon 8", "exon 9", "exon 10", "exon 11", "exon 12"), name = "") + 
  scale_x_discrete(labels = c("FKBP5", "ENSE00002312621.1", "ENSE00001370159.1", "ENSE00000747346.1", "ENSE00001391213.1", "ENSE00000747342.1", "ENSE00003659116.1", "ENSE00000747304.1", "ENSE00000747303.1", "ENSE00000747302.1", "ENSE00000747301.1", "ENSE00001193215.7")) + 
  geom_hline(yintercept = median(fkbp5.plot.data[fkbp5.plot.data$exon_id == "FKBP5", "expression"]), linetype = "dashed") 
dev.off()

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
fkbp.cases <- as.data.frame(fkbp.transcript.exprs[which(names(fkbp.transcript.exprs) %in% cases)])
fkbp.ctrls <- as.data.frame(fkbp.transcript.exprs[which(names(fkbp.transcript.exprs) %in% ctrl)])
fkbp.cases <- cbind(fkbp.cases, rep("Cases", 107))
colnames(fkbp.cases) <- c("Expression", "Dx")
fkbp.ctrls <- cbind(fkbp.ctrls, rep("Controls", 62))
colnames(fkbp.ctrls) <- c("Expression", "Dx")
fkbp <- rbind(fkbp.cases, fkbp.ctrls)
ggplot(fkbp, aes(x = Dx, y = Expression, fill = Dx)) + 
  geom_boxplot() +
  theme_grey(base_size = 18)
