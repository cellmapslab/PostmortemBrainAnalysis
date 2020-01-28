# Load locally installed packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
library(oligo)

# 1. After rma (median polish) to coreSummaries and probesetSummaries
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/2_affymetrix_exon_microarray_data_after_rma_core_summaries.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/2_affymetrix_exon_microarray_data_after_rma_probeset_summaries.RData")

fkbp5.probeset.ids <- c("2951573", "2951574", "2951575", "2951576", "2951577", "2951578", "2951579", "2951581", "2951582", "2951583", "2951587", "2951589", "2951590", "2951593", "2951595", "2951596", "2951601", "2951603", "2951630", "2951635", "2951640", "2951656", "2951657")
fkbp5.transcript.ids <- c("2951567", "2951655", "3401099", "3933039")
fkbp5.gene.id <- "ENSG00000096060.10"

fkbp5.probesets <- exprs(probesetSummaries)[rownames(exprs(probesetSummaries)) %in% fkbp5.probeset.ids, ]
fkbp5.probesets.median <- t(as.data.frame(colMedians(fkbp5.probesets)))
colnames(fkbp5.probesets.median) <- colnames(fkbp5.probesets)
rownames(fkbp5.probesets.median) <- "fkbp5.probesetSummaries"

fkbp5.core <- exprs(coreSummaries)[rownames(exprs(coreSummaries)) %in% fkbp5.transcript.ids, ]
fkbp5.core.median <- t(as.data.frame(colMedians(fkbp5.core)))
colnames(fkbp5.core.median) <- colnames(fkbp5.core)
rownames(fkbp5.core.median) <- "fkbp5.coreSummaries"

identical(colnames(fkbp5.core.median), colnames(fkbp5.probesets.median))
fkbp5.core.probeset <- rbind(fkbp5.core.median, fkbp5.probesets.median)

plot(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# 2. After rma and control/cross-hybridization/chrX/chrY/chrM/multi-mapped IDs filtering
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/3_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/3_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered.RData")

fkbp5.probeset.ids <- c("2951573", "2951574", "2951575", "2951576", "2951577", "2951578", "2951579", "2951581", "2951582", "2951583", "2951587", "2951589", "2951590", "2951593", "2951595", "2951596", "2951601", "2951603", "2951630", "2951635", "2951640", "2951656", "2951657")
fkbp5.transcript.ids <- c("2951567", "2951655", "3401099", "3933039")

fkbp5.probesets <- exprs(probesetSummaries)[rownames(exprs(probesetSummaries)) %in% fkbp5.probeset.ids, ]
fkbp5.probesets.median <- t(as.data.frame(colMedians(fkbp5.probesets)))
colnames(fkbp5.probesets.median) <- colnames(fkbp5.probesets)
rownames(fkbp5.probesets.median) <- "fkbp5.probesetSummaries"

fkbp5.core <- exprs(coreSummaries)[rownames(exprs(coreSummaries)) %in% fkbp5.transcript.ids, ]
fkbp5.core.median <- t(as.data.frame(colMedians(fkbp5.core)))
colnames(fkbp5.core.median) <- colnames(fkbp5.core)
rownames(fkbp5.core.median) <- "fkbp5.coreSummaries"

identical(colnames(fkbp5.core.median), colnames(fkbp5.probesets.median))
fkbp5.core.probeset <- rbind(fkbp5.core.median, fkbp5.probesets.median)

plot(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# 3. After batch correction
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_batch_corrected.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/4_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected.RData")

fkbp5.probeset.ids <- c("2951573", "2951574", "2951575", "2951576", "2951577", "2951578", "2951579", "2951581", "2951582", "2951583", "2951587", "2951589", "2951590", "2951593", "2951595", "2951596", "2951601", "2951603", "2951630", "2951635", "2951640", "2951656", "2951657")
fkbp5.transcript.ids <- c("2951567", "2951655", "3401099", "3933039")

fkbp5.probesets <- exprs(probesetSummaries)[rownames(exprs(probesetSummaries)) %in% fkbp5.probeset.ids, ]
fkbp5.probesets.median <- t(as.data.frame(colMedians(fkbp5.probesets)))
colnames(fkbp5.probesets.median) <- colnames(fkbp5.probesets)
rownames(fkbp5.probesets.median) <- "fkbp5.probesetSummaries"

fkbp5.core <- exprs(coreSummaries)[rownames(exprs(coreSummaries)) %in% fkbp5.transcript.ids, ]
fkbp5.core.median <- t(as.data.frame(colMedians(fkbp5.core)))
colnames(fkbp5.core.median) <- colnames(fkbp5.core)
rownames(fkbp5.core.median) <- "fkbp5.coreSummaries"

identical(colnames(fkbp5.core.median), colnames(fkbp5.probesets.median))
fkbp5.core.probeset <- rbind(fkbp5.core.median, fkbp5.probesets.median)

plot(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# 4. After SVA
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/5_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_batch_corrected_SVA.RData")
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/5_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA.RData")

fkbp5.probeset.ids <- c("2951573", "2951574", "2951575", "2951576", "2951577", "2951578", "2951579", "2951581", "2951582", "2951583", "2951587", "2951589", "2951590", "2951593", "2951595", "2951596", "2951601", "2951603", "2951630", "2951635", "2951640", "2951656", "2951657")
fkbp5.transcript.ids <- c("2951567", "2951655", "3401099", "3933039")

fkbp5.probesets <- exprs(probesetSummaries)[rownames(exprs(probesetSummaries)) %in% fkbp5.probeset.ids, ]
fkbp5.probesets.median <- t(as.data.frame(colMedians(fkbp5.probesets)))
colnames(fkbp5.probesets.median) <- colnames(fkbp5.probesets)
rownames(fkbp5.probesets.median) <- "fkbp5.probesetSummaries"

fkbp5.core <- exprs(coreSummaries)[rownames(exprs(coreSummaries)) %in% fkbp5.transcript.ids, ]
fkbp5.core.median <- t(as.data.frame(colMedians(fkbp5.core)))
colnames(fkbp5.core.median) <- colnames(fkbp5.core)
rownames(fkbp5.core.median) <- "fkbp5.coreSummaries"

identical(colnames(fkbp5.core.median), colnames(fkbp5.probesets.median))
fkbp5.core.probeset <- rbind(fkbp5.core.median, fkbp5.probesets.median)

plot(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# 5. After median to genes
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_batch_corrected_SVA_gene_summaries.RData")
coreSummaries <- exprs(geneSummaries)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/6_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries.RData")
probesetSummaries <- exprs(geneSummaries)

fkbp5.probesets <- t(as.data.frame(probesetSummaries[rownames(probesetSummaries) == "ENSG00000096060.10", ]))
fkbp5.core <- t(as.data.frame(coreSummaries[rownames(coreSummaries) == "ENSG00000096060.10", ]))

identical(colnames(fkbp5.core), colnames(fkbp5.probesets))
fkbp5.core.probeset <- rbind(fkbp5.core, fkbp5.probesets)

plot(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])


# 6. After annotation
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_core_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
coreSummaries <- exprs(geneSummaries)
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEX/RData_Files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_gene_summaries_annotated.RData")
probeSummaries <- exprs(geneSummaries)
rm(geneSummaries)

# FKBP5
fkbp5.probeset <- as.data.frame(probeSummaries[rownames(probeSummaries) == "ENSG00000096060.10", ])
fkbp5.core <- as.data.frame(coreSummaries[rownames(coreSummaries) == "ENSG00000096060.10", ])
rownames(fkbp5.core) <- "FKBP5.core"

identical(colnames(fkbp5.core), colnames(fkbp5.probeset))
fkbp5.core.probeset <- rbind(fkbp5.core, fkbp5.probeset)
fkbp5.core.probeset <- t(t(fkbp5.core.probeset))

plot(as.vector(fkbp5.core.probeset[2,]), as.vector(fkbp5.core.probeset[1,]))
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# SERTAD4
fkbp5.probeset <- as.data.frame(probeSummaries[rownames(probeSummaries) == "ENSG00000082497.7", ])
fkbp5.core <- as.data.frame(coreSummaries[rownames(coreSummaries) == "ENSG00000082497.7", ])
rownames(fkbp5.core) <- "FKBP5.core"

identical(colnames(fkbp5.core), colnames(fkbp5.probeset))
fkbp5.core.probeset <- rbind(fkbp5.core, fkbp5.probeset)
fkbp5.core.probeset <- t(t(fkbp5.core.probeset))

plot(as.vector(fkbp5.core.probeset[2,]), as.vector(fkbp5.core.probeset[1,]))
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])

# EPM2A
fkbp5.probeset <- as.data.frame(probeSummaries[rownames(probeSummaries) == "ENSG00000112425.9", ])
fkbp5.core <- as.data.frame(coreSummaries[rownames(coreSummaries) == "ENSG00000112425.9", ])
rownames(fkbp5.core) <- "FKBP5.core"

identical(colnames(fkbp5.core), colnames(fkbp5.probeset))
fkbp5.core.probeset <- rbind(fkbp5.core, fkbp5.probeset)
fkbp5.core.probeset <- t(t(fkbp5.core.probeset))

plot(as.vector(fkbp5.core.probeset[2,]), as.vector(fkbp5.core.probeset[1,]))
cor(fkbp5.core.probeset[1,], fkbp5.core.probeset[2,])
