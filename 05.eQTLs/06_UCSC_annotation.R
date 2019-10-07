# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(ChIPseeker)
  library(org.Hs.eg.db)
  library(dplyr)
  library(ggplot2)
})

### 1. LOAD eQTL results ###
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
length(eqtls.gene$SNP)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)

# cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

# get cases only unique eQTLs
unique.eqtl.ids <- setdiff(eqtls.exon.cases[, c("SNP", "Exon")], eqtls.exon[, c("SNP", "Exon")])
eqtls.exon.cases <- as.data.frame(eqtls.exon.cases)

eqtls.exon.cases.unique <- data.frame()
for (i in 1:length(unique.eqtl.ids$SNP)) {
  eqtls.exon.cases.unique <- rbind(eqtls.exon.cases.unique, filter(eqtls.exon.cases, SNP %in% unique.eqtl.ids[i, "SNP"] & Exon %in% unique.eqtl.ids[i, "Exon"]))
}
head(eqtls.exon.cases.unique)

### 2. Transform eQTL data into SNP-based GRanges objects ###
# Gene eQTLs
eqtls.gene.gr <- GRanges(seqnames = eqtls.gene$SNP_location, snp_id = eqtls.gene$SNP)
eqtls.gene.gr <- eqtls.gene.gr[!duplicated(eqtls.gene.gr)]  # 36,637 SNPs
seqlevelsStyle(eqtls.gene.gr) <- "UCSC"

# Gene cases eQTLs
eqtls.gene.cases.gr <- GRanges(seqnames = eqtls.gene.cases$SNP_location, snp_id = eqtls.gene.cases$SNP)
eqtls.gene.cases.gr <- eqtls.gene.cases.gr[!duplicated(eqtls.gene.cases.gr)]  # 12,605 SNPs
seqlevelsStyle(eqtls.gene.cases.gr) <- "UCSC"

# Transcript eQTLs
eqtls.transcript.gr <- GRanges(seqnames = eqtls.transcript$SNP_location, snp_id = eqtls.transcript$SNP)
eqtls.transcript.gr <- eqtls.transcript.gr[!duplicated(eqtls.transcript.gr)]  # 43,039 SNPs
seqlevelsStyle(eqtls.transcript.gr) <- "UCSC"

# Transcript cases eQTLs
eqtls.transcript.cases.gr <- GRanges(seqnames = eqtls.transcript.cases$SNP_location, snp_id = eqtls.transcript.cases$SNP)
eqtls.transcript.cases.gr <- eqtls.transcript.cases.gr[!duplicated(eqtls.transcript.cases.gr)]  # 16,237 SNPs
seqlevelsStyle(eqtls.transcript.cases.gr) <- "UCSC"

# Exon eQTLs
eqtls.exon.gr <- GRanges(seqnames = eqtls.exon$SNP_location, snp_id = eqtls.exon$SNP)
eqtls.exon.gr <- eqtls.exon.gr[!duplicated(eqtls.exon.gr)]  # 190,570 SNPs
seqlevelsStyle(eqtls.exon.gr) <- "UCSC"

# Exon cases eQTLs
eqtls.exon.cases.gr <- GRanges(seqnames = eqtls.exon.cases$SNP_location, snp_id = eqtls.exon.cases$SNP)
eqtls.exon.cases.gr <- eqtls.exon.cases.gr[!duplicated(eqtls.exon.cases.gr)]  # 79,108 SNPs
seqlevelsStyle(eqtls.exon.cases.gr) <- "UCSC"

# Exon cases unique eQTLs
eqtls.exon.cases.unique.gr <- GRanges(seqnames = eqtls.exon.cases.unique$SNP_location, snp_id = eqtls.exon.cases.unique$SNP)
eqtls.exon.cases.unique.gr <- eqtls.exon.cases.gr[!duplicated(eqtls.exon.cases.unique.gr)]  # 79,108 SNPs
seqlevelsStyle(eqtls.exon.cases.unique.gr) <- "UCSC"

# remove unnecessary data
rm(eqtls.gene, eqtls.gene.cases, eqtls.transcript, eqtls.transcript.cases, eqtls.exon, eqtls.exon.cases)

# 3. LOAD public data
# CMC (Fromer et al.) gene eQTLs
cmc.gene <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc.gene <- cmc.gene[cmc.gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.gene.gr <- GRanges(seqnames = cmc.gene$SNP_location,
                       snp_id = cmc.gene$SNP)  # generate GRanges object 
cmc.gene.gr <- cmc.gene.gr[!duplicated(cmc.gene.gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
rm(cmc.gene)

# CMC (Fromer et al.) isoforms eQTLs
cmc.iso <- fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt")
cmc.iso <- cmc.iso[cmc.iso$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc.iso.gr <- GRanges(seqnames = cmc.iso$SNP_location,
                      snp_id = cmc.iso$SNP)  # generate GRanges object 
cmc.iso.gr <- cmc.iso.gr[!duplicated(cmc.iso.gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
rm(cmc.iso)

### GTEx Whole Blood
gtex.wb <- fread("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/Whole_Blood.v7.egenes.txt")
gtex.wb <- gtex.wb[gtex.wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
gtex.wb.gr <- GRanges(seqnames = gtex.wb$chr, 
                      ranges = IRanges(start = as.numeric(as.character(gtex.wb$pos)), 
                                       end = as.numeric(as.character(gtex.wb$pos))),
                      snp_id = gtex.wb$rs_id_dbSNP147_GRCh37p13)  # generate GRanges object only with gene ID, name and location
gtex.wb.gr <- gtex.wb.gr[!duplicated(gtex.wb.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
seqlevelsStyle(gtex.wb.gr) <- "UCSC"
rm(gtex.wb)

### GTEx Brain (all brain tissues)
brains <- c("Brain_Amygdala.v7.egenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v7.egenes.txt", "Brain_Caudate_basal_ganglia.v7.egenes.txt", "Brain_Cerebellar_Hemisphere.v7.egenes.txt", "Brain_Cerebellum.v7.egenes.txt", "Brain_Cortex.v7.egenes.txt", "Brain_Frontal_Cortex_BA9.v7.egenes.txt", "Brain_Hippocampus.v7.egenes.txt", "Brain_Hypothalamus.v7.egenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v7.egenes.txt", "Brain_Putamen_basal_ganglia.v7.egenes.txt", "Brain_Spinal_cord_cervical_c-1.v7.egenes.txt", "Brain_Substantia_nigra.v7.egenes.txt")
gtex.brains <- data.frame()
for (i in brains) {
  filename <- paste0("/storage/groups/ccm01/datasets/2018_GTExv7_cis_eQTL/GTEx_Analysis_v7_eQTL/", i)
  data <- fread(filename)
  data <- data[data$qval < 0.05, ]  # filter for significant hits
  gtex.brains <- rbind(gtex.brains, data) # merge different brain files together
}
gtex.brains.gr <- GRanges(seqnames = gtex.brains$chr, 
                          ranges = IRanges(start = as.numeric(as.character(gtex.brains$pos)), 
                                           end = as.numeric(as.character(gtex.brains$pos))),
                          snp_id = gtex.brains$rs_id_dbSNP147_GRCh37p13)  # generate GRanges object only with gene ID, name and location
gtex.brains.gr <- gtex.brains.gr[!duplicated(gtex.brains.gr)]  # remove duplicated rows: 58,589 => 17,373 genes
seqlevelsStyle(gtex.brains.gr) <- "UCSC"
rm(brains, filename, data, i, gtex.brains)


# 4. Annotate data
as.gene <- annotatePeak(unique(eqtls.gene.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.gene.cases <- annotatePeak(unique(eqtls.gene.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.transcript <- annotatePeak(unique(eqtls.transcript.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.transcript.cases <- annotatePeak(unique(eqtls.transcript.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.exon <- annotatePeak(unique(eqtls.exon.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.exon.cases <- annotatePeak(unique(eqtls.exon.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.exon.cases.unique <- annotatePeak(unique(eqtls.exon.cases.unique.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.cmc <- annotatePeak(unique(cmc.gene.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.cmc.iso <- annotatePeak(unique(cmc.iso.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.gtex.wb <- annotatePeak(unique(gtex.wb.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.gtex.brain <- annotatePeak(unique(gtex.brains.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)

# 1.4 Plot the annotations
# 1.4.1 Combine all annotations
genomic.annotation <- bind_rows(Gene = as.gene@annoStat,
                                Gene_cases = as.gene.cases@annoStat,
                                Transcript = as.transcript@annoStat,
                                Transcript_cases = as.transcript.cases@annoStat,
                                Exon = as.exon@annoStat,
                                Exon_cases = as.exon.cases@annoStat,
                                Exon_cases_unique = as.exon.cases.unique@annoStat,
                                CMC = as.cmc@annoStat,
                                CMC_iso = as.cmc.iso@annoStat,
                                GTEx_blood = as.gtex.wb@annoStat,
                                GTEx_brain = as.gtex.brain@annoStat,
                                .id = 'Type')
genomic.annotation$Type <- as.factor(genomic.annotation$Type)
# relevel to gene, transcript, exon order
genomic.annotation$Type <- factor(genomic.annotation$Type, levels = c("GTEx_brain", "GTEx_blood", "CMC_iso", "CMC", "Exon_cases_unique", "Exon_cases", "Exon", "Transcript_cases", "Transcript", "Gene_cases", "Gene"))

# 1.4.2 Plot Annotations as stacked barplot
ggplot(genomic.annotation, aes(Type, Frequency, fill = Feature)) + 
  geom_bar(stat = "identity", position = 'fill') + 
  scale_fill_brewer(palette = "Paired", direction = -1) +
  coord_flip() + 
  labs(title = "UCSC Genomic Features", y = "Frequency", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 22))
ggsave("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/02_Barplot_stacked_eqtl_annotations.pdf", height = 8, width = 11)

genomic.annotation

