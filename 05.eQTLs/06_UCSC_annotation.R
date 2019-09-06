# Load locally installed packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  library(ChIPseeker)
  library(org.Hs.eg.db)
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

# remove unnecessary data
rm(eqtls.gene, eqtls.gene.cases, eqtls.transcript, eqtls.transcript.cases, eqtls.exon, eqtls.exon.cases)

# 3. LOAD public data
# BIOS whole blood gene eQTL SNPs
bios.gene <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/gene_level_eQTLs.txt", sep = "\t")
bios.gene[, length(bios.gene)] <- NULL # remove last column
bios.gene <- bios.gene[bios.gene$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.gene.gr <- GRanges(seqnames = bios.gene$SNPChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.gene$SNPChrPos)), 
                                         end = as.numeric(as.character(bios.gene$SNPChrPos))),
                        snp_id = bios.gene$SNPName)  # generate GRanges object
bios.gene.gr <- bios.gene.gr[!duplicated(bios.gene.gr)]  # remove duplicated rows: 5,526,636 => 1,957,252 SNPs
seqlevelsStyle(bios.gene.gr) <- "UCSC"
rm(bios.gene)

# BIOS whole blood exon eQTL SNPs
bios.exon <- fread("/storage/groups/ccm01/datasets/2019_BIOS_QTL_browser/eQTL_results/exon_level_eQTLs_independent_effects.txt", sep = "\t")
bios.exon[, length(bios.exon)] <- NULL # remove last column
bios.exon <- bios.exon[bios.exon$FDR < 0.05, ]  # filter for eQTLs with FDR < 0.05
bios.exon.gr <- GRanges(seqnames = bios.exon$SNPChr,
                        ranges = IRanges(start = as.numeric(as.character(bios.exon$SNPChrPos)), 
                                         end = as.numeric(as.character(bios.exon$SNPChrPos))),
                        snp_id = bios.exon$SNPName)  # generate GRanges object
bios.exon.gr <- bios.exon.gr[!duplicated(bios.exon.gr)]  # remove duplicated rows: 5,526,636 => 1,957,252 SNPs
seqlevelsStyle(bios.exon.gr) <- "UCSC"
rm(bios.exon)

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


# 4. Annotate data
as.gene <- annotatePeak(unique(eqtls.gene.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.gene.cases <- annotatePeak(unique(eqtls.gene.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.transcript <- annotatePeak(unique(eqtls.transcript.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.transcript.cases <- annotatePeak(unique(eqtls.transcript.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.exon <- annotatePeak(unique(eqtls.exon.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.exon.cases <- annotatePeak(unique(eqtls.exon.cases.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.bios <- annotatePeak(unique(bios.gene.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.bios.exon <- annotatePeak(unique(bios.exon.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.cmc <- annotatePeak(unique(cmc.gene.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)
as.cmc.iso <- annotatePeak(unique(cmc.iso.gr), TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene, annoDb = "org.Hs.eg.db", verbose = F)

# 1.4 Plot the annotations
# 1.4.1 Combine all annotations
genomic.annotation <- bind_rows(Gene = as.gene@annoStat,
                                Gene_cases = as.gene.cases@annoStat,
                                Transcript = as.transcript@annoStat,
                                Transcript_cases = as.transcript.cases@annoStat,
                                Exon = as.exon@annoStat,
                                Exon_cases = as.exon.cases@annoStat,
                                Bios = as.bios@annoStat,
                                Bios_exon = as.bios.exon@annoStat,
                                CMC = as.cmc@annoStat,
                                CMC_iso = as.cmc.iso@annoStat,
                                .id = 'Type')
genomic.annotation$Type <- as.factor(genomic.annotation$Type)
# relevel to gene, transcript, exon order
genomic.annotation$Type <- factor(genomic.annotation$Type, levels = c("CMC_iso", "CMC", "Bios_exon", "Bios", "Exon_cases", "Exon", "Transcript_cases", "Transcript", "Gene_cases", "Gene"))

# 1.4.2 Plot Annotations as stacked barplot
ggplot(genomic.annotation, aes(Type, Frequency, fill = Feature)) + 
  geom_bar(stat = "identity", position = 'fill') + 
  scale_fill_brewer(palette = "Paired", direction = -1) +
  coord_flip() + 
  labs(title = "UCSC Genomic Features", y = "Frequency", x = "") +
  theme_minimal() +
  theme(text = element_text(size = 22))
ggsave("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/figures/02_Barplot_stacked_eqtl_annotations.pdf", height = 8, width = 11)

