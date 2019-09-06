# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(tidyr)
  library(arules)
})

### 1. Load Background data ###
# We first need to get MAF values for each SNP using PLINK
# /storage/groups/ccm01/tools/PLINK/plink  --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --freq
# Now read in the MAF values for each SNP
mafs <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/plink.frq", header = TRUE)

# generate bins of MAF in 0.05 steps
boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
mafs$bin <- discretize(mafs$MAF, method = "fixed", breaks = boundaries, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# load background SNPs
snps <- read.table("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", header = F, stringsAsFactors = FALSE) # snp chr pos
snps <- snps[, c("V2", "V4")]
colnames(snps) <- c("SNP", "POS")

# merge maf with background SNPs
background <- merge(mafs, snps, by = "SNP")

# save into GRanges object
background.gr <- GRanges(seqnames = background$CHR,
                         ranges = IRanges(start = as.numeric(as.character(background$POS)), 
                                          end = as.numeric(as.character(background$POS))),
                         snp_id = background$SNP, 
                         bin = background$bin)
background.gr <- background.gr[!duplicated(background.gr)]  # 6,798,048 SNPs
rm(snps, background)

### 2. LOAD eQTL results ###
# cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

### 3. SNP-based anotation ###
### 3.1 Transform eQTL data into SNP-based GRanges objects ###
# Gene eQTLs
eqtls.gene <- merge(eqtls.gene, mafs, by = "SNP")
eqtls.gene.gr <- GRanges(seqnames = eqtls.gene$SNP_location, snp_id = eqtls.gene$SNP, gene_name = eqtls.gene$Gene_symbol, maf = eqtls.gene$MAF, bin = eqtls.gene$bin)
eqtls.gene.gr <- eqtls.gene.gr[!duplicated(eqtls.gene.gr)]  # 36,637 SNPs
# Gene cases eQTLs
eqtls.gene.cases <- merge(eqtls.gene.cases, mafs, by = "SNP")
eqtls.gene.cases.gr <- GRanges(seqnames = eqtls.gene.cases$SNP_location, snp_id = eqtls.gene.cases$SNP, gene_name = eqtls.gene.cases$Gene_symbol, maf = eqtls.gene.cases$MAF, bin = eqtls.gene.cases$bin)
eqtls.gene.cases.gr <- eqtls.gene.cases.gr[!duplicated(eqtls.gene.cases.gr)]  # 12,605 SNPs
# Transcript eQTLs
eqtls.transcript <- merge(eqtls.transcript, mafs, by = "SNP")
eqtls.transcript.gr <- GRanges(seqnames = eqtls.transcript$SNP_location, snp_id = eqtls.transcript$SNP, transcript_id = eqtls.transcript$Transcript, maf = eqtls.transcript$MAF, bin = eqtls.transcript$bin)
eqtls.transcript.gr <- eqtls.transcript.gr[!duplicated(eqtls.transcript.gr)]  # 43,039 SNPs
# Transcript cases eQTLs
eqtls.transcript.cases <- merge(eqtls.transcript.cases, mafs, by = "SNP")
eqtls.transcript.cases.gr <- GRanges(seqnames = eqtls.transcript.cases$SNP_location, snp_id = eqtls.transcript.cases$SNP, transcript_id = eqtls.transcript.cases$Transcript, maf = eqtls.transcript.cases$MAF, bin = eqtls.transcript.cases$bin)
eqtls.transcript.cases.gr <- eqtls.transcript.cases.gr[!duplicated(eqtls.transcript.cases.gr)]  # 16,237 SNPs
# Exon eQTLs
eqtls.exon <- merge(eqtls.exon, mafs, by = "SNP")
eqtls.exon.gr <- GRanges(seqnames = eqtls.exon$SNP_location, snp_id = eqtls.exon$SNP, exon_id = eqtls.exon$Exon, maf = eqtls.exon$MAF, bin = eqtls.exon$bin)
eqtls.exon.gr <- eqtls.exon.gr[!duplicated(eqtls.exon.gr)]  # 190,570 SNPs
# Exon cases eQTLs
eqtls.exon.cases <- merge(eqtls.exon.cases, mafs, by = "SNP")
eqtls.exon.cases.gr <- GRanges(seqnames = eqtls.exon.cases$SNP_location, snp_id = eqtls.exon.cases$SNP, exon_id = eqtls.exon.cases$Exon, maf = eqtls.exon.cases$MAF, bin = eqtls.exon.cases$bin)
eqtls.exon.cases.gr <- eqtls.exon.cases.gr[!duplicated(eqtls.exon.cases.gr)]  # 79,108 SNPs
# remove unnecessary data
rm(eqtls.gene, eqtls.gene.cases, eqtls.transcript, eqtls.transcript.cases, eqtls.exon, eqtls.exon.cases)

### 2.2 Public data ###
### 2.2.1 ChromHMM
hmm <- readRDS("/storage/groups/ccm01/datasets/2019_chromHMM/chromHMM_all_states.Rds")

# remove chromosome X, Y and M
seqlevels(hmm, pruning.mode="coarse") <- setdiff(seqlevels(hmm), c("chrX", "chrY", "chrM"))

# get ChromHMM brain data
hmm.brain <- hmm[hmm$code == "E073"]
seqlevelsStyle(hmm.brain) <- "NCBI"

# get ChromHMM blood data
hmm.blood <- hmm[hmm$code == "E062"]
seqlevelsStyle(hmm.blood) <- "NCBI"

# Test frequency of eQTLs on chromosomes
#eqtls.chr <- as.data.frame(table(seqnames(eqtls.transcript.cases.gr)))
#eqtls.hmm <- as.data.frame(table(seqnames(hmm.blood)))
#chromosomes <- merge(eqtls.chr, eqtls.hmm, by = "Var1")
#wilcox.test(chromosomes$Freq.x, chromosomes$Freq.y, alternative = "two.sided")

# 2.4 Annotation
#################
# function, to calculate the gene annotation to a public dataset, but with respect to the MAF bins
# own: own cis-eQTL data set as GRanges object
# public: a public cis-eQTL data set as GRanges object
# nperm: number of permutations (usually a 1000 times)
# background: all background genes coming from the gene expression (here: 14,048 unqiue genes)
annotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  # resample nperm times
  resampling <- lapply(1:3, function(x) sum(sapply(1:11, function(x) sum(overlapsAny(sample(background[background$bin == x], length(own[own$bin == x])), public)))))
  # IF 0 => echt kleiner als 1e^(-5/6)
  p.value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p.value), as.character(or))) # return the resampling p-value
}

results <- data.frame()
perm <- 10000

# ChromHMM states
states <- c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies")
#interesting.codes <- unique(hmm.blood$name)

# 2.4.1 annotation to ChromHMM Blood
#for (state in states) {
#  overlap <- paste0("ChromHMM_blood_", state, "_SNP_overlap")
#  pval <- paste0("ChromHMM_blood_", state, "_SNP_pVal")
#  or <- paste0("ChromHMM_blood_", state, "_SNP_or")
#  results[c(overlap, pval, or), "Gene"] <- annotatePublicData(eqtls.gene.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done gene: ", state))
#  results[c(overlap, pval, or), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done gene cases only: ", state))
#  results[c(overlap, pval, or), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done gtranscript: ", state))
#  results[c(overlap, pval, or), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done transcript cases only: ", state))
#  results[c(overlap, pval, or), "Exon"] <- annotatePublicData(eqtls.exon.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done exon: ", state))
#  results[c(overlap, pval, or), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, hmm.blood[hmm.blood$name == state], perm, background.gr)
#  print(paste0("Done exon cases only: ", state))
#}

# 2.4.6 write results into a file
#write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_ChromHMM_blood.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

results <- data.frame()
perm <- 10000

### 2.5 Cortex meQTL ###
# 2.5.1 Load meQTLs
cortex.meqtl <- fread("/storage/groups/ccm01/datasets/2019_cortex_meQTL/GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv")
cortex.meqtl <- cortex.meqtl[cortex.meqtl$p < 0.05, ]
cortex.meqtl.gr <- GRanges(seqnames = cortex.meqtl$snpChr, 
                       ranges = IRanges(start = as.numeric(as.character(cortex.meqtl$snpPos)), 
                                        end = as.numeric(as.character(cortex.meqtl$snpPos))),
                       snp_id = cortex.meqtl$snpRsNum)  # generate GRanges object only with gene ID, name and location
cortex.meqtl.gr <- cortex.meqtl.gr[!duplicated(cortex.meqtl.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
seqlevelsStyle(cortex.meqtl.gr) <- "NCBI"
rm(cortex.meqtl)

# 2.5.2 annotate meQTLs
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene"] <- annotatePublicData(eqtls.gene.gr, cortex.meqtl.gr, perm, background.gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Gene_cases"] <- annotatePublicData(eqtls.gene.cases.gr, cortex.meqtl.gr, perm, background.gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript"] <- annotatePublicData(eqtls.transcript.gr, cortex.meqtl.gr, perm, background.gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Transcript_cases"] <- annotatePublicData(eqtls.transcript.cases.gr, cortex.meqtl.gr, perm, background.gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon"] <- annotatePublicData(eqtls.exon.gr, cortex.meqtl.gr, perm, background.gr)
results[c("Cortex_meQTL_SNP_overlap", "Cortex_meQTL_SNP_pVal", "Cortex_meQTL_SNP_or"), "Exon_cases"] <- annotatePublicData(eqtls.exon.cases.gr, cortex.meqtl.gr, perm, background.gr)

# 2.5.3 write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/pmb_cis-eQTL_annotation_results_10000_permutations_meQTLs.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")




# load SCZ2 GWAS
gwas.scz <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/SCZ2/scz2.snp.results.txt")
gwas.scz <- gwas.scz[gwas.scz$p < 0.05, ]
gwas.scz.gr <- GRanges(seqnames = gwas.scz$hg19chrc, 
                       ranges = IRanges(start = as.numeric(as.character(gwas.scz$bp)), 
                                        end = as.numeric(as.character(gwas.scz$bp))),
                       snp_id = gwas.scz$snpid)  # generate GRanges object only with gene ID, name and location
gwas.scz.gr <- gwas.scz.gr[!duplicated(gwas.scz.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
seqlevelsStyle(gwas.scz.gr) <- "NCBI"
rm(gwas.scz)

sum(overlapsAny(eqtls.exon.cases.gr, gwas.scz.gr))
sum(overlapsAny(eqtls.exon.cases.gr, gwas.scz.gr))/(length(eqtls.exon.cases.gr) / 100)


# load BIP-2018 GWAS
gwas.bpd <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/BIP-2018/daner_PGC_BIP32b_mds7a_0416a")
gwas.bpd <- gwas.bpd[gwas.bpd$P < 0.05, ]
gwas.bpd.gr <- GRanges(seqnames = gwas.bpd$CHR, 
                       ranges = IRanges(start = as.numeric(as.character(gwas.bpd$BP)), 
                                        end = as.numeric(as.character(gwas.bpd$BP))),
                       snp_id = gwas.bpd$SNP)  # generate GRanges object only with gene ID, name and location
gwas.bpd.gr <- gwas.bpd.gr[!duplicated(gwas.bpd.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
rm(gwas.bpd)

sum(overlapsAny(eqtls.exon.cases.gr, gwas.bpd.gr))
sum(overlapsAny(eqtls.exon.cases.gr, gwas.bpd.gr))/(length(eqtls.exon.cases.gr) / 100)


# load MDD2 GWAS
gwas.mdd <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/MDD2/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe")
gwas.mdd <- gwas.mdd[gwas.mdd$P < 0.05, ]
gwas.mdd.gr <- GRanges(seqnames = gwas.mdd$CHR, 
                       ranges = IRanges(start = as.numeric(as.character(gwas.mdd$BP)), 
                                        end = as.numeric(as.character(gwas.mdd$BP))),
                       snp_id = gwas.mdd$SNP)  # generate GRanges object only with gene ID, name and location
gwas.mdd.gr <- gwas.mdd.gr[!duplicated(gwas.mdd.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
rm(gwas.mdd)

sum(overlapsAny(eqtls.exon.cases.gr, gwas.mdd.gr))
sum(overlapsAny(eqtls.exon.cases.gr, gwas.mdd.gr))/(length(eqtls.exon.cases.gr) / 100)


# load Cross-disorder GWAS
gwas.cd <- fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Cross_Disorder_Full/pgc.cross.full.2013-03.txt")
gwas.cd <- gwas.cd[gwas.cd$pval < 0.05, ]
gwas.cd.gr <- GRanges(seqnames = gwas.cd$hg18chr, 
                       ranges = IRanges(start = as.numeric(as.character(gwas.cd$bp)), 
                                        end = as.numeric(as.character(gwas.cd$bp))),
                       snp_id = gwas.cd$snpid)  # generate GRanges object only with gene ID, name and location
gwas.cd.gr <- gwas.cd.gr[!duplicated(gwas.cd.gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
rm(gwas.cd)

sum(overlapsAny(eqtls.exon.gr, gwas.cd.gr))
sum(overlapsAny(eqtls.exon.gr, gwas.cd.gr))/(length(eqtls.exon.gr) / 100)




