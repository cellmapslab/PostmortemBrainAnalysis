### R-Script to generate GRanges object for our eQTL data as well as several public datasets (CMC, GTEx, VEP) ###
# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)
library(GenomicRanges)
library(GenomeInfoDb)
library(rtracklayer)
library(IRanges)
library(arules)

# function to generate an eQTL GRanges object (with MAF bins)
# @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
# @param mdata Merging data (e.g. Background SNPs)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObject <- function(input, mdata, ofile){
    # merge with MAF info
    input <- merge(input, mdata, by = "SNP")
    # generate GRanges object
    input_gr <- GenomicRanges::GRanges(seqnames = input$CHR, 
                                       ranges = IRanges::IRanges(start = as.numeric(as.character(input$POS)), 
                                                                 end = as.numeric(as.character(input$POS))),
                                       snp_id = input$SNP, 
                                       bin = input$bin)
    # remove duplicates
    input_gr <- input_gr[!duplicated(input_gr)] 
    # save GRanges object to an rds file
    saveRDS(input_gr, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/", ofile))
    # return finalized GRanges object
    return(input_gr)
}

### Background SNPs
# We first need to get MAF values for each SNP using PLINK
# /storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged --freq --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/annotation/pmb_info06_all_chr_maf0.05_merged_maf_information
# Now read in the MAF values for each SNP
mafs <- data.table::fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/annotation/pmb_info06_all_chr_maf0.05_merged_maf_information.frq", data.table = F)

# generate bins of MAF in 0.05 steps
boundaries <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 1.0)
mafs$bin <- arules::discretize(mafs$MAF, method = "fixed", breaks = boundaries, labels=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11))

# load background SNPs
snps <- data.table::fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged.bim", data.table = F, select = c(2, 4), col.names = c("SNP", "POS")) # snp chr pos

# Generate background SNP data as GRanges object
background_gr <- GenerateGrangesObject(input = snps, 
                                       mdata = mafs, 
                                       ofile = "background_SNPs.rds")
length(background_gr$snp_id) # 6,798,048
rm(snps)

# Generate merging file for eQTL data from background SNPs 
background <- as.data.frame(background_gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")
head(background)

### Gene eQTLs
eqtls_gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls_gene_gr <- GenerateGrangesObject(input = eqtls_gene,
                                       mdata = background,
                                       ofile = "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(eqtls_gene_gr$snp_id) # 52,496 SNPs
rm(eqtls_gene, eqtls_gene_gr)

### Gene cases only eQTLs
eqtls_gene_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls_gene_cases_gr <- GenerateGrangesObject(input = eqtls_gene_cases,
                                             mdata = background,
                                             ofile = "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds")
length(eqtls_gene_cases_gr$snp_id) # 18,944 SNPs
rm(eqtls_gene_cases, eqtls_gene_cases_gr)

### Transcript eQTLs
eqtls_transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls_transcript_gr <- GenerateGrangesObject(input = eqtls_transcript,
                                             mdata = background,
                                             ofile = "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(eqtls_transcript_gr$snp_id) # 112,861 SNPs
rm(eqtls_transcript, eqtls_transcript_gr)

### Transcript cases only eQTLs
eqtls_transcript_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls_transcript_cases_gr <- GenerateGrangesObject(input = eqtls_transcript_cases,
                                                   mdata = background,
                                                   ofile = "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds")
length(eqtls_transcript_cases_gr$snp_id) # 40,791 SNPs
rm(eqtls_transcript_cases, eqtls_transcript_cases_gr)

### Exon eQTLs
eqtls_exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls_exon_gr <- GenerateGrangesObject(input = eqtls_exon,
                                       mdata = background,
                                       ofile = "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(eqtls_exon_gr$snp_id) # 190,313 SNPs
rm(eqtls_exon, eqtls_exon_gr)

### Exon cases only eQTLs
eqtls_exon_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls_exon_cases_gr <- GenerateGrangesObject(input = eqtls_exon_cases,
                                             mdata = background,
                                             ofile = "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds")
length(eqtls_exon_cases_gr$snp_id) # 79,524 SNPs
rm(eqtls_exon_cases, eqtls_exon_cases_gr)

### CMC (Fromer et al.) gene eQTLs
cmc_gene <- data.table::fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc_gene <- cmc_gene[cmc_gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc_gene_gr <- GenomicRanges::GRanges(seqnames = cmc_gene$SNP_location, 
                                      snp_id = cmc_gene$SNP)  # generate GRanges object 
cmc_gene_gr <- cmc_gene_gr[!duplicated(cmc_gene_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
GenomeInfoDb::seqlevelsStyle(cmc_gene_gr) <- "NCBI"
saveRDS(cmc_gene_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/CMC_adjustedSVA_eQTL_SNPs.rds")
rm(cmc_gene, cmc_gene_gr)

### CMC (Fromer et al.) isoform eQTLs
cmc_iso <- data.table::fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt")
cmc_iso <- cmc_iso[cmc_iso$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc_iso_gr <- GenomicRanges::GRanges(seqnames = cmc_iso$SNP_location,
                                     snp_id = cmc_iso$SNP)  # generate GRanges object 
cmc_iso_gr <- cmc_iso_gr[!duplicated(cmc_iso_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
GenomeInfoDb::seqlevelsStyle(cmc_iso_gr) <- "NCBI"
saveRDS(cmc_iso_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/CMC_adjustedSVA_isoQTL_SNPs.rds")
rm(cmc_iso, cmc_iso_gr)

### GTEx v8 is based on the new human genome GRCh38, as we are working on GRCH37/hg19 we need to downlift to this genome
# 3. Load hg38 to hg19 chain (chain downloaded from UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)
hg38_to_hg19 <- rtracklayer::import.chain("/storage/groups/ccm01/datasets/2019_GTEx_v8/hg38ToHg19.over.chain")

### GTEx Whole Blood eQTL
gtex_wb <- data.table::fread("/storage/groups/ccm01/datasets/2019_GTEx_v8/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.egenes.txt")
gtex_wb <- gtex_wb[gtex_wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
gtex_wb_gr <- GenomicRanges::GRanges(seqnames = gtex_wb$chr, 
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_wb$variant_pos)), 
                                                               end = as.numeric(as.character(gtex_wb$variant_pos))),
                                     snp_id = gtex_wb$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(gtex_wb_gr) <- "UCSC"  # UCSC needed for the lift over
gtex_wb_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_wb_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 12,360 => 12,292
gtex_wb_gr_hg19 <- gtex_wb_gr_hg19[!duplicated(gtex_wb_gr_hg19)]  # no duplicated rows found: 12,292 => 11,785 SNPs
GenomeInfoDb::seqlevelsStyle(gtex_wb_gr_hg19) <- "NCBI"
saveRDS(gtex_wb_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds")
rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)

### GTEx Brain (all brain tissues) eQTL
brains <- c("Brain_Amygdala.v8.egenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v8.egenes.txt", "Brain_Caudate_basal_ganglia.v8.egenes.txt", "Brain_Cerebellar_Hemisphere.v8.egenes.txt", "Brain_Cerebellum.v8.egenes.txt", "Brain_Cortex.v8.egenes.txt", "Brain_Frontal_Cortex_BA9.v8.egenes.txt", "Brain_Hippocampus.v8.egenes.txt", "Brain_Hypothalamus.v8.egenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v8.egenes.txt", "Brain_Putamen_basal_ganglia.v8.egenes.txt", "Brain_Spinal_cord_cervical_c-1.v8.egenes.txt", "Brain_Substantia_nigra.v8.egenes.txt")
gtex_brains <- data.frame()
for (i in brains) {
  filename <- paste0("/storage/groups/ccm01/datasets/2019_GTEx_v8/GTEx_Analysis_v8_eQTL/", i)
  data <- data.table::fread(filename)
  data <- data[data$qval < 0.05, ]  # filter for significant results
  gtex_brains <- rbind(gtex_brains, data) # merge different brain files together
}
gtex_brains_gr <- GenomicRanges::GRanges(seqnames = gtex_brains$chr, 
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_brains$variant_pos)), 
                                                                   end = as.numeric(as.character(gtex_brains$variant_pos))),
                                         snp_id = gtex_brains$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr) <- "UCSC"
gtex_brains_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_brains_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 89,312 => 88,885
gtex_brains_gr_hg19 <- gtex_brains_gr_hg19[!duplicated(gtex_brains_gr_hg19)]  # no duplicated rows found: 88,885 => 61,935 SNPs
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr_hg19) <- "NCBI"
saveRDS(gtex_brains_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Brains_hg19_eQTL_SNPs.rds")
rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19)

### GTEx Whole Blood splice QTL (sQTLs)
gtex_wb <- data.table::fread("/storage/groups/ccm01/datasets/2019_GTEx_v8/GTEx_Analysis_v8_sQTL/Whole_Blood.v8.sgenes.txt")
gtex_wb <- gtex_wb[gtex_wb$qval < 0.05, ] # filter for FDR < 0.05 (qvalue)
gtex_wb_gr <- GenomicRanges::GRanges(seqnames = gtex_wb$chr,
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_wb$variant_pos)),
                                                               end = as.numeric(as.character(gtex_wb$variant_pos))),
                                     snp_id = gtex_wb$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(gtex_wb_gr) <- "UCSC"  # UCSC needed for the lift over
gtex_wb_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_wb_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 3,013 => 3,004
gtex_wb_gr_hg19 <- gtex_wb_gr_hg19[!duplicated(gtex_wb_gr_hg19)]  # no duplicated rows found: 3,004 => 2,837 SNPs
GenomeInfoDb::seqlevelsStyle(gtex_wb_gr_hg19) <- "NCBI"
saveRDS(gtex_wb_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds")
rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)

### GTEx Brain (all brain tissues) splice QTLs (sQTLs)
brains <- c("Brain_Amygdala.v8.sgenes.txt", "Brain_Anterior_cingulate_cortex_BA24.v8.sgenes.txt", "Brain_Caudate_basal_ganglia.v8.sgenes.txt", "Brain_Cerebellar_Hemisphere.v8.sgenes.txt", "Brain_Cerebellum.v8.sgenes.txt", "Brain_Cortex.v8.sgenes.txt", "Brain_Frontal_Cortex_BA9.v8.sgenes.txt", "Brain_Hippocampus.v8.sgenes.txt", "Brain_Hypothalamus.v8.sgenes.txt", "Brain_Nucleus_accumbens_basal_ganglia.v8.sgenes.txt", "Brain_Putamen_basal_ganglia.v8.sgenes.txt", "Brain_Spinal_cord_cervical_c-1.v8.sgenes.txt", "Brain_Substantia_nigra.v8.sgenes.txt")
gtex_brains <- data.frame()
for (i in brains) {
    filename <- paste0("/storage/groups/ccm01/datasets/2019_GTEx_v8/GTEx_Analysis_v8_sQTL/", i)
    data <- data.table::fread(filename)
    data <- data[data$qval < 0.05, ]  # filter for significant results
    gtex_brains <- rbind(gtex_brains, data) # merge different brain files together
}
gtex_brains_gr <- GenomicRanges::GRanges(seqnames = gtex_brains$chr,
                                         ranges = IRanges::IRanges(start = as.numeric(as.character(gtex_brains$variant_pos)),
                                                                   end = as.numeric(as.character(gtex_brains$variant_pos))),
                                         snp_id = gtex_brains$rs_id_dbSNP151_GRCh38p7)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr) <- "UCSC"
gtex_brains_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_brains_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 20,642 => 20,527
gtex_brains_gr_hg19 <- gtex_brains_gr_hg19[!duplicated(gtex_brains_gr_hg19)]  # no duplicated rows found: 20,527 => 11,804 SNPs
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr_hg19) <- "NCBI"
saveRDS(gtex_brains_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Brains_hg19_sQTL_SNPs.rds")
rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19, hg38_to_hg19)

### VEP: Variant Effect Predictor
#vep <- data.table::fread("WGS_Feature_overlap_collapsed_VEP_short_4torus.MAF01.txt", data.table = F)
#vep$SNP.Position <- as.numeric(sapply(as.character(vep$SNP),function(x)unlist(strsplit(x,"_"))[2]))
#vep$SNP.Chr <- as.character(sapply(as.character(vep$SNP),function(x)unlist(strsplit(x,"_"))[1]))
#vep_melt <- melt(vep, id.vars = c("SNP", "SNP.Chr", "SNP.Position"))
#vepGR <- makeGRangesFromDataFrame(data.frame(chr = vep_melt$SNP.Chr, start = vep_melt$SNP.Position, end = vep_melt$SNP.Position, code = as.character(vep_melt$variable)), keep.extra.columns = TRUE) 
##lift down to hg19
#hg38ToHg19 <- import.chain("hg38ToHg19.over.chain")
#vepGR_hg19 <- unlist(liftOver(vepGR, hg38ToHg19))#10719405
#seqlevelsStyle(vepGR_hg19) <- "NCBI"
#saveRDS(vepGR_hg19, "VEP_hg19_collapsed_MAF01.rds")

