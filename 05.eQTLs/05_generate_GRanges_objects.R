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
library(biomaRt) # to get SNP positions

# function to generate an eQTL GRanges object (with MAF bins)
# @param input Input data of which an GRanges object should be made of (e.g. eQTL hits)
# @param mdata Merging data (e.g. Background SNPs)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObjecteQTL <- function(input, mdata, ofile){
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
    GenomeInfoDb::seqlevelsStyle(input_gr) <- "NCBI"
    # save GRanges object to an rds file
    saveRDS(input_gr, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData_files/", ofile))
    # return finalized GRanges object
    return(input_gr)
}

# function to generate a GRanges object (for all datasets with Chromosome, base pair position and rs SNP ID information)
# @param input Input data of which an GRanges object should be made of (needs CHR, BP, P and SNP as column names)
# @param ofile Name for the rds output file of the generated GRanges object 
GenerateGrangesObject <- function(input, ofile){
  # filter for p-value (P) < 0.05
  input <- input[input$P < 0.05, ]
  # generate GRanges object
  input_gr <- GenomicRanges::GRanges(seqnames = input$CHR,
                                     ranges = IRanges::IRanges(start = as.numeric(as.character(input$BP)),
                                                              end = as.numeric(as.character(input$BP))),
                                     snp_id = input$SNP,
                                     p_value = input$P)  # generate GRanges object only with chr, position, rs SNP ID and p-value
  # remove duplicates
  input_gr <- input_gr[!duplicated(input_gr)] 
  GenomeInfoDb::seqlevelsStyle(input_gr) <- "NCBI"
  print(length(input_gr)) # print length
  # save GRanges object to an rds file
  saveRDS(input_gr, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/", ofile))
  rm(input_gr)
}

### Load BiomaRt SNP informations
snpmart <- useMart(biomart = "ENSEMBL_MART_SNP", 
                   host = "grch37.ensembl.org", 
                   path = "/biomart/martservice", 
                   dataset = "hsapiens_snp")


#######################
### Background SNPs ###
#######################
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
background_gr <- GenerateGrangesObjecteQTL(input = snps, 
                                       mdata = mafs, 
                                       ofile = "background_SNPs.rds")
table(seqnames(background_gr))  # no chr X/Y/M (remove in dataset if included there)
length(background_gr$snp_id) # 6,798,048
rm(snps)

# Generate merging file for eQTL data from background SNPs 
background <- as.data.frame(background_gr)[, c("seqnames", "start", "snp_id", "bin")]
colnames(background) <- c("CHR", "POS", "SNP", "bin")
head(background)


#################
### PMB eQTLs ###
#################
### geSNPs ###
# Gene eQTLs
geSNPs <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
geSNPs_gr <- GenerateGrangesObjecteQTL(input = geSNPs,
                                       mdata = background,
                                       ofile = "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(geSNPs_gr$snp_id) # 52,496 SNPs
rm(geSNPs, geSNPs_gr)

# Clumped gene eQTLs
geSNPs_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_vs_controls_clumped.txt", header = T, stringsAsFactors = F)
geSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = geSNPs_clumped,
                                       mdata = background,
                                       ofile = "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
length(geSNPs_clumped_gr$snp_id) # 1,401 SNPs
rm(geSNPs_clumped, geSNPs_clumped_gr)

# Gene cases only eQTLs
geSNPs_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
geSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = geSNPs_cases,
                                             mdata = background,
                                             ofile = "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds")
length(geSNPs_cases_gr$snp_id) # 18,944 SNPs
rm(geSNPs_cases, geSNPs_cases_gr)

# Clumped gene cases only eQTLs
geSNPs_cases_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_only_clumped.txt", header = T, stringsAsFactors = F)
geSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = geSNPs_cases_clumped,
                                             mdata = background,
                                             ofile = "pmb_gene_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
length(geSNPs_cases_clumped_gr$snp_id) # 485 SNPs
rm(geSNPs_cases_clumped, geSNPs_cases_clumped_gr)

### teSNPs ###
# Transcript eQTLs
teSNPs <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
teSNPs_gr <- GenerateGrangesObjecteQTL(input = teSNPs,
                                             mdata = background,
                                             ofile = "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(teSNPs_gr$snp_id) # 124,888 SNPs
rm(teSNPs, teSNPs_gr)

# Clumped transcript eQTLs
teSNPs_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_vs_controls_clumped.txt", header = T, stringsAsFactors = F)
teSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = teSNPs_clumped,
                                    mdata = background,
                                    ofile = "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
length(teSNPs_clumped_gr$snp_id) # 4,176 SNPs
rm(teSNPs_clumped, teSNPs_clumped_gr)

# Transcript cases only eQTLs
teSNPs_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
teSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = teSNPs_cases,
                                                  mdata = background,
                                                  ofile = "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds")
length(teSNPs_cases_gr$snp_id) # 46,461 SNPs
rm(teSNPs_cases, teSNPs_cases_gr)

# Clumped transcript cases only eQTLs
teSNPs_cases_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_only_clumped.txt", header = T, stringsAsFactors = F)
teSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = teSNPs_cases_clumped,
                                          mdata = background,
                                          ofile = "pmb_transcript_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
length(teSNPs_cases_clumped_gr$snp_id) # 1,481 SNPs
rm(teSNPs_cases_clumped, teSNPs_cases_clumped_gr)

### eeSNPs ###
# Exon eQTLs
eeSNPs <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eeSNPs_gr <- GenerateGrangesObjecteQTL(input = eeSNPs,
                                       mdata = background,
                                       ofile = "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds")
length(eeSNPs_gr$snp_id) # 314,185 SNPs
rm(eeSNPs, eeSNPs_gr)

# Clumped exon eQTLs
eeSNPs_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_vs_controls_clumped.txt", header = T, stringsAsFactors = F)
eeSNPs_clumped_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_clumped,
                                    mdata = background,
                                    ofile = "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs_clumped.rds")
length(eeSNPs_clumped_gr$snp_id) # 10,135 SNPs
rm(eeSNPs_clumped, eeSNPs_clumped_gr)

# Exon cases only eQTLs
eeSNPs_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)
eeSNPs_cases_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_cases,
                                             mdata = background,
                                             ofile = "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds")
length(eeSNPs_cases_gr$snp_id) # 138,662 SNPs
rm(eeSNPs_cases, eeSNPs_cases_gr)

# Clumped exon cases only eQTLs
eeSNPs_cases_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_only_clumped.txt", header = T, stringsAsFactors = F)
eeSNPs_cases_clumped_gr <- GenerateGrangesObjecteQTL(input = eeSNPs_cases_clumped,
                                          mdata = background,
                                          ofile = "pmb_exon_hits_cases_only_cis_eQTL_SNPs_clumped.rds")
length(eeSNPs_cases_clumped_gr$snp_id) # 4,367 SNPs
rm(eeSNPs_cases_clumped, eeSNPs_cases_clumped_gr)


###########
### CMC ###
###########
# CMC (Fromer et al.) gene eQTLs: geSNPs
cmc_gene <- data.table::fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_eQTL-adjustedSVA.txt")
cmc_gene <- cmc_gene[cmc_gene$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc_gene_gr <- GenomicRanges::GRanges(seqnames = cmc_gene$SNP_location, 
                                      snp_id = cmc_gene$SNP)  # generate GRanges object 
cmc_gene_gr <- cmc_gene_gr[!duplicated(cmc_gene_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
GenomeInfoDb::seqlevelsStyle(cmc_gene_gr) <- "NCBI"
seqlevels(cmc_gene_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cmc_gene_gr), c("X", "Y", "M"))
saveRDS(cmc_gene_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/CMC_adjustedSVA_eQTL_SNPs.rds")
rm(cmc_gene, cmc_gene_gr)

# CMC (Fromer et al.) isoform eQTLs: ieSNPs
cmc_iso <- data.table::fread("/storage/groups/ccm01/datasets/2019_CMC_Fromer_et_al/eQTL/CMC_MSSM-Penn-Pitt_DLPFC_Analysis_isoQTL-adjustedSVA.txt")
cmc_iso <- cmc_iso[cmc_iso$FDR < 0.05, ]  # filter for FDR <= 0.05
cmc_iso_gr <- GenomicRanges::GRanges(seqnames = cmc_iso$SNP_location,
                                     snp_id = cmc_iso$SNP)  # generate GRanges object 
cmc_iso_gr <- cmc_iso_gr[!duplicated(cmc_iso_gr)]  # remove duplicated rows: 2,199,772 => 1,342,767 SNPs
GenomeInfoDb::seqlevelsStyle(cmc_iso_gr) <- "NCBI"
seqlevels(cmc_iso_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cmc_iso_gr), c("X", "Y", "M"))
saveRDS(cmc_iso_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/CMC_adjustedSVA_isoQTL_SNPs.rds")
rm(cmc_iso, cmc_iso_gr)


############
### GTEx ###
############
### eQTL SNPs (eSNPs) ###
# GTEx v8 is based on the new human genome GRCh38, as we are working on GRCH37/hg19 we need to downlift to this genome
# Load hg38 to hg19 chain (chain downloaded from UCSC: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz)
hg38_to_hg19 <- rtracklayer::import.chain("/storage/groups/ccm01/datasets/2019_GTEx_v8/hg38ToHg19.over.chain")

# GTEx Whole Blood
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
seqlevels(gtex_wb_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_wb_gr_hg19), c("X", "Y", "M"))
saveRDS(gtex_wb_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds")
rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)

# GTEx Brain (all brain tissues)
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
                                         snp_id = gtex_brains$rs_id_dbSNP151_GRCh38p7,
                                         state = gtex_brains)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr) <- "UCSC"
gtex_brains_gr_hg19 <- unlist(rtracklayer::liftOver(gtex_brains_gr, hg38_to_hg19))  # lift over from GRCh38 to GRCh37 (hg19): 89,312 => 88,885
gtex_brains_gr_hg19 <- gtex_brains_gr_hg19[!duplicated(gtex_brains_gr_hg19)]  # no duplicated rows found: 88,885 => 61,935 SNPs
GenomeInfoDb::seqlevelsStyle(gtex_brains_gr_hg19) <- "NCBI"
seqlevels(gtex_brains_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_brains_gr_hg19), c("X", "Y", "M"))
saveRDS(gtex_brains_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Brains_hg19_eQTL_SNPs.rds")
rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19)


### splice QTL (sQTL) SNPs (seSNPs) ###
# GTEx Whole Blood
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
seqlevels(gtex_wb_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_wb_gr_hg19), c("X", "Y", "M"))
saveRDS(gtex_wb_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds")
rm(gtex_wb, gtex_wb_gr, gtex_wb_gr_hg19)

# GTEx Brain (all brain tissues)
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
seqlevels(gtex_brains_gr_hg19, pruning.mode = "coarse") <- setdiff(seqlevels(gtex_brains_gr_hg19), c("X", "Y", "M"))
saveRDS(gtex_brains_gr_hg19, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GTEx_v8_Brains_hg19_sQTL_SNPs.rds")
rm(brains, filename, data, i, gtex_brains, gtex_brains_gr, gtex_brains_gr_hg19, hg38_to_hg19)


####################
### Cortex meQTL ###
####################
cortex_meqtl <- fread("/storage/groups/ccm01/datasets/2019_cortex_meQTL/GSE74193_jaffe_onlineTable2_meQTLs_allPairs.csv")
cortex_meqtl <- cortex_meqtl[cortex_meqtl$p < 0.05, ]
cortex_meqtl_gr <- GRanges(seqnames = cortex_meqtl$snpChr, 
                           ranges = IRanges(start = as.numeric(as.character(cortex_meqtl$snpPos)), 
                                            end = as.numeric(as.character(cortex_meqtl$snpPos))),
                           snp_id = cortex_meqtl$snpRsNum)  # generate GRanges object only with gene ID, name and location
cortex_meqtl_gr <- cortex_meqtl_gr[!duplicated(cortex_meqtl_gr)]  # remove duplicated rows: 8,663 genes (no duplicates)
GenomeInfoDb::seqlevelsStyle(cortex_meqtl_gr) <- "NCBI"
seqlevels(cortex_meqtl_gr, pruning.mode = "coarse") <- setdiff(seqlevels(cortex_meqtl_gr), c("X", "Y", "M"))
saveRDS(cortex_meqtl_gr, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/Cortex_meQTL_SNPs.rds")
rm(cortex_meqtl, cortex_meqtl_gr)


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


###############################
### GWAS Summary Statistics ###
###############################
### AD-2019 ###
ad <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/AD-2019/AD_sumstats_Jansenetal_2019sept.txt", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(ad, ofile = "GWAS_sumstats_AD-2019_p0.05.rds")
rm(ad)

### ADHD-2017 ###
adhd <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/ADHD-2017/adhd_jul2017", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(adhd, ofile = "GWAS_sumstats_ADHD-2017_p0.05.rds")
rm(adhd)

### ALCDEP-2018 ###
alcdep <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/ALCDEP-2018/pgc_alcdep.discovery.aug2018_release.txt", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
alcdep <- alcdep[alcdep$P < 0.05, ]  # filter for p <= 0.05 (as otherwise it takes too long to change the SNP IDs)
alcdep$SNP <- sapply(alcdep$SNP, function(x) if(startsWith(x, "rs")) strsplit(x, ":")[[1]][1] else x) # remove CHR, BP and Alleles from the rs SNP ID
GenerateGrangesObject(alcdep, ofile = "GWAS_sumstats_ALCDEP-2018_p0.05.rds")
rm(alcdep)

### AN-2019 ###
an <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/AN-2019/pgcAN2.2019-07.vcf.tsv", header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"), skip = 70)
colnames(an) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(an, ofile = "GWAS_sumstats_AN-2019_p0.05.rds")
rm(an)

### ASD-2017 ###
asd <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/ASD-2017/iPSYCH-PGC_ASD_Nov2017", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(asd, ofile = "GWAS_sumstats_ASD-2017_p0.05.rds")
rm(asd)

### AUDIT-2018 ###
audit <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/AUDIT-2018/AUDIT_UKB_2018_AJP.txt", header = T, stringsAsFactors = F, select = c("chr", "rsid", "p_T"))
audit <- audit[audit$p_T < 0.05, ]  # filter for p <= 0.05 (as otherwise it takes too long to change the SNP IDs)
# split into IDs without rs and with rs SNP ID (to get positions)
audit_chr_id <- audit[!(startsWith(audit$rsid, "rs")), ] 
audit_chr_id$BP <- sapply(audit_chr_id$rsid, function(x) strsplit(strsplit(x, ":")[[1]][2], "_")[[1]][1])
audit_chr_id <- audit_chr_id[, c("chr", "BP", "rsid", "p_T")]
colnames(audit_chr_id) <- c("CHR", "BP", "SNP", "P")
audit_rs_id <- audit[startsWith(audit$rsid, "rs"), ] 
# get positions for the rs IDs
rsids <- getBM(attributes = c("refsnp_id", "chrom_end"), 
               filters = "snp_filter",
               values = audit_rs_id$rsid,
               mart = snpmart,
               uniqueRows = TRUE)
# merge with AUDIT rs SNP IDs
audit_rs_id <- merge(audit_rs_id, rsids, by.x = "rsid", by.y = "refsnp_id", all.x = T)
audit_rs_id <- audit_rs_id[, c("chr", "chrom_end", "rsid", "p_T")]
colnames(audit_rs_id) <- c("CHR", "BP", "SNP", "P")
audit_rs_id$BP <- as.character(audit_rs_id$BP)
audit <- rbind(audit_chr_id, audit_rs_id)
audit <- audit[complete.cases(audit)] # remove SNPs without a position
GenerateGrangesObject(audit, ofile = "GWAS_sumstats_AUDIT-2018_p0.05.rds")
rm(audit, audit_chr_id, audit_rs_id, rsids)

### BIP-2018 ###
bip <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/BIP-2018/daner_PGC_BIP32b_mds7a_0416a", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(bip, ofile = "GWAS_sumstats_BIP-2018_p0.05.rds")
rm(bip)

### BMI-2015 ###
bmi <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/BMI-2015/ENGAGE1000G_BMI.txt", header = T, stringsAsFactors = F, select = c("chromosome", "position", "SNP", "P"))
colnames(bmi) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(bmi, ofile = "GWAS_sumstats_BMI-2015_p0.05.rds")
rm(bmi)

### CAD-2017 ###
# needs downlift from hg38 to hg19!
cad <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/CAD-2017/UKBB.GWAS1KG.EXOME.CAD.SOFT.META.PublicRelease.300517.txt", header = T, stringsAsFactors = F, select = c("chr", "bp_hg19", "snptestid", "p-value_gc"))
cad <- cad[cad$`p-value_gc` < 0.05, ] # filter for p-value (pval) <= 0.05
cad <- cad[complete.cases(cad)] # remove rows with empty chromosome and bp position
# generate GRanges object
cad_gr <- GenomicRanges::GRanges(seqnames = cad$chr,
                                 ranges = IRanges::IRanges(start = as.numeric(as.character(cad$bp_hg19)),
                                                           end = as.numeric(as.character(cad$bp_hg19))),
                                 snp_id = cad$snptestid)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(cad_gr) <- "UCSC"
# load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
hg38_to_hg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg38ToHg19.over.chain")
# uplift from hg18 to hg19
cad_gr_hg19 <- unlist(liftOver(cad_gr, hg38_to_hg19))
# remove duplicates, change to NCBI format and print length
cad_gr_hg19 <- cad_gr_hg19[!duplicated(cad_gr_hg19)] 
GenomeInfoDb::seqlevelsStyle(cad_gr_hg19) <- "NCBI"
print(length(cad_gr_hg19)) # print length
# save GRanges object to an rds file
saveRDS(cad_gr_hg19, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GWAS_sumstats_CAD-2017_p0.05.rds")
rm(cad, cad_gr, hg38_to_hg19, cad_gr_hg19)

### Cannabis-2018 ###
cannabis <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Cannabis-2018/Cannabis_ICC_UKB_het.txt", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
cannabis <- cannabis[complete.cases(cannabis)] # remove rows with missing information 
GenerateGrangesObject(cannabis, ofile = "GWAS_sumstats_Cannabis-2018_p0.05.rds")
rm(cannabis)

### CDG1-2013 ###
# needs uplift from hg18 to hg19!
cdg1 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/CDG1-2013/pgc.cross.full.2013-03.txt", header = T, stringsAsFactors = F, select = c("hg18chr", "bp", "snpid", "pval"))
cdg1 <- cdg1[cdg1$pval < 0.05, ] # filter for p-value (pval) <= 0.05
# generate GRanges object
cdg1_gr <- GenomicRanges::GRanges(seqnames = cdg1$hg18chr,
                                  ranges = IRanges::IRanges(start = as.numeric(as.character(cdg1$bp)),
                                                            end = as.numeric(as.character(cdg1$bp))),
                                  snp_id = cdg1$snpid)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(cdg1_gr) <- "UCSC"
# load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
hg18_to_hg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg18ToHg19.over.chain")
# uplift from hg18 to hg19
cdg1_gr_hg19 <- unlist(liftOver(cdg1_gr, hg18_to_hg19))
# remove duplicates, change to NCBI format and print length
cdg1_gr_hg19 <- cdg1_gr_hg19[!duplicated(cdg1_gr_hg19)] 
GenomeInfoDb::seqlevelsStyle(cdg1_gr_hg19) <- "NCBI"
print(length(cdg1_gr_hg19)) # print length
# save GRanges object to an rds file
saveRDS(cdg1_gr_hg19, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GWAS_sumstats_CDG1-2013_p0.05.rds")
rm(cdg1, cdg1_gr, hg18_to_hg19, cdg1_gr_hg19)

### CDG2-2019 ###
cdg2 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/CDG2-2019/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt", header = T, stringsAsFactors = F, select = c("CHROM", "POS", "ID", "PVAL"))
colnames(cdg2) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(cdg2, ofile = "GWAS_sumstats_CDG2-2019_p0.05.rds")
rm(cdg2)

### EA-2018 ###
ea <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/EA-2018/GWAS_EA_excl23andMe.txt", header = T, stringsAsFactors = F, select = c("CHR", "POS", "MarkerName", "Pval"))
colnames(ea) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(ea, ofile = "GWAS_sumstats_EA-2018_p0.05.rds")
rm(ea)

### ED-2017 ###
ed <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/ED-2017/pgc.ed.freeze1.summarystatistics.July2017.txt", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(ed, ofile = "GWAS_sumstats_ED-2017_p0.05.rds")
rm(ed)

### Height-2014 ###
# needs uplift from hg18 to hg19!
height <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Height-2014/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt", header = T, stringsAsFactors = F, select = c("MarkerName", "p"))
height <- height[height$p < 0.05, ] # filter for p-value (p) <= 0.05
# load HapMap to get chromosome and position for all SNP IDs
hapmap <- data.table::fread(file = "/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Height-2014/hapmap36.map", header = T, stringsAsFactors = F)
hapmap$Zero <- NULL
# merge HapMap and Height 
height <- merge(height, hapmap, by.x = "MarkerName", by.y = "SNPID")
# generate GRanges object
height_gr <- GenomicRanges::GRanges(seqnames = height$Chromosome,
                                    ranges = IRanges::IRanges(start = as.numeric(as.character(height$Position)),
                                                             end = as.numeric(as.character(height$Position))),
                                    snp_id = height$MarkerName)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(height_gr) <- "UCSC"
# load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
hg18_to_hg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg18ToHg19.over.chain")
# uplift from hg18 to hg19
height_gr_hg19 <- unlist(liftOver(height_gr, hg18_to_hg19))
# remove duplicates, change to NCBI format and print length
height_gr_hg19 <- height_gr_hg19[!duplicated(height_gr_hg19)] 
GenomeInfoDb::seqlevelsStyle(height_gr_hg19) <- "NCBI"
print(length(height_gr_hg19)) # print length
# save GRanges object to an rds file
saveRDS(height_gr_hg19, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GWAS_sumstats_Height-2014_p0.05.rds")
rm(height, hapmap, height_gr, hg18_to_hg19, height_gr_hg19)

### Intelligence-2018 ###
intelligence <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Intelligence-2018/SavageJansen_2018_intelligence_metaanalysis.txt", header = T, stringsAsFactors = F, select = c("CHR", "POS", "SNP", "P"))
colnames(intelligence) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(intelligence, ofile = "GWAS_sumstats_Intelligence-2018_p0.05.rds")
rm(intelligence)

### MDD1-2012 ###
# needs uplift from hg18 to hg19!
mdd1 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/MDD1-2012/pgc.mdd.full.2012-04.txt", header = T, stringsAsFactors = F, select = c("hg18chr", "bp", "snpid", "pval"))
mdd1 <- mdd1[mdd1$pval < 0.05, ] # filter for p-value (pval) <= 0.05
# generate GRanges object
mdd1_gr <- GenomicRanges::GRanges(seqnames = mdd1$hg18chr,
                                    ranges = IRanges::IRanges(start = as.numeric(as.character(mdd1$bp)),
                                                              end = as.numeric(as.character(mdd1$bp))),
                                    snp_id = mdd1$snpid)  # generate GRanges object only with chr, position and rs SNP ID
GenomeInfoDb::seqlevelsStyle(mdd1_gr) <- "UCSC"
# load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
hg18_to_hg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg18ToHg19.over.chain")
# uplift from hg18 to hg19
mdd1_gr_hg19 <- unlist(liftOver(mdd1_gr, hg18_to_hg19))
# remove duplicates, change to NCBI format and print length
mdd1_gr_hg19 <- mdd1_gr_hg19[!duplicated(mdd1_gr_hg19)] 
GenomeInfoDb::seqlevelsStyle(mdd1_gr_hg19) <- "NCBI"
print(length(mdd1_gr_hg19)) # print length
# save GRanges object to an rds file
saveRDS(mdd1_gr_hg19, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/GWAS_sumstats_MDD1-2012_p0.05.rds")
rm(mdd1, mdd1_gr, hg18_to_hg19, mdd1_gr_hg19)

### MDD2-2018 ###
mdd2 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/MDD2-2018/daner_pgc_mdd_meta_align_170228_filtered_incl23andMe", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(mdd2, ofile = "GWAS_sumstats_MDD2-2018_p0.05.rds")
rm(mdd2)

### MDD3-2019 ###
mdd3 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/MDD3-2019/PGC_UKB_depression_genome-wide.txt", header = T, stringsAsFactors = F, select = c("MarkerName", "P"))
mdd3 <- mdd3[mdd3$P < 0.05, ]  # filter for p <= 0.05 (as otherwise it takes too long to get the positions)
# get positions for the rs IDs
rsids <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_end"), 
               filters = "snp_filter",
               values = mdd3$MarkerName,
               mart = snpmart,
               uniqueRows = TRUE)
# merge with MDD3 rs SNP IDs
mdd3 <- merge(mdd3, rsids, by.x = "MarkerName", by.y = "refsnp_id", all.x = T)
# reorder columns, change column names and remove all rows with NA's
mdd3 <- mdd3[, c("chr_name", "chrom_end", "MarkerName", "P")]
colnames(mdd3) <- c("CHR", "BP", "SNP", "P")
mdd3 <- mdd3[complete.cases(mdd3)]
# create GRanges object and save into .rds file
GenerateGrangesObject(mdd3, ofile = "GWAS_sumstats_MDD3-2019_p0.05.rds")
rm(mdd3, rsids)

### OCD-2017 ###
ocd <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/OCD-2017/ocd_aug2017", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(ocd, ofile = "GWAS_sumstats_OCD-2017_p0.05.rds")
rm(ocd)

### PTSD-2019 ###
ptsd <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/PTSD-2019/pts_all_freeze2_overall.results", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(ptsd, ofile = "GWAS_sumstats_PTSD-2019_p0.05.rds")
rm(ptsd)

### SCZ2-2014 ###
scz2 <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/SCZ2-2014/scz2.snp.results.txt", header = T, stringsAsFactors = F, select = c("hg19chrc", "bp", "snpid", "p"))
colnames(scz2) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(scz2, ofile = "GWAS_sumstats_SCZ2-2014_p0.05.rds")
rm(scz2)

### T2D-2017 ###
t2d <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/T2D-2017/METAANALYSIS_DIAGRAM_SE1.txt", header = T, stringsAsFactors = F, select = c("Chr:Position", "P-value"))
t2d <- t2d[t2d$`P-value` < 0.05, ]  # filter for p <= 0.05 (as otherwise it takes too long to get the chr/positions)
t2d$"CHR" <- sapply(t2d$`Chr:Position`, function(x) strsplit(x, ":")[[1]][1]) # get chr
t2d$"BP" <- sapply(t2d$`Chr:Position`, function(x) strsplit(x, ":")[[1]][2]) # get bp position
t2d <- t2d[, c("CHR", "BP", "Chr:Position", "P-value")]
colnames(t2d) <- c("CHR", "BP", "SNP", "P")
GenerateGrangesObject(t2d, ofile = "GWAS_sumstats_T2D-2017_p0.05.rds")
rm(t2d)

### TS-2018 ###
ts <- data.table::fread("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/TS-2018/TS_Oct2018", header = T, stringsAsFactors = F, select = c("CHR", "BP", "SNP", "P"))
GenerateGrangesObject(ts, ofile = "GWAS_sumstats_TS-2018_p0.05.rds")
rm(ts)
