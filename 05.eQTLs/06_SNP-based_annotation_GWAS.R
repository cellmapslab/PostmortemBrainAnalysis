# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

# load packages
library(data.table)
library(GenomicRanges)
library(tidyr)
library(arules)
library(S4Vectors)
library(reshape2)
library(rtracklayer)
library(foreach)
library(doSNOW)

# function, to calculate the SNP annotation to a public dataset, but with respect to the MAF bins
# @param own Own cis-eQTL data set as GRanges object
# @param public A public cis-eQTL data set as GRanges object
# @param nperm Number of permutations (usually a 1000 times)
# @param background All background genes coming from the gene expression (here: 14,048 unqiue genes)
# @return vector including the number of overlapping SNPs (1 = overlap), the resampling p-value (2 = p_value) and the resampling odds ratio (3 = or)
AnnotatePublicData <- function(own, public, nperm, background){ 
  overlap <- sum(overlapsAny(own, public))  # number of overlaps between the two datasets
  background_bins <- lapply(1:11, function(x) background[background$bin == x])  # get all background SNPs in MAF bin 1 to 11
  own_bin_lengths <- sapply(1:11, function(x) length(own[own$bin == x]))  # get length of all 11 MAF bins
  # resample nperm times
  resampling <- lapply(1:nperm, function(x)
    sum(sapply(1:11, function(x)
      sum(
        IRanges::overlapsAny(sample(background_bins[[x]], own_bin_lengths[x]), public)
      ))))
  # IF 0 => actually smaller than 1e^(-5/6)
  p_value <- if((sum(resampling >= overlap) / nperm) == 0){ paste("<=", 1/nperm, sep = " ") } 
  else {sum(resampling >= overlap) / nperm}
  or <- overlap/mean(unlist(resampling))
  return(c(as.character(overlap), as.character(p_value), as.character(or))) # return the resampling overlapping SNPs, p-value and odds ratio
}

# path, where all datasets are stored (as .rds files)
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/06_RData/"


###########
### PMB ###
###########
### Load eQTL results ###
eqtls_gene <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_transcript <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_exon <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))
eqtls_gene_cases <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_transcript_cases <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))
eqtls_exon_cases <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))

### Load Background data ###
background <- readRDS(paste0(ipath,"background_SNPs.rds"))
length(background@ranges) #[1] 6,798,048


###########
### CMC ###
###########
# CMC (Fromer et al.) gene eQTLs and isoQTLs
cmc_eqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_eQTL_SNPs.rds"))
cmc_isoqtl <- readRDS(paste0(ipath, "CMC_adjustedSVA_isoQTL_SNPs.rds"))


############
### GTEx ###
############
### GTEx Whole Blood eQTLs and sQTLs
gtex_wb_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_eQTL_SNPs.rds"))
gtex_wb_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Whole_Blood_hg19_sQTL_SNPs.rds"))

### GTEx Brain (all brain tissues) eQTLs and sQTLs
gtex_brains_eqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_eQTL_SNPs.rds"))
gtex_brains_sqtl <- readRDS(paste0(ipath, "GTEx_v8_Brains_hg19_sQTL_SNPs.rds"))


#############
### GWAS ###
#############
### Load public data: GWAS summary statistics ###
# datasets are already filtered for p-value <= 0.05
ad <- readRDS(paste0(ipath, "GWAS_sumstats_AD-2019_p0.05.rds"))  # GWAS Alzheimer's Disease (AD-2019)
adhd <- readRDS(paste0(ipath, "GWAS_sumstats_ADHD-2017_p0.05.rds"))  # GWAS Attention-Deficit Hyperactivity Disorder (ADHD-2017)
alcdep <- readRDS(paste0(ipath, "GWAS_sumstats_ALCDEP-2018_p0.05.rds"))  # GWAS Alcohol Dependence (ALCDEP-2018) <- similar to AUDIT-2018
an <- readRDS(paste0(ipath, "GWAS_sumstats_AN-2019_p0.05.rds"))  # GWAS Anorexia Nervosa (AN-2019) <- similar to ED-2017
asd <- readRDS(paste0(ipath, "GWAS_sumstats_ASD-2017_p0.05.rds"))  # GWAS Autism Spectrum Disorder (ASD-2017)
audit <- readRDS(paste0(ipath, "GWAS_sumstats_AUDIT-2018_p0.05.rds"))  # GWAS Alcohol Use Disorders Identification Test (AUDIT-2018) <- similar to ALCDEP-2018
bip <- readRDS(paste0(ipath, "GWAS_sumstats_BIP-2018_p0.05.rds"))  # GWAS Bipolar Disorder (BIP-2018)
bmi <- readRDS(paste0(ipath, "GWAS_sumstats_BMI-2015_p0.05.rds"))  # GWAS Body Mass Index / Weight (BMI-2015)
cad <- readRDS(paste0(ipath, "GWAS_sumstats_CAD-2017_p0.05.rds"))  # GWAS Coronary Artery Disease (CAD-2017)
cannabis <- readRDS(paste0(ipath, "GWAS_sumstats_Cannabis-2018_p0.05.rds"))  # GWAS Cannabis Use (Cannabis-2018)
cdg1 <- readRDS(paste0(ipath, "GWAS_sumstats_CDG1-2013_p0.05.rds"))  # GWAS Cross-Disorder Group (CDG1-2013)
cdg2 <- readRDS(paste0(ipath, "GWAS_sumstats_CDG2-2019_p0.05.rds"))  # GWAS Cross-Disorder Group (CDG2-2019)
ea <- readRDS(paste0(ipath, "GWAS_sumstats_EA-2018_p0.05.rds"))  # GWAS Educational Attainment (EA-2018) <- similar to Intelligence-2018
ed <- readRDS(paste0(ipath, "GWAS_sumstats_ED-2017_p0.05.rds"))  # GWAS Eating Disorder (ED-2017) <- similar to AN-2019
height <- readRDS(paste0(ipath, "GWAS_sumstats_Height-2014_p0.05.rds"))  # GWAS Height (Height-2014)
intelligence <- readRDS(paste0(ipath, "GWAS_sumstats_Intelligence-2018_p0.05.rds"))  # GWAS Intelligence (Intelligence-2018) <- similar to EA-2018
mdd1 <- readRDS(paste0(ipath, "GWAS_sumstats_MDD1-2012_p0.05.rds"))  # GWAS Major Depression Disorder  (MDD-2012)
mdd2 <- readRDS(paste0(ipath, "GWAS_sumstats_MDD2-2018_p0.05.rds"))  # GWAS Major Depression Disorder  (MDD-2018)
mdd3 <- readRDS(paste0(ipath, "GWAS_sumstats_MDD3-2019_p0.05.rds"))  # GWAS Major Depression Disorder (MDD-2019)
ocd <- readRDS(paste0(ipath, "GWAS_sumstats_OCD-2017_p0.05.rds"))  # GWAS Obsessive Compulsive Disorder (OCD-2017)
ptsd <- readRDS(paste0(ipath, "GWAS_sumstats_PTSD-2019_p0.05.rds"))  # GWAS Posttraumatic Stress Disorder (PTSD-2019)
scz2 <- readRDS(paste0(ipath, "GWAS_sumstats_SCZ2-2014_p0.05.rds"))  # GWAS Schizophrenia (SCZ2-2014)
t2d <- readRDS(paste0(ipath, "GWAS_sumstats_T2D-2017_p0.05.rds"))  # GWAS Type 2 Diabetes (T2D-2017)
ts <- readRDS(paste0(ipath, "GWAS_sumstats_TS-2018_p0.05.rds"))  # GWAS Tourette Syndrom (TS-2018)

# filter GWAS datasets for genome-wide (gw) hits (p-value <= 5e-8)  
ad_gw <- ad[ad$p_value <= 5e-8, ]
adhd_gw <- adhd[adhd$p_value <= 5e-8, ]
alcdep_gw <- alcdep[alcdep$p_value <= 5e-8, ]
an_gw <- an[an$p_value <= 5e-8, ]
asd_gw <- asd[asd$p_value <= 5e-8, ]
audit_gw <- audit[audit$p_value <= 5e-8, ]
bip_gw <- bip[bip$p_value <= 5e-8, ]
bmi_gw <- bmi[bmi$p_value <= 5e-8, ]
cad_gw <- cad[cad$p_value <= 5e-8, ]
cannabis_gw <- cannabis[cannabis$p_value <= 5e-8, ]
cdg1_gw <- cdg1[cdg1$p_value <= 5e-8, ]
cdg2_gw <- cdg2[cdg2$p_value <= 5e-8, ]
ea_gw <- ea[ea$p_value <= 5e-8, ]
ed_gw <- ed[ed$p_value <= 5e-8, ]
height_gw <- height[height$p_value <= 5e-8, ]
intelligence_gw <- intelligence[intelligence$p_value <= 5e-8, ]
mdd1_gw <- mdd1[mdd1$p_value <= 5e-8, ]
mdd2_gw <- mdd2[mdd2$p_value <= 5e-8, ]
mdd3_gw <- mdd3[mdd3$p_value <= 5e-8, ]
ocd_gw <- ocd[ocd$p_value <= 5e-8, ]
ptsd_gw <- ptsd[ptsd$p_value <= 5e-8, ]
scz2_gw <- scz2[scz2$p_value <= 5e-8, ]
t2d_gw <- t2d[t2d$p_value <= 5e-8, ]
ts_gw <- ts[ts$p_value <= 5e-8, ]

# save names of all GWAS (nominal p < 0.05) in the order they should be running
gwas <- GRangesList("GWAS_SCZ-2014" = scz2, "GWAS_BIP-2018" = bip, "GWAS_MDD1-2012" = mdd1, "GWAS_MDD2-2018" = mdd2, "GWAS_MDD3-2019" = mdd3, "GWAS_CDG1-2013" = cdg1, "GWAS_CDG2-2019" = cdg2, "GWAS_ADHD-2017" = adhd, "GWAS_ASD-2017" = asd, "GWAS_EA-2018" = ea, "GWAS_BMI-2015" = bmi, "GWAS_AD-2019" = ad, "GWAS_ALCDEP-2018" = alcdep, "GWAS_AN-2019" = an, "GWAS_AUDIT-2018" = audit, "GWAS_CAD_2017" = cad, "GWAS_Cannabis-2018" = cannabis, "GWAS_ED-2017" = ed, "GWAS_Height-2014" = height, "GWAS_Intelligence-2018" = intelligence, "GWAS_OCD-2017" = ocd, "GWAS_PTSD-2019" = ptsd, "GWAS_T2D-2017" = t2d, "GWAS_TS-2018" = ts)
rm(ad, adhd, alcdep, an, asd, audit, bip, bmi, cad, cannabis, cdg1, cdg2, ea, ed, height, intelligence, mdd1, mdd2, mdd3, ocd, ptsd, scz2, t2d, ts)
# save names of all GWAS (genome-wide filtered: p < 5e-8) in the order they should be running
gwas_gw <- GRangesList("GWAS_SCZ-2014" = scz2_gw, "GWAS_BIP-2018" = bip_gw, "GWAS_MDD1-2012" = mdd1_gw, "GWAS_MDD2-2018" = mdd2_gw, "GWAS_MDD3-2019" = mdd3_gw, "GWAS_CDG1-2013" = cdg1_gw, "GWAS_CDG2-2019" = cdg2_gw, "GWAS_ADHD-2017" = adhd_gw, "GWAS_ASD-2017" = asd_gw, "GWAS_EA-2018" = ea_gw, "GWAS_BMI-2015" = bmi_gw, "GWAS_AD-2019" = ad_gw, "GWAS_ALCDEP-2018" = alcdep_gw, "GWAS_AN-2019" = an_gw, "GWAS_AUDIT-2018" = audit_gw, "GWAS_CAD_2017" = cad_gw, "GWAS_Cannabis-2018" = cannabis_gw, "GWAS_ED-2017" = ed_gw, "GWAS_Height-2014" = height_gw, "GWAS_Intelligence-2018" = intelligence_gw, "GWAS_OCD-2017" = ocd_gw, "GWAS_PTSD-2019" = ptsd_gw, "GWAS_T2D-2017" = t2d_gw, "GWAS_TS-2018" = ts_gw)
rm(ad_gw, adhd_gw, alcdep_gw, an_gw, asd_gw, audit_gw, bip_gw, bmi_gw, cad_gw, cannabis_gw, cdg1_gw, cdg2_gw, ea_gw, ed_gw, height_gw, intelligence_gw, mdd1_gw, mdd2_gw, mdd3_gw, ocd_gw, ptsd_gw, scz2_gw, t2d_gw, ts_gw)


### Get overlap of all GWAS (p < 0.05) to CMC, GTEx wb and brain
results <- data.frame()
for (i in seq_along(gwas)) {
  gwas_selected <- gwas[[i]]
  print(names(gwas)[[i]])
  results[names(gwas)[[i]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, gwas_selected))
  results[names(gwas)[[i]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, gwas_selected))
  results[names(gwas)[[i]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, gwas_selected))
  results[names(gwas)[[i]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, gwas_selected))
  results[names(gwas)[[i]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, gwas_selected))
  results[names(gwas)[[i]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, gwas_selected))
}
### write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/CMC_GTEx_SNP_overlap_to_all_GWAS_p0.05.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
rm(results)

### Get overlap of all GWAS (p < 0.05) to CMC, GTEx wb and brain
results <- data.frame()
for (i in seq_along(gwas_gw)) {
  gwas_selected <- gwas_gw[[i]]
  print(names(gwas_gw)[[i]])
  results[names(gwas_gw)[[i]], "CMC_eQTL"] <- sum(overlapsAny(cmc_eqtl, gwas_selected))
  results[names(gwas_gw)[[i]], "CMC_isoQTL"] <- sum(overlapsAny(cmc_isoqtl, gwas_selected))
  results[names(gwas_gw)[[i]], "GTEx_eQTL_wb"] <- sum(overlapsAny(gtex_wb_eqtl, gwas_selected))
  results[names(gwas_gw)[[i]], "GTEx_sQTL_wb"] <- sum(overlapsAny(gtex_wb_sqtl, gwas_selected))
  results[names(gwas_gw)[[i]], "GTEx_eQTL_brain"] <- sum(overlapsAny(gtex_brains_eqtl, gwas_selected))
  results[names(gwas_gw)[[i]], "GTEx_sQTL_brain"] <- sum(overlapsAny(gtex_brains_sqtl, gwas_selected))
}
### write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/CMC_GTEx_SNP_overlap_to_all_GWAS_p5e-8.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
rm(results)

### Annotation to genome-wide GWAS hits ###
# First, set parameters
perm <- 10000

# 2.4.1 annotation to GWAS summary statistics (parallelized)
cl <- makeCluster(12)
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above

results2 <- foreach (i = seq_along(gwas_gw), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                       results <- data.frame()
                       overlap <- paste0(names(gwas_gw)[i], "_SNP_overlap")
                       pval <- paste0(names(gwas_gw)[i], "_SNP_pVal")
                       or <- paste0(names(gwas_gw)[i], "_SNP_or")
                       gwas_selected <- gwas_gw[[i]]
                       results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases, gwas_selected, perm, background)
                       write.table(results, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas_gw)[i],"_p5e-8.tsv"), row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
                       results
                     }
stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p5e-8.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
rm(results, results2)


### Annotation to nominal GWAS hits ###
# First, set parameters
perm <- 10000

# 2.4.1 annotation to GWAS summary statistics (parallelized)
cl <- makeCluster(12)
registerDoSNOW(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())  # load packages from above

results2 <- foreach (i = seq_along(gwas), .combine = rbind,
                     .packages = c("GenomicRanges", "S4Vectors")) %dopar% {
                       results <- data.frame()
                       overlap <- paste0(names(gwas)[i], "_SNP_overlap")
                       pval <- paste0(names(gwas)[i], "_SNP_pVal")
                       or <- paste0(names(gwas)[i], "_SNP_or")
                       gwas_selected <- gwas[[i]]
                       results[c(overlap, pval, or), "Gene"] <- AnnotatePublicData(eqtls_gene, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Gene_cases"] <- AnnotatePublicData(eqtls_gene_cases, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Transcript"] <- AnnotatePublicData(eqtls_transcript, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Transcript_cases"] <- AnnotatePublicData(eqtls_transcript_cases, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Exon"] <- AnnotatePublicData(eqtls_exon, gwas_selected, perm, background)
                       results[c(overlap, pval, or), "Exon_cases"] <- AnnotatePublicData(eqtls_exon_cases, gwas_selected, perm, background)
                       write.table(results, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_", names(gwas)[i],"_p0.05.tsv"), row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
                       results
                     }
stopCluster(cl)

# 2.4.6 write results into a file
write.table(results2, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/07_annotations/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_ALL_GWAS_hits_p0.05.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")

print(sessionInfo())

