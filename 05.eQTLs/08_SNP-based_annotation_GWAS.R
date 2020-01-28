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

# function, to calculate the annotation for one public dataset to all eQTL datasets
# @param input  Public dataset (e.g. CMC eQTLs, GWAS AD-2019)
# @param name   Name of the public dataset
# @param res    Data frame in which the results should be stored
# @return       the results data frame including the newly calculated results
RunAnnotationForAllEQTLs <- function(input, name, results){
  res <- data.frame()
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene"] <- AnnotatePublicData(gene_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Gene_cases"] <- AnnotatePublicData(gene_cases_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript"] <- AnnotatePublicData(transcript_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Transcript_cases"] <- AnnotatePublicData(transcript_cases_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon"] <- AnnotatePublicData(exon_eqtl, input, perm, background)
  res[c(paste0(name, "_overlap"), paste0(name, "_pVal"), paste0(name, "_odds_ratio")), "Exon_cases"] <- AnnotatePublicData(exon_cases_eqtl, input, perm, background)
  # Write single runs into a file
  write.table(res, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/GWAS_CMC_GTEx/pmb_cis-eQTL_annotation_results_10000_permutations_", name, ".tsv"), row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
  results <- rbind(results, res)
  print(paste0("Done ", name))
  return(results)
}

# path, where all datasets are stored (as .rds files)
ipath <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/RData/"

### Load eQTL results ###
gene_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # gene based eQTLs
gene_cases_eqtl <- readRDS(paste0(ipath, "pmb_gene_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only gene based eQTLs
transcript_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # transcript based eQTLs (similar to splice QTLs)
transcript_cases_eqtl <- readRDS(paste0(ipath, "pmb_transcript_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only transcript eQTLs (similar to splice QTLs)
exon_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_vs_controls_cis_eQTL_SNPs.rds"))  # exon based eQTLs (eeQTLs)
exon_cases_eqtl <- readRDS(paste0(ipath, "pmb_exon_hits_cases_only_cis_eQTL_SNPs.rds"))  # cases only exon based eQTLs (eeQTLs)

### Load Background data ###
background <- readRDS(paste0(ipath, "background_SNPs.rds"))


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

# save names of all GWAS in the order they should be running
#gwas <- GRangesList(scz2, bip, cdg1, cdg2, mdd1, mdd2, mdd3, adhd, asd, ea, bmi, ad, alcdep, an, audit, cad, cannabis, ed, height, intelligence, ocd, ptsd, t2d, ts)
#gwas <- c("scz2", "bip", "cdg1", "cdg2", "mdd1", "mdd2", "mdd3", "adhd", "asd", "ea", "bmi", "ad", "alcdep", "an", "audit", "cad", "cannabis", "ed", "height", "intelligence", "ocd", "ptsd", "t2d", "ts")
gwas_gw <- GRangesList(scz2_gw, bip_gw, cdg1_gw, cdg2_gw, mdd1_gw, mdd2_gw, mdd3_gw, adhd_gw, asd_gw, ea_gw, bmi_gw, ad_gw, alcdep_gw, an_gw, audit_gw, cad_gw, cannabis_gw, ed_gw, height_gw, intelligence_gw, ocd_gw, ptsd_gw, t2d_gw, ts_gw)
#gwas_gw <- c("scz2_gw", "bip_gw", "cdg1_gw", "cdg2_gw", "mdd1_gw", "mdd2_gw", "mdd3_gw", "adhd_gw", "asd_gw", "ea_gw", "bmi_gw", "ad_gw", "alcdep_gw", "an_gw", "audit_gw", "cad_gw", "cannabis_gw", "ed_gw", "height_gw", "intelligence_gw", "ocd_gw", "ptsd_gw", "t2d_gw", "ts_gw")

### Annotation to GWAS summary statistics ###
# Set parameters
results <- setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Gene", "Gene_cases", "Transcript", "Transcript_cases", "Exon", "Exon_cases"))
perm <- 10000

# Run annotation for all GWAS SNPs
results <- RunAnnotationForAllEQTLs(ad, "GWAS_AD-2019", results)
results <- RunAnnotationForAllEQTLs(adhd, "GWAS_ADHD-2017", results)
results <- RunAnnotationForAllEQTLs(alcdep, "GWAS_ALCDEP-2018", results)
results <- RunAnnotationForAllEQTLs(an, "GWAS_AN-2019", results)
results <- RunAnnotationForAllEQTLs(asd, "GWAS_ASD-2017", results)
results <- RunAnnotationForAllEQTLs(audit, "GWAS_AUDIT-2018", results)
results <- RunAnnotationForAllEQTLs(bip, "GWAS_BIP-2018", results)
results <- RunAnnotationForAllEQTLs(bmi, "GWAS_BMI-2015", results)
results <- RunAnnotationForAllEQTLs(cad, "GWAS_CAD_2017", results)
results <- RunAnnotationForAllEQTLs(cannabis, "GWAS_Cannabis-2018", results)
results <- RunAnnotationForAllEQTLs(cdg1, "GWAS_CDG1-2013", results)
results <- RunAnnotationForAllEQTLs(cdg2, "GWAS_CDG2-2019", results)
results <- RunAnnotationForAllEQTLs(ea, "GWAS_EA-2018", results)
results <- RunAnnotationForAllEQTLs(ed, "GWAS_ED-2017", results)
results <- RunAnnotationForAllEQTLs(height, "GWAS_Height-2014", results)
results <- RunAnnotationForAllEQTLs(intelligence, "GWAS_Intelligence-2018", results)
results <- RunAnnotationForAllEQTLs(mdd1, "GWAS_MDD1-2012", results)
results <- RunAnnotationForAllEQTLs(mdd2, "GWAS_MDD2-2018", results)
results <- RunAnnotationForAllEQTLs(mdd3, "GWAS_MDD3-2019", results)
results <- RunAnnotationForAllEQTLs(ocd, "GWAS_OCD-2017", results)
results <- RunAnnotationForAllEQTLs(ptsd, "GWAS_PTSD-2019", results)
results <- RunAnnotationForAllEQTLs(scz2, "GWAS_SCZ-2014", results)
results <- RunAnnotationForAllEQTLs(t2d, "GWAS_T2D-2017", results)
results <- RunAnnotationForAllEQTLs(ts, "GWAS_TS-2018", results)

# Write results into a file
write.table(results, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_annotation_results_10000_permutations_GWAS_summary_statistics.tsv", row.names = T, col.names = T, quote = F, dec = ".", sep = "\t")
