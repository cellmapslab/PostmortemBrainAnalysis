# add library path
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
library(rtracklayer)
library(GenomicRanges)

##############
### Height ###
##############
# print("Start with Height summary statistics file")
# # 1. Read in GWAS summary statistics file and add "chr" to each chromosome number (item) in the hg18chr column
# height <- read.csv2(file = "/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Height/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt", 
#                     header = TRUE, 
#                     sep = "\t", 
#                     stringsAsFactors = F)
#height[!grepl("rs", height$MarkerName), "MarkerName"]

# # 2. Load HapMap dataset (Downloaded from: https://homepages.uni-regensburg.de/~wit59712/easystrata/hapmap36.map)
# hapmap <- read.csv2(file = "/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Height/hapmap36.map", 
#                     header = TRUE, 
#                     sep = "\t", 
#                     stringsAsFactors = F)
# hapmap$Zero <- NULL

# # 3. Merge HapMap and Height GWAS to get the chromosome and position for each SNP ID
# height.info <- merge(height, hapmap, by.x = "MarkerName", by.y = "SNPID")
# height.info$Chromosome <- paste("chr", height.info$Chromosome, sep="")
#length(height$MarkerName) - length(height.info$MarkerName)

# # 4. Generate a GRanges object from the merged height dataset
# heightGR <- makeGRangesFromDataFrame(height.info, 
#                                      keep.extra.columns = T, 
#                                      seqnames.field = "Chromosome", 
#                                      ignore.strand = T, 
#                                      start.field = "Position", 
#                                      end.field = "Position")
# head(heightGR)

# # 5. Load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
# hg18ToHg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg18ToHg19.over.chain")

# # 6. Do liftover from hg18 to hg19
# heightGR.ranges_hg19 <- unlist(liftOver(heightGR, hg18ToHg19))
# #head(heightGR.ranges_hg19)

# # 7. Save data to a dataframe and reorder to a similar appearance then the other GWAS summary statistic files
# heightDF <- data.frame(Chr = seqnames(heightGR.ranges_hg19),
#                        MarkerName = heightGR.ranges_hg19$MarkerName,
#                        BP = start(heightGR.ranges_hg19),
#                        Allele1 = heightGR.ranges_hg19$Allele1,
#                        Allele2 = heightGR.ranges_hg19$Allele2,
#                        Freq.Allele1.HapMapCEU = heightGR.ranges_hg19$Freq.Allele1.HapMapCEU,
#                        b = heightGR.ranges_hg19$b,
#                        SE = heightGR.ranges_hg19$SE,
#                        p = heightGR.ranges_hg19$p,
#                        N = heightGR.ranges_hg19$N)
# heightDF$Chr <- substring(heightDF$Chr, 4)  # remove chr before the number
# heightDF <- heightDF[order(heightDF$Chr), ]  # order by chromosome
# #heightDF <- heightDF[, c(2, 1, 3, 4:ncol(heightDF))]  # reorder columns to: chr, SNP ID, pos, a1, a2, ...
# #colnames(heightDF) <- c("CHR", "SNP", "BP", "A1", "A2", "FreqAllele1HapMapCEU", "B", "SE", "P", "N")  # rename columns
# head(heightDF)

# # 8. Save data into a file
# write.table(heightDF, 
#             file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/genotype_data/2.1_sumstats_liftover_and_removed_D_I/Height_hg19.sumstats", 
#             sep = "\t",
#             row.names = F,
#             quote = F)

# # 9. Check, how many SNPs we loose in the process
# # we loose 4643 SNP IDs (if we just merge by SNP ID, we remove over 40,000 SNPs)
# lostSNPs <- length(height$MarkerName) - length(heightDF$MarkerName)
# print("Number of lost SNPs in the genotype file: ")
# print(lostSNPs)


######################
### Cross-Disorder ###
######################
print("Start with Cross-Disorder summary statistics file")
# 1. Read in GWAS summary statistics file and add "chr" to each chromosome number (item) in the hg18chr column
crossd <- read.csv2(file = "/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/Cross_Disorder_Full/pgc.cross.full.2013-03.txt", 
                    header = TRUE, 
                    sep = "\t", 
                    stringsAsFactors = F)

# 2. Generate a GRanges object from the cross-disorder dataset
crossdGR <- makeGRangesFromDataFrame(crossd, 
                                     keep.extra.columns = T, 
                                     seqnames.field = "hg18chr", 
                                     ignore.strand = T, 
                                     start.field = "bp", 
                                     end.field = "bp")
# Change chromosome from numbers (1, 2, etc.) to UCSC format (chr1, chr2, etc.)
seqlevelsStyle(crossdGR) <- "UCSC"
#head(crossdGR)

# 3. Load hg18 to hg19 chain (chain downloaded from UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz)
hg18ToHg19 <- import.chain("/storage/groups/ccm01/datasets/2018_GWAS_summary_statistics/hg18ToHg19.over.chain")
#head(hg18ToHg19)

# 4. Do liftover from hg18 to hg19
crossdGR.ranges_hg19 <- unlist(liftOver(crossdGR, hg18ToHg19))
#head(crossdGR.ranges_hg19)

# 5. Save data to a dataframe and reorder to a similar appearance then the other GWAS summary statistic files
crossdDF <- data.frame(hg19chr = seqnames(crossdGR.ranges_hg19),
                       snpid = crossdGR.ranges_hg19$snpid,
                       a1 = crossdGR.ranges_hg19$a1,
                       a2 = crossdGR.ranges_hg19$a2,
                       bp = start(crossdGR.ranges_hg19),
                       info = crossdGR.ranges_hg19$info,
                       or = crossdGR.ranges_hg19$or,
                       se = crossdGR.ranges_hg19$se,
                       p = crossdGR.ranges_hg19$pval,
                       ngt = crossdGR.ranges_hg19$ngt,
                       CEUaf = crossdGR.ranges_hg19$CEUaf)
crossdDF <- crossdDF[order(crossdDF$hg19chr), ]  # order by chromosome
#head(crossdDF)

# 8. Save data into a file
write.table(crossdDF, 
            file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/preprocessed_sumstats/2.1_sumstats_liftover_and_removed_D_I/Cross_disorder_hg19.sumstats", 
            sep = "\t",
            row.names = F,
            quote = F)

# 9. Check, how many SNPs we loose in the process
# we loose 4643 SNP IDs (if we just merge by SNP ID, we remove over 40,000 SNPs)
lostSNPs <- length(crossd$MarkerName) - length(crossdDF$MarkerName)
print("Number of lost SNPs in the genotype file: ")
print(lostSNPs)

# 10. Remove all saved variables
rm(list = ls())
