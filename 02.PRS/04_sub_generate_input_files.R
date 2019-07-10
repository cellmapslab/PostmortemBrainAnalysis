# add path for locally installed R-packages for R version 3.5
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

library(data.table) 

# get command line parameters (path + filename for summary statistic file and genotype file for input and output)
args <- commandArgs(trailingOnly = T)
sumstats.in <- args[1]
sumstats.out <- args[2]
geno.in <- args[3]
cohort.in <- args[4]
snps <- args[5]

# read in .bim file of the samples SNPs file
gbim <- fread(geno.in, data.table = F)

# read in .bim file of the cohort SNPs file
cbim <- fread(cohort.in, data.table = F)

# read in summary statistics file
sum.stats <- fread(sumstats.in, data.table = F)

# # get overlapping SNPs (using the ID, as we transformed them into chr:pos)
overlapping.snps <- intersect(intersect(cbim$V2, gbim$V2), sum.stats[, 2]) 
# print(head(sum.stats.snps))

# save only the overlapping SNPs into LDpred input summary statistics file
fwrite(sum.stats[sum.stats[, 2] %in% overlapping.snps, ], sumstats.out, row.names = F, col.names = T, sep = "\t", quote = F)

# save overlapping SNPs into an extra file, which can be used in plink to get only the overlapping SNPs for the genotype data as well
write.table(overlapping.snps, snps, row.names = F, col.names = F, sep = "\t", quote = F)

