.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table) 

args <- commandArgs(trailingOnly = T)
input <- args[1]
output <- args[2]
tools <- args[3]

# 1. Load genotype file (one file for each chromosome)
gen.file <- fread(paste0(input, '.out'), data.table = F)

# recode all SNP IDs to chr:pos (e.g. chr1:1234) <- without a1 and a2 as LDpred checks this in the coord_genotypes.py file and we have the problem that we do not know which one is a1/a2
gen.file[, 2] <- paste0("chr", gen.file[, 1], ":", gen.file[, 3]) # with a1 (which is on position two in the MPIP datasets) => gen.file[, 5]
fwrite(gen.file, paste0(output, '_renamed.gen'), row.names = F, col.names = F, sep = " ", quote = F)

# multi-allelic SNPs - keep SNP with highest MAF
system(paste0(tools, "/PLINK/plink  --gen ", output, "_renamed.gen --sample ", input, ".sample --allow-no-sex --freq --out ", output, "_renamed"))
freq <- fread(paste0(output, "_renamed.frq"), data.table = F)
duplicated.snps <- freq[duplicated(freq$SNP), ]$SNP
keep.gen.file <- c()
for(j in 1:length(duplicated.snps)){
	  tmp.sel <- freq[freq$SNP %in% duplicated.snps[j], ]
	  tmp.sel.maxMAF <- tmp.sel[which.max(tmp.sel$MAF), ]
	  keep.gen.file <- rbind(keep.gen.file, tmp.sel.maxMAF)
}

# write not duplicated + one of the duplicated (the one with highest MAF) to output file
gen.file.nodup <- gen.file[!(duplicated(gen.file$V3)), ]
gen.file.keep <- gen.file[row.names(keep.gen.file), ]
fwrite(rbind(gen.file.nodup, gen.file.keep), paste0(output, '_renamed_noMultiAllelicSNPs.gen'), row.names = F, col.names = F, sep = " ", quote = F)

print(sessionInfo())
