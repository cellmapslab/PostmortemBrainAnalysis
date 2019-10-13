# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)

### Functions ###
# function to get all unique SNPs of the eQTL hits
GetUniqueSnps <- function(eqtls, file){
  unique.snps <- as.data.frame(unique(eqtls$SNP))
  colnames(unique.snps) <- "SNP"
  fwrite(unique.snps, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", file), col.names = F, sep = "\n")
  return(unique.snps)
}

# function to get the chromosome and minimal p-value for all unique SNPs of eQTL hits
GetChrAndPvalueOfUniqueSnps <- function(eqtls, snps, file){
  min.pval.snps <- eqtls[ , .SD[which.min(`p-value`)], by = `SNP`]  # get minimal p-value for each SNP
  unique.snps <- merge(snps, min.pval.snps[, c("SNP", "SNP_location", "p-value")], by = "SNP")
  unique.snps$CHR <- sapply(unique.snps$SNP_location, function(x) strsplit(x, ":")[[1]][1])
  unique.snps$P <- unique.snps$`p-value`
  unique.snps <- unique.snps[, c("SNP", "CHR", "P")]
  print(head(unique.snps))
  print(dim(unique.snps))
  fwrite(unique.snps, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", file), sep = "\t")
}

# function to load the clumping results, reduce eQTL hits to top LD SNPs and save the resulting eQTLs into an extra file
GetClumpedSnpsEqtls <- function(input, eqtls, output){
  clumped.snps <- fread(paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", input), header = T, select = "SNP", stringsAsFactors = F)
  eqtls <- eqtls[SNP %in% clumped.snps$SNP, ]
  print(dim(eqtls))
  write.table(eqtls, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/", output), row.names = F, quote = F, dec = ".", sep = " ")
}

### 1. LOAD eQTL results ###
# 1.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls.gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls.exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# 1.2 cis-eQTLs of cases only in genes, transcripts and exons
eqtls.gene.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.transcript.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls.exon.cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/hits/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

# 2. Generate list of all SNP rs IDs, which we want to do clumping for
# 2.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
gene.snps <- GetUniqueSnps(eqtls.gene, file = "01_pmb_cis-eQTL_gene_SNP_IDs.txt")
transcript.snps <- GetUniqueSnps(eqtls.transcript, file = "01_pmb_cis-eQTL_transcript_SNP_IDs.txt")
exon.snps <- GetUniqueSnps(eqtls.exon, file = "01_pmb_cis-eQTL_exon_SNP_IDs.txt")
# 2.2 cis-eQTLs of cases only in genes, transcripts and exons
gene.cases.snps <- GetUniqueSnps(eqtls.gene.cases, file = "01_pmb_cis-eQTL_gene_cases_SNP_IDs.txt")
transcript.cases.snps <- GetUniqueSnps(eqtls.transcript.cases, file = "01_pmb_cis-eQTL_transcript_cases_SNP_IDs.txt")
exon.cases.snps <- GetUniqueSnps(eqtls.exon.cases, file = "01_pmb_cis-eQTL_exon_cases_SNP_IDs.txt")

# 3. Generate list of all SNP rs IDs, which we want to do clumping for including their corrresponding chromosome and minimal eQTL p-value
# 3.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
GetChrAndPvalueOfUniqueSnps(eqtls.gene, gene.snps, file = "02_pmb_cis-eQTL_gene_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls.transcript, transcript.snps, file = "02_pmb_cis-eQTL_transcript_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls.exon, exon.snps, file = "02_pmb_cis-eQTL_exon_clumping.txt")
# 3.2 cis-eQTLs of cases only in genes, transcripts and exons
GetChrAndPvalueOfUniqueSnps(eqtls.gene.cases, gene.cases.snps, file = "02_pmb_cis-eQTL_gene_cases_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls.transcript.cases, transcript.cases.snps, file = "02_pmb_cis-eQTL_transcript_cases_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls.exon.cases, exon.cases.snps, file = "02_pmb_cis-eQTL_exon_cases_clumping.txt")

rm(gene.snps, transcript.snps, exon.snps, gene.cases.snps, transcript.cases.snps, exon.cases.snps)
gc()

# 3. Run clumping using PLINK
# 3.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_gene_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_gene_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_gene")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_transcript_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_transcript_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_transcript")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_exon_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_exon_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_exon")
# 3.2 cis-eQTLs of cases only in genes, transcripts and exons
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_gene_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_gene_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_gene_cases")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_transcript_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_transcript_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_transcript_cases")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.2_merged_chr/pmb_info06_all_chr_maf0.05_merged_noDup_final --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_exon_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_exon_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_exon_cases")

# 4. Get clumped SNP eQTLs
# 4.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_gene.clumped", eqtls.gene, output = "pmb_cis-eQTL_gene_hits_cases_vs_controls_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_transcript.clumped", eqtls.transcript, output = "pmb_cis-eQTL_transcript_hits_cases_vs_controls_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_exon.clumped", eqtls.exon, output = "pmb_cis-eQTL_exon_hits_cases_vs_controls_clumped.txt")
# 4.2 cis-eQTLs of cases only in genes, transcripts and exons
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_gene_cases.clumped", eqtls.gene.cases, output = "pmb_cis-eQTL_gene_hits_cases_only_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_transcript_cases.clumped", eqtls.transcript.cases, output = "pmb_cis-eQTL_transcript_hits_cases_only_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_exon_cases.clumped", eqtls.exon.cases, output = "pmb_cis-eQTL_exon_hits_cases_only_clumped.txt")

