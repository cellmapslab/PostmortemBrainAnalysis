# Load locally installed packages for R version 3.6
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)

### Functions ###
# function to get all unique SNPs of the eQTL results
GetUniqueSnps <- function(eqtls, file){
  unique_snps <- as.data.frame(unique(eqtls$SNP))
  colnames(unique_snps) <- "SNP"
  fwrite(unique_snps, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", file), col.names = F, sep = "\n")
  return(unique_snps)
}

# function to get the chromosome and minimal p-value for all unique SNPs of eQTL results
GetChrAndPvalueOfUniqueSnps <- function(eqtls, snps, file){
  min_pval_snps <- eqtls[ , .SD[which.min(`p-value`)], by = `SNP`]  # get minimal p-value for each SNP
  unique_snps <- merge(snps, min_pval_snps[, c("SNP", "SNP_location", "p-value")], by = "SNP")
  unique_snps$CHR <- sapply(unique_snps$SNP_location, function(x) strsplit(x, ":")[[1]][1])
  unique_snps$P <- unique_snps$`p-value`
  unique_snps <- unique_snps[, c("SNP", "CHR", "P")]
  print(head(unique_snps))
  print(dim(unique_snps))
  fwrite(unique_snps, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", file), sep = "\t")
}

# function to load the clumping results, reduce eQTL results to top LD SNPs and save the resulting eQTLs into an extra file
GetClumpedSnpsEqtls <- function(input, eqtls, output){
  clumped_snps <- fread(paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/", input), header = T, select = "SNP", stringsAsFactors = F)
  eqtls <- eqtls[SNP %in% clumped_snps$SNP, ]
  print(dim(eqtls))
  write.table(eqtls, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/", output), row.names = F, quote = F, dec = ".", sep = " ")
}

### 1. LOAD eQTL results ###
# 1.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
eqtls_gene <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls_transcript <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
eqtls_exon <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_vs_controls.txt", header = T, stringsAsFactors = F)
# 1.2 cis-eQTLs of cases only in genes, transcripts and exons
eqtls_gene_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls_transcript_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_transcript_hits_cases_only.txt", header = T, stringsAsFactors = F)
eqtls_exon_cases <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_exon_hits_cases_only.txt", header = T, stringsAsFactors = F)

# 2. Generate list of all SNP rs IDs, which we want to do clumping for
# 2.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
gene_snps <- GetUniqueSnps(eqtls_gene, file = "01_pmb_cis-eQTL_gene_SNP_IDs.txt")
transcript_snps <- GetUniqueSnps(eqtls_transcript, file = "01_pmb_cis-eQTL_transcript_SNP_IDs.txt")
exon_snps <- GetUniqueSnps(eqtls_exon, file = "01_pmb_cis-eQTL_exon_SNP_IDs.txt")
# 2.2 cis-eQTLs of cases only in genes, transcripts and exons
gene_cases_snps <- GetUniqueSnps(eqtls_gene_cases, file = "01_pmb_cis-eQTL_gene_cases_SNP_IDs.txt")
transcript_cases_snps <- GetUniqueSnps(eqtls_transcript_cases, file = "01_pmb_cis-eQTL_transcript_cases_SNP_IDs.txt")
exon_cases_snps <- GetUniqueSnps(eqtls_exon_cases, file = "01_pmb_cis-eQTL_exon_cases_SNP_IDs.txt")

# 3. Generate list of all SNP rs IDs, which we want to do clumping for including their corrresponding chromosome and minimal eQTL p-value
# 3.1 cis-eQTLs of cases and controls (ALL) in genes, transcripts and exons
GetChrAndPvalueOfUniqueSnps(eqtls_gene, gene_snps, file = "02_pmb_cis-eQTL_gene_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls_transcript, transcript_snps, file = "02_pmb_cis-eQTL_transcript_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls_exon, exon_snps, file = "02_pmb_cis-eQTL_exon_clumping.txt")
# 3.2 cis-eQTLs of cases only in genes, transcripts and exons
GetChrAndPvalueOfUniqueSnps(eqtls_gene_cases, gene_cases_snps, file = "02_pmb_cis-eQTL_gene_cases_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls_transcript_cases, transcript_cases_snps, file = "02_pmb_cis-eQTL_transcript_cases_clumping.txt")
GetChrAndPvalueOfUniqueSnps(eqtls_exon_cases, exon_cases_snps, file = "02_pmb_cis-eQTL_exon_cases_clumping.txt")

rm(gene_snps, transcript_snps, exon_snps, gene_cases_snps, transcript_cases_snps, exon_cases_snps)
gc()

# 3. Run clumping using PLINK
# 3.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_gene_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_gene_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_gene")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_transcript_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_transcript_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_transcript")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_exon_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_exon_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_exon")
# 3.2 cis-eQTLs of cases only in genes, transcripts and exons
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_gene_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_gene_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_gene_cases")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_transcript_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_transcript_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_transcript_cases")
system("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/preprocessed_genotypes/1.3_duplicated_SNPs_removed/pmb_info06_all_chr_maf0.05_merged_noDup --extract /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/01_pmb_cis-eQTL_exon_cases_SNP_IDs.txt --clump /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/02_pmb_cis-eQTL_exon_cases_clumping.txt --clump-r2 0.2 --clump-p1 0.05 --clump-p2 1 --clump-kb 1000 --out /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/clumping/03_pmb_cis-eQTL_exon_cases")

# 4. Get clumped SNP eQTLs
# 4.1 cis-eQTLs of cases and controls (ALL) of genes, transcripts and exons
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_gene.clumped", eqtls_gene, output = "pmb_cis-eQTL_gene_results_cases_vs_controls_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_transcript.clumped", eqtls_transcript, output = "pmb_cis-eQTL_transcript_results_cases_vs_controls_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_exon.clumped", eqtls_exon, output = "pmb_cis-eQTL_exon_results_cases_vs_controls_clumped.txt")
# 4.2 cis-eQTLs of cases only in genes, transcripts and exons
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_gene_cases.clumped", eqtls_gene_cases, output = "pmb_cis-eQTL_gene_results_cases_only_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_transcript_cases.clumped", eqtls_transcript_cases, output = "pmb_cis-eQTL_transcript_results_cases_only_clumped.txt")
GetClumpedSnpsEqtls(input = "03_pmb_cis-eQTL_exon_cases.clumped", eqtls_exon_cases, output = "pmb_cis-eQTL_exon_results_cases_only_clumped.txt")

#eqtls_gene_clumped <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/results/pmb_cis-eQTL_gene_results_cases_only_clumped.txt")
#length(unique(sort(eqtls_gene_clumped$SNP)))
#length(unique(sort(eqtls_gene_clumped$Gene)))
#length(unique(sort(eqtls_gene_clumped$Transcript)))
#length(eqtls_gene_clumped$SNP)

