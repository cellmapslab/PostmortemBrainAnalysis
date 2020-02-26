.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.6/")

library(data.table)

input_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input/"
output_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Smoothing/input/expression/"

# 1. Generate gene expression input files
gene_expr <- fread(paste0(input_path, "01_pmb_gene_probeset_expr_CASE_CTL.txt"), header = T, stringsAsFactors = F, data.table = F)
count <- 1
for (i in seq(1, nrow(gene_expr), 10000)) {
  if (i < 10000) {
    chunk <- gene_expr[i : (i + 10000 - 1), ]
    # Save gene expression data into a new chunk file
    oname <- paste0(output_path, "gene/01_pmb_gene_probeset_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- gene_expr[i : nrow(gene_expr), ]
    # Save gene expression data into a new chunk file
    oname <- paste0(output_path, "gene/01_pmb_gene_probeset_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}

# 2. Generate transcript expression input files
transcript_expr <- fread(paste0(input_path, "01_pmb_transcript_expr_CASE_CTL.txt"), header = T, stringsAsFactors = F, data.table = F)
count <- 1
for (i in seq(1, nrow(transcript_expr), 10000)) {
  if (i < 100000) {
    chunk <- transcript_expr[i : (i + 10000 - 1), ]
    # Save transcript expression data into a new chunk file
    oname <- paste0(output_path, "transcript/01_pmb_transcript_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- transcript_expr[i : nrow(transcript_expr), ]
    # Save transcript expression data into a new chunk file
    oname <- paste0(output_path, "transcript/01_pmb_transcript_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}

# 3. Generate exon expression input files
exon_expr <- fread(paste0(input_path, "01_pmb_exon_expr_CASE_CTL.txt"), header = T, stringsAsFactors = F, data.table = F)
count <- 1
for (i in seq(1, nrow(exon_expr), 50000)) {
  if (i < 200000) {
    chunk <- exon_expr[i : (i + 50000 - 1), ]
    # Save exon expression data into a new chunk file
    oname <- paste0(output_path, "exon/01_pmb_exon_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  } else {
    chunk <- exon_expr[i : nrow(exon_expr), ]
    # Save exon expression data into a new chunk file
    oname <- paste0(output_path, "exon/01_pmb_exon_expr_CASE_CTL_chunk", count, ".txt")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count <- count + 1
  }
}

