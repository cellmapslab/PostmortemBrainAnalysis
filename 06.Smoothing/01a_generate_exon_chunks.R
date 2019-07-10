.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")

library(data.table)

# 1. Load exon expression input file
exon_expr <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/input_files/01_pmb_exon_expr_SCZ_CTL.txt", header = T, stringsAsFactors = F, data.table = F)
count = 1
for (i in seq(1, nrow(exon_expr), 10000)) {
  if (i < 530000) {
    chunk = exon_expr[i : (i + 10000 - 1), ]
    # Save gene expression data into a new chunk file
    oname = paste("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/input/expression/01_pmb_exon_expr_SCZ_CTL_chunk", count, ".txt", sep = "")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count = count + 1
  } else {
    chunk = exon_expr[i : nrow(exon_expr), ]
    # Save gene expression data into a new chunk file
    oname = paste("/storage/scratch/users/karolina.worf/MPIP/06.Smoothing/input/expression/01_pmb_exon_expr_SCZ_CTL_chunk", count, ".txt", sep = "")
    fwrite(x = chunk, file = oname, row.names = F, col.names = T, sep = "\t", quote = F)
    # increase count
    count = count + 1
  }
}
