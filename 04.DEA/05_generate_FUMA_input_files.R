#DEA
exon = read.delim("1.3_diff_expr_exons_for_diagnosis_new_pvalue_and_fdr.csv")
#dim(exon)
#[1]17522     
exon.sig = exon[ exon$FDR.estimate<=0.1,]
dim(exon.sig)#2224

trans= read.delim("1.3_diff_expr_transcripts_for_diagnosis_new_pvalue_and_fdr.csv")
#dim(trans) 17507
trans.sig = trans[ trans$FDR.estimate<=0.1,]
dim(trans.sig)#6

gene = read.delim("1.1_diff_expr_genes_for_diagnosis.csv")
gene.sig = trans[gene$FDR.estimate<=0.1,]
dim(gene.sig)#0

intersect(exon.sig$Gene.Symbol, trans.sig$Gene.Symbol)

#generate bg and Fuma input
write.table(unique(exon.sig$Gene.Symbol), file="1.3_exon_genes.txt", row.names = F,col.names = F,quote = F)
write.table(unique(exon$Gene.Symbol), file="1.3_exon_bg.txt", row.names = F,col.names = F,quote = F)

write.table(trans.sig$Gene.Symbol, file="1.3_trans_genes.txt", row.names = F,col.names = F,quote = F)
write.table(unique(trans$Gene.Symbol), file="1.3_trans_bg.txt", row.names = F,col.names = F,quote = F)

write.table(unique(gene$Gene.Symbol), file="1.3_gene_bg.txt", row.names = F,col.names = F,quote = F)
