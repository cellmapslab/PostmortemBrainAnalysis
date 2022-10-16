### Code for figure 5 D & E ####
library(data.table)
library(oligo)
library(tidyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(viridis)

# function to load the PRS all score file, change the colnames into right order and return only IID and 0.01 PRS
GetPRS <- function(filename){
  all_scores <- suppressWarnings(fread(filename, stringsAsFactors = F))
  # set colnames right 
  names <- colnames(all_scores)[-1] # remove first colname and save all others
  all_scores$`1` <- NULL # remove last column with NAs
  colnames(all_scores) <- names # set right colnames
  # return only IID and PRS for 0.01
  score <- all_scores[, c("IID", "0.01")]
  return(score)  
}

# load exon expression data
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/01_RData_files/7_affymetrix_exon_microarray_data_after_rma_probeset_summaries_filtered_batch_corrected_SVA_exon_summaries_annotated.RData")
exprs(exon_summaries) <- as.matrix(exprs(exon_summaries))

# load SCZ2-2014 PRS
prs_scz2 <- GetPRS("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/02.PRS/04_PRSice-2/SCZ2-2014/pmb_info06_all_chr_maf0.01_IDs_changed_noDup_final_SCZ2-2014_PRSice-2_PRS.all.score")

# Generate pheno
pData(exon_summaries) <- plyr::join(pData(exon_summaries), as.data.frame(prs_scz2))
colnames(pData(exon_summaries))[length(pData(exon_summaries)[1,])] <- "PRS.SCZ2"
pheno <- pData(exon_summaries)[, c("ShortName", "PRS.SCZ2", "Dx")]

# load SCZ2-2019 eQTS expression
load("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results/3.5.14_diff_expr_exons_for_SCZ2-2014_PRS_new_pvalue_and_fdr.RData")

# load exon eQTLs
exon_eqtl <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/04_result_tables/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt")

#############
### USP24 ###
#############
# Scatterplot - USP24 exon expression vs. SCZ GWAS
usp24_exons <- row.names(fData(exon_summaries)[grep("USP24", fData(exon_summaries)$gene_name), ])

# load SCZ2-2019 eQTS expression
eexp <- exprs(exon_summaries)
exon <- rownames(res)[rownames(res) %in% usp24_exons]
prs <- eexp[rownames(eexp) == exon, ]
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "USP24_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = USP24_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) +
  geom_point(size = 2) + 
  xlab(paste0(exon, " expression")) + 
  ylab("SCZ PRS") + 
  ggtitle("USP24") + 
  #stat_ellipse(type = "t") + # draw ellipses
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = 7, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_USP24_", exon, "exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)

# Boxplot - Now we generate the boxplot with the expression of exon ENSE00001139439.2|ENSE00001811443.1 of USP24 (ID is joined with other exon) with SNP rs1475703
# get exon expression
exon_expr <- t(t(exprs(exon_summaries)[grep("ENSE00001139439.2|ENSE00001811443.1", rownames(exon_summaries)), ]))
colnames(exon_expr) <- "ENSE00001139439.2"

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --snp rs1475703 --recode --tab --out USP24_SNPs"))
snp <- read.delim("USP24_SNPs.ped", head = T)
pheno2 <- pData(exon_summaries)[, c("IID","Dx", "ShortName")]

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/USP24_", colnames(exon_expr), "_expression_vs_rs1475703_snp.rda"))

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "T.T", measure.vars = c("ENSE00001139439.2", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$T.T != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(T.T, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs1475703") + 
  ylab(paste0(colnames(exon_expr), " expression")) + 
  annotate(geom = "text", x = 0.8, y = 8, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00001139439.2|ENSE00001811443.1" & exon_eqtl$SNP == "rs1475703", `FDR estimate`], 43))) +
  ggtitle("USP24")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_USP24_ENSE00001139439_expression_vs_SNP_rs1475703.pdf", width = 6, height = 5)

# get SCZ2-2019 eQTS expression (residuals)
exon <- rownames(res)[rownames(res) %in% usp24_exons]
prs <- res[rownames(res)[rownames(res) %in% usp24_exons], ] # use residuals instead of scores
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "USP24_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = USP24_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) +
  geom_point(size = 2) + 
  xlab(paste0(exon, " expression")) + 
  ylab("SCZ PRS") + 
  ggtitle("USP24") + 
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = -0.5, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_USP24_", exon, "_residuals_exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)


<<<<<<< HEAD
# Boxplot with residuals
# get exon expression
results_path = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results"
load(paste0(results_path, "/3.1_diff_expr_exons_for_diagnosis.RData"))
exon_expr <- t(t(res[grep("ENSE00001139439.2|ENSE00001811443.1", rownames(res)), ]))
colnames(exon_expr) <- "ENSE00001139439.2"

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --snp rs1475703 --recode --tab --out USP24_SNPs"))
snp <- read.delim("USP24_SNPs.ped", head = T)
pheno2 <- pData(exon_summaries)[, c("IID","Dx", "ShortName")]

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/USP24_", colnames(exon_expr), "_expression_vs_rs1475703_snp.rda"))

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "T.T", measure.vars = c("ENSE00001139439.2", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$T.T != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(T.T, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs1475703") + 
  ylab(paste0(colnames(exon_expr), " expression")) + 
  annotate(geom = "text", x = 0.8, y = 3.2, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00001139439.2|ENSE00001811443.1" & exon_eqtl$SNP == "rs1475703", `FDR estimate`], 43))) +
  ggtitle("USP24")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_USP24_ENSE00001139439_residuals_exon_expr_vs_SNP_rs1475703.pdf", width = 6, height = 5)


=======
>>>>>>> ffd06f5 (Minor updates)
############
### GAD1 ###
############
# Scatterplot - GAD1 exon expression vs. SCZ GWAS
gad1_exons <- row.names(fData(exon_summaries)[grep("GAD1", fData(exon_summaries)$gene_name), ])

# load SCZ2-2019 eQTS expression
eexp <- exprs(exon_summaries)
exon <- rownames(res)[rownames(res) %in% gad1_exons]
prs <- eexp[rownames(eexp) == exon, ]
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "GAD1_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = GAD1_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) + 
  geom_point(size = 2) + 
  xlab(paste0(strsplit(exon, "\\|")[[1]][1], " expression")) + 
  ylab("SCZ PRS") + 
  ggtitle("GAD1") + 
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = 8, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_GAD1_", exon, "_exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)

<<<<<<< HEAD

=======
>>>>>>> ffd06f5 (Minor updates)
# Boxplot - Now we generate the boxplot with the expression of exon ENSE00003699191.1|ENSE00001955143.1|ENSE00003696650.1|ENSE00001708604.1 of GAD1 (ID is joined with other exons) with SNP rs1420385
# get exon expression
exon_expr <- t(t(exprs(exon_summaries)[grep("ENSE00003699191.1|ENSE00001955143.1|ENSE00003696650.1|ENSE00001708604.1", rownames(exon_summaries)), ]))
colnames(exon_expr) <- "ENSE00003699191.1"

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --snp rs1420385 --recode --tab --out GAD1_SNPs"))
snp <- read.delim("GAD1_SNPs.ped", head = T)
pheno2 <- pData(exon_summaries)[, c("IID","Dx", "ShortName")]

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/GAD1_ENSE00003699191.1_expression_vs_rs1420385_snp.rda")

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "G.G", measure.vars = c("ENSE00003699191.1", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$G.G != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(G.G, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs1420385") + 
  ylab(paste0(colnames(exon_expr), " expression")) +
  annotate(geom = "text", x = 0.8, y = 9.8, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00003699191.1|ENSE00001955143.1|ENSE00003696650.1|ENSE00001708604.1" & exon_eqtl$SNP == "rs1420385", `FDR estimate`], 15))) +
  ggtitle("GAD1")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_GAD1_ENSE00003699191.1_expression_vs_SNP_rs1420385.pdf", width = 6, height = 5)

<<<<<<< HEAD

# Scatterplot with residuals
=======
>>>>>>> ffd06f5 (Minor updates)
# get SCZ2-2019 eQTS expression (residuals)
prs <- res[rownames(res)[rownames(res) %in% gad1_exons], ] # use residuals instead of scores
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "GAD1_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = GAD1_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) + 
  geom_point(size = 2) + 
  xlab(paste0(strsplit(exon, "\\|")[[1]][1], " expression")) + 
  ylab("SCZ PRS") + 
  ggtitle("GAD1") + 
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = -0.9, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_GAD1_", exon, "_residuals_exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)


<<<<<<< HEAD
# Boxplot with residuals
# get exon expression
results_path = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results"
load(paste0(results_path, "/3.1_diff_expr_exons_for_diagnosis.RData"))
exon_expr <- t(t(res[grep("ENSE00003699191.1|ENSE00001955143.1|ENSE00003696650.1|ENSE00001708604.1", rownames(res)), ]))
colnames(exon_expr) <- "ENSE00003699191.1"

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/GAD1_ENSE00003699191.1_expression_vs_rs1420385_snp.rda")

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "G.G", measure.vars = c("ENSE00003699191.1", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$G.G != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(G.G, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs1420385") + 
  ylab(paste0(colnames(exon_expr), " expression")) +
  annotate(geom = "text", x = 0.8, y = 0.8, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00003699191.1|ENSE00001955143.1|ENSE00003696650.1|ENSE00001708604.1" & exon_eqtl$SNP == "rs1420385", `FDR estimate`], 15))) +
  ggtitle("GAD1")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_GAD1_ENSE00003699191.1_residuals_exon_expr_vs_SNP_rs1420385.pdf", width = 6, height = 5)


=======
>>>>>>> ffd06f5 (Minor updates)
##############
### PCDH15 ###
##############
# Scatterplot - PCDH15 exon expression vs. SCZ GWAS
pheno <- pData(exon_summaries)[, c("ShortName", "PRS.SCZ2", "Dx")]
pcdh15_exons <- row.names(fData(exon_summaries)[grep("PCDH15", fData(exon_summaries)$gene_name), ])

# load SCZ2-2019 eQTS expression
eexp <- exprs(exon_summaries)
exon <- rownames(res)[rownames(res) %in% pcdh15_exons]
prs <- eexp[rownames(eexp) == exon, ]
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "PCDH15_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = PCDH15_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) + 
  geom_point(size = 2) + 
  xlab(paste0(exon, " expression")) + 
  xlab("SCZ PRS") + 
  ggtitle("PCDH15") +
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = 2.3, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_PCDH15_", exon, "_exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)

# Boxplot - Now we generate the boxplot with the expression of exon ENSE00001587529.1 of PCDH15 with SNP rs11523830
# get exon expression
exon_expr <- t(t(exprs(exon_summaries)[grep("ENSE00001587529.1", rownames(exon_summaries)), ]))
colnames(exon_expr) <- "ENSE00001587529.1"

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --snp rs11523830 --recode --tab --out PCDH15_SNPs"))
snp <- read.delim("PCDH15_SNPs.ped", head = T)
pheno2 <- pData(exon_summaries)[, c("IID","Dx", "ShortName")]

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/PCDH15_ENSE00001587529.1_expression_vs_rs11523830_snp.rda")

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "T.T", measure.vars = c("ENSE00001587529.1", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$T.T != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(T.T, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs11523830") + 
  ylab(paste0(colnames(exon_expr), " expression")) + 
  annotate(geom = "text", x = 0.8, y = 5, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00001587529.1" & exon_eqtl$SNP == "rs11523830", `FDR estimate`], 3))) +
  ggtitle("PCDH15")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_PCDH15_ENSE00001587529.1_expression_vs_SNP_rs11523830.pdf", width = 6, height = 5)

# get SCZ2-2019 eQTS expression (residuals)
prs <- res[rownames(res)[rownames(res) %in% pcdh15_exons], ] # use residuals instead of scores
prs <- t(t(as.data.frame(prs)))
colnames(prs) <- "PCDH15_SCZ"
plot_data <- merge(pheno, prs, by.x = "ShortName", by.y = "row.names")

# generate scatterplot
ggplot(plot_data, aes(x = PCDH15_SCZ, y = PRS.SCZ2, col = Dx)) + 
  theme_classic(base_size = 20) + 
  geom_point(size = 2) + 
  ylab(paste0(exon, " expression")) + 
  xlab("SCZ PRS") + 
  ggtitle("PCDH15") +
  scale_color_manual(values = brewer.pal(11, "RdBu")[c(9, 3)], name = "Diagnosis", labels = c("Controls", "Cases")) +
  geom_smooth(method = "lm", se = F) + # draw lines
  annotate(geom = "text", x = 0, y = -0.0022, label = paste0("FDR = ", round(top_hits[top_hits$`Ensembl Exon ID` == exon, "FDR estimate"], 3)))
ggsave(file = paste0("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/3_scatterplot_PCDH15_", exon, "_residuals_exon_expr_vs_SCZ2_PRS.pdf"), width = 10, height = 6)

<<<<<<< HEAD

# Boxplot with residuals
# get exon expression
results_path = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results"
load(paste0(results_path, "/3.1_diff_expr_exons_for_diagnosis.RData"))
exon_expr <- t(t(res[grep("ENSE00001587529.1", rownames(res)), ]))
colnames(exon_expr) <- "ENSE00001587529.1"

# snp data
system(paste0("/storage/groups/ccm01/tools/PLINK/plink --bfile /storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/01_preprocessed_genotypes/1.3_samples_filtered/pmb_info06_all_chr_maf0.05_merged --snp rs11523830 --recode --tab --out PCDH15_SNPs"))
snp <- read.delim("PCDH15_SNPs.ped", head = T)
pheno2 <- pData(exon_summaries)[, c("IID","Dx", "ShortName")]

exp_p <- merge(exon_expr, pheno2, by.x = "row.names", by.y = "ShortName")
exp_p_snp <- merge(exp_p, snp, by.x = "IID", by.y = "X2.1")
save(exp_p_snp, file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/PCDH15_ENSE00001587529.1_expression_vs_rs11523830_snp.rda")

# Generate plot
plot_data <- reshape2::melt(exp_p_snp, snp = "T.T", measure.vars = c("ENSE00001587529.1", grep("ENS", colnames(exp_p_snp), value = T)))
plot_data <- plot_data[plot_data$T.T != "0 0",]
plot_data$variable <- as.character(plot_data$variable)
# plot figure
ggplot(na.omit(plot_data), aes(T.T, value)) +
  geom_boxplot() + 
  theme_classic(base_size = 20) +
  theme(strip.text = element_text(size = 12), axis.text = element_text(size = 9), 
        axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16), 
        plot.title = element_text(hjust = 0.5, size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("rs11523830") + 
  ylab(paste0(colnames(exon_expr), " expression")) + 
  annotate(geom = "text", x = 0.8, y = 1.8, label = paste0("FDR = ", round(exon_eqtl[exon_eqtl$`Ensembl Exon ID` == "ENSE00001587529.1" & exon_eqtl$SNP == "rs11523830", `FDR estimate`], 3))) +
  ggtitle("PCDH15")
ggsave(file = "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/04_figures/4_boxplot_PCDH15_ENSE00001587529.1_residuals_exon_expr_vs_SNP_rs11523830.pdf", width = 6, height = 5)

=======
>>>>>>> ffd06f5 (Minor updates)
