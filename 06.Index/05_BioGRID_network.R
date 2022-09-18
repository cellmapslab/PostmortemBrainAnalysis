### BioGrid network for core genes ###
library(data.table)
library(aliases2entrez) # to map gene symbols to entrez ids
library(igraph)
library(gplots)
library(dplyr)


# load Biogrid data
biogrid <- fread("/storage/groups/ccm01/datasets/2021_BIOGRID-4.4.201/BIOGRID-ORGANISM-Homo_sapiens-4.4.201.tab3.txt")
dim(biogrid)
# remove self interactions
biogrid <- biogrid[biogrid$`Entrez Gene Interactor A` != biogrid$`Entrez Gene Interactor B`, ]
dim(biogrid)

# load 29 core genes
core29 <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/29_core_genes_diff_expr_in_cell_types.tsv", header = FALSE)

# load 832 core genes
core <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/02_overlap/4_overlap_eeGenes_ABMS_rare.txt")

#convert gene symbols to Entrez IDs
# import the correspondence file
file <- system.file("extdata", "HGNC.txt", package = "aliases2entrez")
HGNC <- read.delim(file)
# run the main function
core29 <- convert_symbols(core29$V1, HGNC)
core <- convert_symbols(core$`Gene Symbol`, HGNC)

# now we can get the Biogrid interactions only for these 29 genes
subnetwork_only_core <- biogrid[biogrid$`Entrez Gene Interactor A` %in% core29$entrezID & biogrid$`Entrez Gene Interactor B` %in% core29$entrezID, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database")]
subnetwork_core_to_all <- biogrid[biogrid$`Entrez Gene Interactor A` %in% core29$entrezID | biogrid$`Entrez Gene Interactor B` %in% core29$entrezID, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database")]

# Get all core genes to other gene connections
cores_left <- subnetwork_core_to_all[subnetwork_core_to_all$`Entrez Gene Interactor A` %in% core29$entrezID, `Entrez Gene Interactor A`:`Entrez Gene Interactor B`]
cores_right <- subnetwork_core_to_all[subnetwork_core_to_all$`Entrez Gene Interactor B` %in% core29$entrezID, `Entrez Gene Interactor A`:`Entrez Gene Interactor B`]
setcolorder(cores_right, c("Entrez Gene Interactor B", "Entrez Gene Interactor A")) # set the order of the columns
colnames(cores_right) <- colnames(cores_left) # change name of columns
subnetwork_core_to_all_ids <- merge(cores_left, cores_right, all = TRUE) # merge the two data tables

# check if in column "Entrez Gene Interactor A" includes only core genes
test <- unique(subnetwork_core_to_all_ids$`Entrez Gene Interactor A`)
length(test)
test[!(test %in% as.character(core29$entrezID))]
rm(test)

# check if a connection was found at least twice
subnetwork_core_to_all_ids$`Entrez Gene Interactor C` <- paste(subnetwork_core_to_all_ids$`Entrez Gene Interactor A`, subnetwork_core_to_all_ids$`Entrez Gene Interactor B`, sep = "-")
frequencies <- as.data.frame(table(subnetwork_core_to_all_ids$`Entrez Gene Interactor C`))
subnetwork_core_to_all_ids <- subnetwork_core_to_all_ids[subnetwork_core_to_all_ids$`Entrez Gene Interactor C` %in% frequencies[frequencies$Freq > 1, 'Var1'], ]
biogrid$`Entrez Gene Interactor C` <- paste(biogrid$`Entrez Gene Interactor A`, biogrid$`Entrez Gene Interactor B`, sep = "-")
# get connections to other genes with at least two experiments
subnetwork_core_with_1_inbetween <- biogrid[biogrid$`Entrez Gene Interactor C` %in% subnetwork_core_to_all_ids$`Entrez Gene Interactor C`, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database")]
# add connections between two core genes
subnetwork_core_with_1_inbetween <- merge(subnetwork_core_with_1_inbetween, subnetwork_only_core, all = T)
# add connections between two non core genes (only if freq > 1)
non_core_genes <- unique(subnetwork_core_with_1_inbetween[!(subnetwork_core_with_1_inbetween$`Entrez Gene Interactor B` %in% core29$entrezID), ]$`Entrez Gene Interactor B`)
subnetwork_non_core_genes <- biogrid[biogrid$`Entrez Gene Interactor A` %in% non_core_genes & biogrid$`Entrez Gene Interactor B` %in% non_core_genes, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database", "Entrez Gene Interactor C")]
frequencies <- as.data.frame(table(subnetwork_non_core_genes$`Entrez Gene Interactor C`))
subnetwork_non_core_genes <- subnetwork_non_core_genes[subnetwork_non_core_genes$`Entrez Gene Interactor C` %in% frequencies[frequencies$Freq > 1, 'Var1'], ]
subnetwork_core_with_1_inbetween <- merge(subnetwork_core_with_1_inbetween, subnetwork_non_core_genes, all = T)

# get the edges
subnetwork_core_with_1_inbetween_edges <- list()
subnetwork_core_with_1_inbetween_edges <- unname(unlist(apply(subnetwork_core_with_1_inbetween[, c("Official Symbol Interactor A", "Official Symbol Interactor B")], 1, function(x) append(subnetwork_core_with_1_inbetween_edges, x))))

# plot subnetwork of biogrid
g <- make_empty_graph() + vertices(unique(c(subnetwork_core_with_1_inbetween_edges, as.character(core29$Symbols))))
g <- make_empty_graph() + vertices(unique(subnetwork_core_with_1_inbetween_edges))

g <- g + edges(subnetwork_core_with_1_inbetween_edges)
g <- as.undirected(g, mode = "collapse", edge.attr.comb = igraph_opt("edge.attr.comb"))
g$layout <- layout_with_fr
l <- layout.fruchterman.reingold(g) 
node.size= c(10,10,10)
plot(g)


genes <- as.data.frame(unique(subnetwork_core_with_1_inbetween_edges))
colnames(genes) <- "Symbol"
genes$Group <- "A"
genes$Symbol <- as.character(genes$Symbol)
# set core genes to B and check their synonyms (only one gene found!)
genes[genes$Symbol %in% core$Symbols, "Group"] <- "B"
for (i in 1:length(genes$Symbol)) {
  if(genes[i, 2] == "A"){
    synonyms <- strsplit(unique(biogrid[biogrid$`Official Symbol Interactor A` == genes[i, 1], "Synonyms Interactor A"])$`Synonyms Interactor A`, "\\|")
    core_gene <- core[core$Symbols %in% unlist(synonyms), ]
    if(length(core_gene$Symbols) > 0){
      print(genes[i, 1])
      print(core[core$Symbols %in% unlist(synonyms), ])
    }
  }
}
genes[genes$Symbol == "DYNLL1", "Group"] <- "B"
genes[genes$Symbol %in% core29$Symbols, "Group"] <- "C"

write.table(genes, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/core_genes_subnetwork_biogrid_groups.txt", quote = F, sep = "\t", row.names = F)
write.table(subnetwork_core_with_1_inbetween, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/core_genes_subnetwork_biogrid.txt", quote = F, row.names = F, sep = "\t")

# overlap to 832 core genes
genes$coreALL <- "0"
genes[genes$Symbol %in% core$Symbols, "coreALL"] <- "1"
genes[genes$Symbol == "DYNLL1", "coreALL"] <- "1"

# overlap to 29 core genes
genes$core29 <- "0"
genes[genes$Symbol %in% core29$Symbols, "core29"] <- "1"

# overlap to exon DEGenes
dea_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/04.DEA/03_DEA_results"
eDEGenes <- fread(paste0(dea_path, "/3.4_diff_expr_exons_for_diagnosis_new_pvalue_and_fdr.csv"), data.table = F, stringsAsFactors = F)
eDEGenes <- eDEGenes[eDEGenes$`FDR estimate` < 0.1, ]
eDEGenes$`Gene Symbol` <- sapply(unlist(strsplit(eDEGenes$`Gene Symbol`, "\\|")), function(x) strsplit(x, "\\.")[[1]][1])
genes$exonDEGenes <- "0"
genes[genes$Symbol %in% eDEGenes$`Gene Symbol`, "exonDEGenes"] <- "1"

# overlap to exon eGenes
eqtl_path <- "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/05.eQTLs/04_result_tables"
eqtls_exon <- fread(paste0(eqtl_path, "/01_pmb_cis-eQTL_exon_hits_cases_vs_controls.txt"), header = T, stringsAsFactors = F)
genes$exoneGenes <- "0"
genes[genes$Symbol %in% eqtls_exon$`Gene Symbol`, "exoneGenes"] <- "1"

# overlap to exon eQT-Scores
# Load ADHD-2019 eQTS
adhd <- fread(paste0(dea_path, "/3.5.1_diff_expr_exons_for_ADHD-2019_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
adhd <- adhd[adhd$`FDR estimate` <= 0.1, ]
genes$ADHDeQTS <- "0"
genes[genes$Symbol %in% adhd$`Gene Symbol`, "ADHDeQTS"] <- "1"


# Load ASD-2019 eQTS
asd <- fread(paste0(dea_path, "/3.5.2_diff_expr_exons_for_ASD-2019_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
asd <- asd[asd$`FDR estimate` <= 0.1, ]
genes[genes$Symbol %in% asd$`Gene Symbol`, ]$Symbol
dim(genes[genes$Symbol %in% asd$`Gene Symbol`, ])

# Load BIP-2019 eQTS
bip <- fread(paste0(dea_path, "/3.5.3_diff_expr_exons_for_BIP-2019_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
bip <- bip[bip$`FDR estimate` <= 0.1, ]
genes$BDeQTS <- "0"
genes[genes$Symbol %in% bip$`Gene Symbol`, "BDeQTS"] <- "1"

# Load CDG-2019 eQTS
cdg <- fread(paste0(dea_path, "/3.5.7_diff_expr_exons_for_CDG2-2019_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
cdg <- cdg[cdg$`FDR estimate` <= 0.1, ]
genes$CDGeQTS <- "0"
genes[genes$Symbol %in% cdg$`Gene Symbol`, "CDGeQTS"] <- "1"

# Load MDD-2019 eQTS
mdd <- fread(paste0(dea_path, "/3.5.13_diff_expr_exons_for_MDD3-2019_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
mdd <- mdd[mdd$`FDR estimate` <= 0.1, ]
genes$MDDeQTS <- "0"
genes[genes$Symbol %in% mdd$`Gene Symbol`, "MDDeQTS"] <- "1"

# Load SCZ-2014 eQTS
scz <- fread(paste0(dea_path, "/3.5.14_diff_expr_exons_for_SCZ2-2014_PRS_new_pvalue_and_fdr.csv"), header = T, data.table = F, stringsAsFactors = F)
scz <- scz[scz$`FDR estimate` <= 0.1, ]
genes$SCZeQTS <- "0"
genes[genes$Symbol %in% scz$`Gene Symbol`, "SCZeQTS"] <- "1"

# save data
write.table(genes, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/core_genes_subnetwork_biogrid_groups_overlaps.txt", quote = F, sep = "\t", row.names = F)


### Subnetwork only high throughput & at least 2 experiments with this connection to the core gene ###
# get only high throughput interactions (remove low)
biogrid_high <- biogrid[biogrid$Throughput != "Low Throughput", ]

# check if a connection was found at least twice
subnetwork_core_to_all_ids$`Entrez Gene Interactor C` <- paste(subnetwork_core_to_all_ids$`Entrez Gene Interactor A`, subnetwork_core_to_all_ids$`Entrez Gene Interactor B`, sep = "-")
frequencies <- as.data.frame(table(subnetwork_core_to_all_ids$`Entrez Gene Interactor C`))
subnetwork_core_to_all_ids <- subnetwork_core_to_all_ids[subnetwork_core_to_all_ids$`Entrez Gene Interactor C` %in% frequencies[frequencies$Freq > 1, 'Var1'], ]
biogrid_high$`Entrez Gene Interactor C` <- paste(biogrid_high$`Entrez Gene Interactor A`, biogrid_high$`Entrez Gene Interactor B`, sep = "-")
# get connections to other genes with at least two experiments
subnetwork_core_with_1_inbetween <- biogrid_high[biogrid_high$`Entrez Gene Interactor C` %in% subnetwork_core_to_all_ids$`Entrez Gene Interactor C`, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database")]
# add connections between two core genes (only if high throughput)
subnetwork_core_with_1_inbetween <- merge(subnetwork_core_with_1_inbetween, subnetwork_only_core[subnetwork_only_core$Throughput == "Low Throughput", ], all = T)
# add connections between two non core genes (only if freq > 1)
non_core_genes <- unique(subnetwork_core_with_1_inbetween[!(subnetwork_core_with_1_inbetween$`Entrez Gene Interactor B` %in% core29$entrezID), ]$`Entrez Gene Interactor B`)
subnetwork_non_core_genes <- biogrid_high[biogrid_high$`Entrez Gene Interactor A` %in% non_core_genes & biogrid_high$`Entrez Gene Interactor B` %in% non_core_genes, c("Entrez Gene Interactor A", "Entrez Gene Interactor B", "Official Symbol Interactor A", "Official Symbol Interactor B", "Experimental System", "Experimental System Type", "Author", "Publication Source", "Throughput", "Source Database", "Entrez Gene Interactor C")]
frequencies <- as.data.frame(table(subnetwork_non_core_genes$`Entrez Gene Interactor C`))
subnetwork_non_core_genes <- subnetwork_non_core_genes[subnetwork_non_core_genes$`Entrez Gene Interactor C` %in% frequencies[frequencies$Freq > 1, 'Var1'], ]
subnetwork_core_with_1_inbetween <- merge(subnetwork_core_with_1_inbetween, subnetwork_non_core_genes, all = T)

# get the edges
subnetwork_core_with_1_inbetween_edges <- list()
subnetwork_core_with_1_inbetween_edges <- unname(unlist(apply(subnetwork_core_with_1_inbetween[, c("Official Symbol Interactor A", "Official Symbol Interactor B")], 1, function(x) append(subnetwork_core_with_1_inbetween_edges, x))))

# plot subnetwork of biogrid
g <- make_empty_graph() + vertices(unique(c(subnetwork_core_with_1_inbetween_edges, as.character(core29$Symbols))))
g <- make_empty_graph() + vertices(unique(subnetwork_core_with_1_inbetween_edges))

g <- g + edges(subnetwork_core_with_1_inbetween_edges)
g <- as.undirected(g, mode = "collapse", edge.attr.comb = igraph_opt("edge.attr.comb"))
g$layout <- layout_with_fr
l <- layout.fruchterman.reingold(g) 
node.size= c(10,10,10)
plot(g)

genes <- as.data.frame(unique(subnetwork_core_with_1_inbetween_edges))
colnames(genes) <- "Symbol"
genes$Group <- "A"
genes$Symbol <- as.character(genes$Symbol)
# set core genes to B and check their synonyms (only one gene found!)
genes[genes$Symbol %in% core$Symbols, "Group"] <- "B"
for (i in 1:length(genes$Symbol)) {
  if(genes[i, 2] == "A"){
    synonyms <- strsplit(unique(biogrid[biogrid$`Official Symbol Interactor A` == genes[i, 1], "Synonyms Interactor A"])$`Synonyms Interactor A`, "\\|")
    core_gene <- core[core$Symbols %in% unlist(synonyms), ]
    if(length(core_gene$Symbols) > 0){
      print(genes[i, 1])
      print(core[core$Symbols %in% unlist(synonyms), ])
    }
  }
}
genes[genes$Symbol %in% core29$Symbols, "Group"] <- "C"


write.table(genes, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/core_genes_subnetwork_biogrid_high_throughput_groups.txt", quote = F, sep = "\t", row.names = F)
write.table(subnetwork_core_with_1_inbetween, "/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/06.Index/core_genes_subnetwork_biogrid_high_throughput.txt", quote = F, row.names = F, sep = "\t")

