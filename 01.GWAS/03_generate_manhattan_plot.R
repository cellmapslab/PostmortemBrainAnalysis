############
### GWAS ###
############
# In this R script, we will generate the manhattan plot for the GWAS results of our dataset.

# 1. Load Packages for own R installation
.libPaths("/home/icb/karolina.worf/R/x86_64-pc-linux-gnu-library/3.5/")
# data transformation packages
library(data.table)
library(dplyr)
# color palette
library(viridis) # for colorblind color palette
# packages for plotting
library(qqman)
library(ggplot2)
library(svglite)

pdf.options(useDingbats = TRUE)

# 2. Load Data
gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc.logistic", header = T, stringsAsFactors = F, data.table = F)
gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values

# 3. Generate Manhattan Plot for the GWAS (logistic regression) with covariates
# 3.1 The easy way using the qqman package:
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_logistic_with_covariates.pdf", width = 17, height = 7)
manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P", suggestiveline = -log10(1e-02))
dev.off()

# 3.2 The more complicated way using ggplot2:
#  First of all, we need to compute the cumulative position of SNP.
don <- gwasResults %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len = max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResults, ., by = c("CHR" = "CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
axisdf = don %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Ready to make the plot using ggplot2:
#png("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_logistic_with_covariates2.png", type = "cairo", width = 1700, height = 700)
pdf("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_logistic_with_covariates2.pdf", width = 17, height = 7)
ggplot(don, aes(x = BPcum, y = -log10(P))) +
  # Show all points
  geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
  scale_color_manual(values = rep(c("black", "darkgrey"), 11)) + #c(viridis(5)[2], viridis(5)[4]), 11)) +
  # custom x-axis:
  scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +  # remove space between plot area and x axis
  # add genome-wide significant and suggestive lines
  geom_hline(yintercept = -log10(5e-8), color = viridis(2)[1], size = 2) +  # genome-wide significant threshold line
  geom_hline(yintercept = -log10(0.01), color = viridis(2)[2], size = 2) +  # nominal threshold line of p < 0.01
  # Custom the theme:
  theme_classic() +
  xlab("Chromosome") +
  theme( 
    legend.position="none",
    plot.margin = margin(4, 2, 1, 1, unit = "mm"),
    axis.line.x = element_blank(),
    axis.line.y = element_line(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text = element_text(size = 30),
    axis.title = element_text(size = 36)
  )
dev.off()

# 4. Generate Manhattan Plot for the GWAS without covariates
# Load data
#gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc", header = T, stringsAsFactors = F, data.table = F)
#gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values

# 4.1 The easy way using the qqman package:
#bitmap("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_assoc_without_covariates.png", "png16m")
#manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P")
#dev.off()

# 4.2 The more complicated way using ggplot2:
#  First of all, we need to compute the cumulative position of SNP.
#don <- gwasResults %>% 
  # Compute chromosome size
#  group_by(CHR) %>% 
#  summarise(chr_len = max(BP)) %>% 
  # Calculate cumulative position of each chromosome
#  mutate(tot = cumsum(chr_len) - chr_len) %>%
#  select(-chr_len) %>%
  # Add this info to the initial dataset
#  left_join(gwasResults, ., by = c("CHR" = "CHR")) %>%
  # Add a cumulative position of each SNP
#  arrange(CHR, BP) %>%
#  mutate(BPcum = BP + tot)

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
#axisdf = don %>% group_by(CHR) %>% summarize(center = (max(BPcum) + min(BPcum)) / 2)

# Ready to make the plot using ggplot2:
#ggplot(don, aes(x = BPcum, y = -log10(P))) +
  # Show all points
#  geom_point( aes(color = as.factor(CHR)), alpha = 0.8, size = 2) +
#  scale_color_manual(values = rep(c("black", "darkgrey"), 11)) + # rep(c("darkgrey", "grey"), 11)) + #c(viridis(5)[2], viridis(5)[4]), 11)) +
  # custom X axis:
#  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand = c(0, 0)) +
#  scale_y_continuous(expand = c(0, 0), limits = c(0, 8)) +  # remove space between plot area and x axis
  # add genome-wide significant and suggestive lines
#  geom_hline(yintercept = -log10(5e-8)) +  # genome-wide significant threshold line
#  geom_hline(yintercept = -log10(0.01), linetype = "dashed") +  # nominal threshold line of p < 0.01
  # Custom the theme:
#  theme_classic() +
#  xlab("Chromosome") +
#  theme( 
#    legend.position="none", 
#    plot.margin = margin(4, 2, 1, 1, unit = "mm"),
#    axis.line.x = element_blank(),
#    axis.line.y = element_line(),
#    panel.border = element_blank(),
#    panel.grid.major.x = element_blank(),
#    panel.grid.minor.x = element_blank(),
#    axis.text = element_text(size = 26),
#    axis.title = element_text(size = 30)
#  )
#ggsave("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_assoc_without_covariates2.pdf", width = 12)

