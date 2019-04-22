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

# 2. Load Data
gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc.logistic", header = T, stringsAsFactors = F, data.table = F)
gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values

# 3. Generate Manhattan Plot for the GWAS (logistic regression) with covariates
# 3.1 The easy way using the qqman package:
bitmap("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_logistic_with_covariates.png", "png16m")
manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P")
dev.off()

# 3.2 The more complicated way using ggplot2:
#  First of all, we need to compute the cumulative position of SNP.
don <- gwasResults %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2)

# Ready to make the plot using ggplot2:
png("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_logistic_with_covariates2.png", type="cairo")
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

# 4. Generate Manhattan Plot for the GWAS without covariates
# Load data
gwasResults <- fread("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/GWAS_results/pmb_info06_all_chr_maf0.01_merged_GEX_filtered_GWAS.assoc", header = T, stringsAsFactors = F, data.table = F)
gwasResults <- gwasResults[!is.na(gwasResults$P), ]  # remove NA values

# 4.1 The easy way using the qqman package:
bitmap("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_assoc_without_covariates.png", "png16m")
manhattan(gwasResults, chr = "CHR", bp = "BP", snp = "SNP", p = "P")
dev.off()

# 4.2 The more complicated way using ggplot2:
#  First of all, we need to compute the cumulative position of SNP.
don <- gwasResults %>% 
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(gwasResults, ., by=c("CHR"="CHR")) %>%
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot)

# Then we need to prepare the X axis. Indeed we do not want to display the cumulative position of SNP in bp, but just show the chromosome name instead.
axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2)

# Ready to make the plot using ggplot2:
png("/storage/groups/ccm01/workspace/MPIP-Dean_postmortem_brain/01.GWAS/figures/manhattan_plot_GWAS_assoc_without_covariates2.png", type="cairo")
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

