#######################################
# Monther-infant Study
#
# Figure 1B
#
# Author: HOU Shuwen
#######################################

# load libraries
library(ggplot2)       # For data visualization
library(ggpubr)        # For combining and formatting plots
library(vegan)         # For diversity analysis
library(phyloseq)      # For handling microbiome data
library(data.table)     # For efficient data manipulation
library(readxl)        # For reading Excel files
library(dplyr)         # For data manipulation
library(tidyverse)     # For data manipulation and visualization
library("ape")         # For phylogenetic analysis

# data import
setwd("~/Downloads/data")

# load RDS file
all <- readRDS("all_phyloseq.rds")

# calculate distance
distance_methods <- c("unifrac", "wunifrac", "bray")  # Define distance methods
distances <- lapply(distance_methods, function(method) {
  phyloseq::distance(all, method = method)  # Calculate distances using each method
})

# generate PCoA plot
plots_list <- list()  # Initialize a list to store plots
for (i in seq_along(distance_methods)) {
  # Calculate ordination
  pcoa_results <- ordinate(all, method = "PCoA", distance = distances[[i]])  # Perform PCoA
  # Make Plot
  pcoa_plot <- NULL
  pcoa_plot <- plot_ordination(all, pcoa_results, color = 'type', shape = 'sequencing') + 
    geom_point(size = 2, alpha = 0.8) +
    scale_shape_manual(values = c(0,1,2)) +
    scale_color_manual(values = c("chocolate3","palegreen3", "gold")) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))
  plots_list[[i]] <- pcoa_plot  # Store the plot in the list
}
combined_plots <- ggarrange(plotlist = plots_list, ncol = 3, nrow = 1, 
                            common.legend = TRUE, legend = "right")
combined_plots
ggsave("./plots/1B.png", combined_plots, width = 11, height = 3)

# PERMANOVA
library(pairwiseAdonis)
meta <- data.frame(all@sam_data)
pairwise.adonis(distances[[1]], meta$type)
pairwise.adonis(distances[[1]], meta$sequencing)
adonis2(distances[[1]] ~ type, data = meta)
adonis2(distances[[1]] ~ sequencing, data = meta)
