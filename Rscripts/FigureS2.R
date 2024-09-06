#######################################
# Monther-infant Study
#
# Figure S2
#
# Author: HOU Shuwen
#######################################

# load libraries
library(ggplot2)
library(ggpubr)
library(scales)

# load RDS file
setwd("~/Downloads/data")
all <- readRDS("all_phyloseq.rds")

# Species-level comparison
# subset the large phyloseq object
RAD <- subset_samples(all, (sequencing == "2bRAD") & (type == "feces"))
WMS <- subset_samples(all, (sequencing == "WMS") & (type == "feces"))
amplicon <- subset_samples(all, (sequencing == "16S") & (type == "feces"))

# extract feature table
subset1 <- WMS@otu_table %>% data.frame %>% select(sort(names(.)))
subset2 <- RAD@otu_table %>% data.frame %>% select(sort(names(.)))
subset3 <- amplicon@otu_table %>% data.frame %>% select(sort(names(.)))

# 2bRAD vs WMS (speces level)
# Create an empty list to store the plots
plots <- list()

# Loop through each column pair from 1 to 33
for (i in 1:33) {
  data <- data.frame(cbind(subset1[, i], subset2[, i]))
  colnames(data) <- c("X1", "X2")
  name <- gsub(".WMS", "", colnames(subset1)[i])
  
  zero_rows <- apply(data[, -1, drop = FALSE], 1, function(row) all(row == 0))
  data <- data[!zero_rows, ]
  
  p <- ggplot(data, aes(x = X1, y = X2)) +
    geom_point(size = 3, color = "skyblue", alpha = 0.6) + 
    labs(title = name,
         x = "Relative Abundance (WMS)",
         y = "Relative Abundance (2bRAD)") +
    geom_smooth(method = "lm", color = "lightblue", fill = "grey90", size = 1, alpha = 1) +
    scale_x_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
    scale_y_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))
  
  # Store the plot in the list
  plots[[i]] <- p
}
comnine_1 <- ggarrange(plotlist = plots, ncol = 4, nrow = 9)

# Optionally, save the plots to files
ggsave("./plots/S2.png", comnine_1, width = 12, height = 23)


# 16S vs WMS (speces level)
# Create an empty list to store the plots
plots <- list()

# Loop through each column pair from 1 to 33
for (i in 1:33) {
  data <- data.frame(cbind(subset1[, i], subset3[, i]))
  colnames(data) <- c("X1", "X2")
  name <- gsub(".WMS", "", colnames(subset1)[i])
  
  zero_rows <- apply(data[, -1, drop = FALSE], 1, function(row) all(row == 0))
  data <- data[!zero_rows, ]
  
  p <- ggplot(data, aes(x = X1, y = X2)) +
    geom_point(size = 3, color = "gold", alpha = 0.6) + 
    labs(title = name,
         x = "Relative Abundance (WMS)",
         y = "Relative Abundance (2bRAD)") +
    geom_smooth(method = "lm", color = "gold", fill = "grey90", size = 1, alpha = 1) +
    scale_x_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
    scale_y_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))
  
  # Store the plot in the list
  plots[[i]] <- p
}
comnine_2 <- ggarrange(plotlist = plots, ncol = 4, nrow = 9)

# Optionally, save the plots to files
ggsave("./plots/S3.png", comnine_2, width = 12, height = 23)

# 16S vs MWS (at genus level)
# data import
  feces_genus_WMS <- read.table("WMS_2b/feces_genus.txt",header=TRUE)
  feces_genus_16S <- read.table("16S/genus-feces.txt",header=TRUE)
  
  # Identify common columns
  common_columns <- intersect(names(feces_genus_WMS), names(feces_genus_16S))
  feces_genus_WMS <- feces_genus_WMS[, common_columns]
  feces_genus_16S <- feces_genus_16S[, common_columns]
  
  # outer join
  merged_feces_genus <- merge(feces_genus_WMS, feces_genus_16S, by = "Genus", all = TRUE)
  merged_feces_genus[is.na(merged_feces_genus)] <- 0
  
  # Delete rows where every column is zero
  zero_rows <- apply(merged_feces_genus[, -1, drop = FALSE], 1, function(row) all(row == 0))
  merged_feces_genus <- merged_feces_genus[!zero_rows, ]
  
  # Recalculate relative abundance
  sums <- numeric(length = ncol(merged_feces_genus)-1)
  for (i in 2:ncol(merged_feces_genus)) {
    sums[i] <- sum(merged_feces_genus[[i]])
  }
  for (i in 2:ncol(merged_feces_genus)) {
    mutate(merged_feces_genus[i] <- merged_feces_genus[i] / sums[i])
  }
  
  
  # distance
  subset1 <- merged_feces_genus[, 2:35]
  subset2 <- merged_feces_genus[, 36:69]
  
  for (i in 1:33) {
    data <- data.frame(cbind(subset1[, i], subset2[, i]))
    colnames(data) <- c("X1", "X2")
    name <- gsub(".WMS", "", colnames(subset1)[i])
    
    zero_rows <- apply(data[, -1, drop = FALSE], 1, function(row) all(row == 0))
    data <- data[!zero_rows, ]
    
    p <- ggplot(data, aes(x = X1, y = X2)) +
      geom_point(size = 3, color = "palegreen3", alpha = 0.6) + 
      labs(title = name,
           x = "Relative Abundance (WMS)",
           y = "Relative Abundance (2bRAD)") +
      geom_smooth(method = "lm", color = "palegreen3", fill = "grey90", size = 1, alpha = 1) +
      scale_x_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
      scale_y_log10(limits = c(0.00001, 0.1), labels = scales::percent_format()) +
      theme(panel.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.5, color = "black"))
    
    # Store the plot in the list
    plots[[i]] <- p
  }
  comnine_3 <- ggarrange(plotlist = plots, ncol = 4, nrow = 9)
  
  # Optionally, save the plots to files
  ggsave("./plots/S4.png", comnine_3, width = 12, height = 23)
