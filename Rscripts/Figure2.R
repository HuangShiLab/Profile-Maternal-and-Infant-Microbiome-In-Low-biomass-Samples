#######################################
# Monther-infant Study
#
# Figure 2
#
# Author: HOU Shuwen
#######################################

# load libraries
library(ggplot2)
library(corrplot)
library(vegan)
library(data.table)
library(readxl)
library(tidyr)
library(dplyr)
library(tidyverse)
library("ape")
library(data.table)
library(vioplot)
library(vegan)
library(phyloseq)
library("patchwork")

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
{
# calculate R square, L2 distance and Bray Curtis distance
r_squared <- numeric()
L2 <- numeric()
bray_curtis <- numeric()
for (i in 1:33) {
  r_squared[i] <- summary(lm(subset1[, i] ~ subset2[, i]))$r.squared
  L2[i] <- vegdist(rbind(subset1[, i],subset2[, i]),method = "euclidean")
  bray_curtis[i] <- vegdist(rbind(subset1[, i],subset2[, i]), method = "bray")
}

# print median values
median(r_squared)
median(L2)
median(bray_curtis)

# for split violin plot
R_violin <- data.frame(r_squared)

# violin plot
data <- data.frame(r_squared = r_squared, L2 = L2, bray_curtis = bray_curtis)
plots_list_1 <- list()
columns_name <- c('r_squared','L2','bray_curtis')
for (i in seq_along(columns_name)) {
  col <- columns_name[i]
  plot <- ggplot(data, aes(x = "", y = .data[[col]])) +
    geom_violin(trim = FALSE, fill = "skyblue") +
    geom_boxplot(width = 0.1, fill = "black") +
    ylim(0, 1) + labs(x = NULL, y = "2bRAD vs WMS (Species)") +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))
  plots_list_1[[col]] <- plot
}


# Merge two feature tables
data <- data.frame(cbind(subset1[,1], subset2[,1]))
zero_rows <- apply(data[, -1, drop = FALSE], 1, function(row) all(row == 0))
data <- data[!zero_rows, ]

# Create a scatter plot with log-transformed data
R2_1 <- paste("R² =", round(r_squared[1], 3))
p1 <- ggplot(data, aes(x = X1, y = X2)) +
  geom_point(size = 3, color = "skyblue", alpha = 0.6) + 
  labs(title = "One Representative Pair",
       x = "Relative Abundance (WMS)",
       y = "Relative Abundance (2bRAD)") +
  geom_smooth(method = "lm", color = "lightblue", fill = "grey90", size = 1, alpha = 1) +
  scale_x_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
  scale_y_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black")) +
  geom_text(aes(x = 0.001, y = 0.0001, label = R2_1),
            color = "black", hjust = 0, vjust = 1, size = 4)
p1
ggsave("./plots/2A.png", p1, width = 4, height = 3)
}

# 16S vs WMS (speces level)
{
  # calculate R square, L2 distance and Bray Curtis distance
  r_squared <- numeric()
  L2 <- numeric()
  bray_curtis <- numeric()
  jaccard <- numeric()
  for (i in 1:33) {
    r_squared[i] <- summary(lm(subset1[, i] ~ subset3[, i]))$r.squared
    L2[i] <- vegdist(rbind(subset1[, i],subset3[, i]),method = "euclidean")
    bray_curtis[i] <- vegdist(rbind(subset1[, i],subset3[, i]), method = "bray")
  }
  
  # print median values
  median(r_squared)
  median(L2)
  median(bray_curtis)
  
  # for split violin plot
  R_violin <- data.frame(r_squared)
  
  # violin plot
  data <- data.frame(r_squared = r_squared, L2 = L2, bray_curtis = bray_curtis)
  plots_list_2 <- list()
  columns_name <- c('r_squared','L2','bray_curtis')
  for (i in seq_along(columns_name)) {
    col <- columns_name[i]
    plot <- ggplot(data, aes(x = "", y = .data[[col]])) +
      geom_violin(trim = FALSE, fill = "gold") +
      geom_boxplot(width = 0.1, fill = "black") +
      ylim(0, 1) + labs(x = NULL, y = "16S vs WMS (Species)") +
      theme(panel.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            axis.line = element_line(size = 0.5, color = "black"))
    plots_list_2[[col]] <- plot
  }

  
  # Merge two feature tables
  data <- data.frame(cbind(subset1[,1], subset3[,1]))
  zero_rows <- apply(data[, -1, drop = FALSE], 1, function(row) all(row == 0))
  data <- data[!zero_rows, ]
  
  # Create a scatter plot with log-transformed data
  R2_1 <- paste("R² =", round(r_squared[6], 3))
  p2 <- ggplot(data, aes(x = X1, y = X2)) +
    geom_point(size = 3, color = "gold", alpha = 0.8) + 
    labs(title = "One Representative Pair",
         x = "Relative Abundance (WMS)",
         y = "Relative Abundance (16S)") +
    geom_smooth(method = "lm", color = "gold", fill = "grey90", size = 1, alpha = 0.7) +
    scale_x_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
    scale_y_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black")) +
    geom_text(aes(x = 0.001, y = 0.0001, label = R2_1),
              color = "black", hjust = 0, vjust = 1, size = 4)
  p2
  ggsave("./plots/2C.png", p2, width = 4, height = 3)
rownames_to_column(var = "RowNames")}

# 16S vs MWS (at genus level)
{# data import
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

r_squared <- numeric()
L2 <- numeric()
bray_curtis <- numeric()

for (i in 1:34) {
  r_squared[i] <- summary(lm(subset1[, i] ~ subset2[, i]))$r.squared
  L2[i] <- vegdist(rbind(subset1[, i],subset2[, i]),method = "euclidean")
  bray_curtis[i] <- vegdist(rbind(subset1[, i],subset2[, i]), method = "bray")
}

# print median values
median(r_squared)
median(L2)
median(bray_curtis)

# for split violin plot
# R_violin[2] <- data.frame(r_squared)

# Draw scatter plot
R2_2 <- paste("R² =", round(r_squared[1], 3))
p3 <- ggplot(merged_feces_genus, aes(x = MST.T2.02.x, y = MST.T2.02.y)) +
  geom_point(size = 3, color = "palegreen3", alpha = 0.8) + 
  labs(title = "One Representative Pair",
       x = "Relative Abundance (WMS)",
       y = "Relative Abundance (16S)") +
  geom_smooth(method = "lm", color = "palegreen3", fill = "lightgrey", size = 1, alpha = 0.7) +
  scale_x_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
  scale_y_log10(limits = c(0.00001, 0.1),labels = scales::percent_format()) +
  theme(panel.background = element_rect(fill = "white", color = NA),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 0.5, color = "black")) +
  geom_text(aes(x = 0.001, y = 0.0001, label = R2_2),
            color = "black", hjust = 0, vjust = 1, size = 4)
p3
ggsave("./plots/2E.png", p3, width = 4, height = 3)

# violin plot
data <- data.frame(r_squared, L2, bray_curtis)
columns_name <- c('r_squared','L2','bray_curtis')
plots_list_3 <- list()
for (i in seq_along(columns_name)) {
  col <- columns_name[i]
  plot <- ggplot(data, aes(x = "", y = .data[[col]])) +
    geom_violin(trim = FALSE, fill = "palegreen3") +
    geom_boxplot(width = 0.1, fill = "black") +
    ylim(0, 1) + labs(x = NULL, y = "16S vs WMS (Genus)") +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))
  plots_list_3[[col]] <- plot
}

r_square <- plots_list_1[[1]] + plots_list_2[[1]] + plots_list_3[[1]] +
  plot_annotation(title = "    All Sample Pairs R-Squared Distribution (N = 33)")
L2 <- plots_list_1[[2]] + plots_list_2[[2]] + plots_list_3[[2]] +
  plot_annotation(title = "    All Sample Pairs Euclidean Distance Distribution (N = 33)")
BC <- plots_list_1[[3]] + plots_list_2[[3]] + plots_list_3[[3]]+
  plot_annotation(title = "    All Sample Pairs BC Dissimilarity Distribution (N = 33)")
}

ggsave("./plots/2B.png", r_square, width = 5, height = 3)
ggsave("./plots/2D.png", L2, width = 5, height = 3)
ggsave("./plots/2F.png", BC, width = 5, height = 3)

# Split violin plot for Poster presentation
{# split violin plot of two methods
library(introdataviz)

# Reshape the data from wide to long format
long_df <- pivot_longer(R_violin, everything(), names_to = "Method", values_to = "R2")

p <- ggplot(long_df,aes(x = "", y = R2, fill = Method))+
  geom_split_violin(alpha = 0.8, trim = FALSE) +
  geom_boxplot(width = .1, fatten = NULL, show.legend = FALSE) +
  theme_minimal() +
  ylim(0, 1) +
  labs(title = "Comparison of Value1 and Value2", x = NULL, y = "Values") +
  scale_fill_manual(values = c("lightblue", "pink"),
                    name = "Method", labels = c("2bRAD", "16S"))
p}