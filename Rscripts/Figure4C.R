#######################################
# Monther-infant Study
#
# Figure 4C
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
library(ggh4x)

# data import
setwd("~/Downloads")

# OTU table
milk_2_4 <- as.data.frame(read_excel("milk_T2_4.xlsx", sheet = "feature_table"))
rownames(milk_2_4) <- milk_2_4[, 1]
milk_2_4 <- milk_2_4[, -1]
OTU <- otu_table(milk_2_4,taxa_are_rows = TRUE)

# Meta table
META <- as.data.frame(read_excel("milk_T2_4.xlsx", sheet = "meta"))
rownames(META) <- META[, 1]
META <- sample_data(META)  # Prepare metadata as sample data

# TAX table
temp <- rownames(OTU)
TAX <- data.frame(Species = rep("A", length(temp)))
TAX[,1] <- temp
TAX <- tax_table(TAX)
row.names(TAX) <- temp
colnames(TAX) <- c("Species")

# generate phyloseq object
milk_T2_4 <- phyloseq(OTU, META,TAX)

library(metagMisc)
milk_T2_4_filter <- phyloseq_filter_sample_wise_abund_trim(
  milk_T2_4,
  minabund = 0.001,
  relabund = TRUE,
  rm_zero_OTUs = TRUE
)

milk_T2_4_filter <- phyloseq_filter_prevalence(
  milk_T2_4_filter,
  prev.trh = 0.1)


library(microbiomeMarker)
lefse <- run_lefse(milk_T2_4_filter, group = "time",
                   kw_cutoff = 0.01, wilcoxon_cutoff = 0.01,
                   taxa_rank = "Species", lda_cutoff = 3.5)

milk_T24_diff <- lefse@marker_table
milk_T24_diff$feature <- factor(milk_T24_diff$feature, levels = milk_T24_diff$feature)
milk_T24_diff$ef_lda <- milk_T24_diff$ef_lda * ifelse(milk_T24_diff$enrich_group=="T2", -1, 1)

plot2 <- ggplot(milk_T24_diff, aes(x = feature, y = ef_lda, fill = enrich_group)) +
  geom_bar(stat = "identity", position="identity") +
  coord_flip() + 
  scale_fill_manual(values = c("T2" = "skyblue","T4"="palegreen3"))+
  labs(title = "Milk-lactation",
       x = "Species", y = "LDA Value") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line = element_line(size = 0.5, color = "black"),
        axis.text.y = element_text(hjust = ifelse(milk_T24_diff$ef_lda < 0, 1.8, 1.4)))
plot2
ggsave("./data/plots/4C.png", plot2, width = 11, height = 7)
