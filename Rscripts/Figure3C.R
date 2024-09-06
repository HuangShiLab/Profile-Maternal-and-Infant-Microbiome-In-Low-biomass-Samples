#######################################
# Monther-infant Study
#
# Figure 3C
#
# Author: HOU Shuwen
#######################################

# delivery in meconium
# need tax table
temp <- rownames(as.data.frame(meconium@otu_table))
TAX <- data.frame(Species = rep("A", length(temp)))
TAX[,1] <- gsub("s__", "", temp)
TAX <- tax_table(TAX)
row.names(TAX) <- temp
colnames(TAX) <- c("Species")
meconium@tax_table <- TAX

# differential abundance
library(microbiomeMarker)
meco_filter <- subset_samples(meconium, (delivery == "eutocia") | (delivery == "cesarean"))
meco_filter <- phyloseq_filter_sample_wise_abund_trim(
  meco_filter,
  minabund = 0.001,
  relabund = TRUE,
  rm_zero_OTUs = TRUE
)

# clr transformation
meco_transformed <- microbiome::transform(meco_filter, 'clr')

lefse <- run_lefse(meco_filter, group = "delivery", transform = "log10p",
                   kw_cutoff = 0.01, wilcoxon_cutoff = 0.01,
                   taxa_rank = "Species", lda_cutoff = 3)
# aldx <- run_aldex(meco_transformed, group = "delivery")
ancom <- run_ancom(meco_filter, norm = "none", group = "delivery")

meco_diff <- lefse@marker_table
meco_diff_ancom <- ancom@marker_table

plot1 <- ggplot(meco_diff, aes(x = feature, y = ef_lda, fill = enrich_group)) +
  geom_bar(stat = "identity") +
  coord_flip() + 
  scale_fill_manual(values = c("cesarean" = "skyblue", "eutocia" = "palegreen3")) +
  labs(x = "Species", y = "LDA Value") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, color = "black"))

# starchy vegetables in milk
# need tax table
temp <- rownames(as.data.frame(milk@otu_table))
TAX <- data.frame(Species = rep("A", length(temp)))
TAX[,1] <- gsub("s__", "", temp)
TAX <- tax_table(TAX)
row.names(TAX) <- temp
colnames(TAX) <- c("Species")
milk@tax_table <- TAX


# filter and transformation
library(metagMisc)

milk_transformed <- microbiome::transform(milk, 'clr')
milk_filter <- phyloseq_filter_sample_wise_abund_trim(
  milk,
  minabund = 0.001,
  relabund = TRUE,
  rm_zero_OTUs = TRUE
)

lefse <- run_lefse(milk_filter, group = "cornbeans", taxa_rank = "Species", lda_cutoff = 3)
# aldx <- run_aldex(milk_transformed, group = "cornbeans")
ancom <- run_ancom(milk_filter, group = "cornbeans")

# unload package
detach("package:microbiomeMarker", unload = TRUE)

milk_diff <- lefse@marker_table

# Lefse plot
plot2 <- ggplot(milk_diff, aes(x = feature, y = ef_lda, fill = enrich_group)) +
    geom_bar(stat = "identity") +
    coord_flip() + 
    scale_fill_manual(values = c("Sometimes" = "gold"))+
    labs(title = "Milk-Starchy veges",
         x = "Species", y = "LDA Value") +
    theme_minimal()
  
library(ggpubr)
combined_plots <- ggarrange(plot1, plot2, ncol = 2, nrow = 1)
  
# box plot to show differential abundance
meco_filter <- subset_samples(meconium, (delivery == "eutocia") | (delivery == "cesarean"))
meco <- as.data.frame(meco_filter@otu_table)
meta <- as.data.frame(meco_filter@sam_data)

species_names <- c("s__Escherichia_coli_D",
                   "s__Escherichia_dysenteriae",
                   "s__Rothia_sp902373285",
                   "s__Pseudomonas_aeruginosa",
                   "s__Staphylococcus_capitis",
                   "s__Staphylococcus_epidermidis")
species_names_cleaned <- gsub("s__", "", species_names)
plot_list <- list()

# Loop through the species names and bind the corresponding rows
for (i in seq_along(species_names)) {
  species <- species_names[i]
  # build data frame for abundance
  abundance <- rbind(meco[species, ], meta$"delivery")
  rownames(abundance) <- c("Abundance","Delivery")
  # data transformation
  abundance <- as.data.frame(t(abundance))
  abundance$Abundance <- as.numeric(abundance$Abundance)
  abundance$Abundance <- log10(abundance$Abundance + 1)
  abundance$Delivery <- as.factor(abundance$Delivery)
  abundance <- remove_outliers(abundance, "Abundance","Delivery")
  # plot
  plot_list[[i]] <- ggplot(abundance, aes(x = Delivery, y = Abundance, fill = Delivery)) +
    geom_boxplot() + geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(alpha = 0.5, color = "black",size = 1, cex = 0.5) +
    scale_fill_manual(values = c("skyblue","palegreen3"),
                      labels = c("cesarean" = "Cesarean", "eutocia" = "Vaginal")) +
    labs(title = paste0(species_names_cleaned[i]), x = "",
         y = "log(Relative Abundance + 1)", color = "Delivery Mode") +
    theme(panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"),
          axis.text.x = element_blank(),
          plot.title = element_text(face = "italic"))
}
combined_plots <- ggarrange(plotlist = plot_list, ncol = 3, nrow = 2, 
                            common.legend = TRUE, legend = "right")
combined_plots
ggsave("./plots/3C.png", combined_plots, width = 10, height = 4.5)

# Define the function to remove outliers
remove_outliers <- function(df, value_col, group_col) {
  df_clean <- df %>%
    group_by(!!sym(group_col)) %>%
    mutate(
      Q1 = quantile(!!sym(value_col), 0.25),
      Q3 = quantile(!!sym(value_col), 0.75),
      IQR = Q3 - Q1,
      Lower = Q1 - 1.5 * IQR,
      Upper = Q3 + 1.5 * IQR
    ) %>%
    filter(!!sym(value_col) >= Lower & !!sym(value_col) <= Upper) %>%
    ungroup() %>%
    select(-Q1, -Q3, -IQR, -Lower, -Upper)
  
  return(df_clean)
}

