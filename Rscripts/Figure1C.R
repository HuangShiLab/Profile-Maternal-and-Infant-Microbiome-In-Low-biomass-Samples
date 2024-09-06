#######################################
# Monther-infant Study
#
# Figure 1C
#
# Author: HOU Shuwen
#######################################

library(btools)
library(metagMisc)
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)
library(rstatix)

# subset feces sample
feces <- subset_samples(all, (type == "feces"))
feces <- phyloseq_filter_taxa_rel_abund(
  feces,  frac = 1e-04)
meco <- subset_samples(all, (type == "meco"))

# have a glance at alpha diversity
measures <- c("Shannon", "Simpson")
plot_richness(feces, x="sequencing", measures = measures)
# pd <- estimate_pd(feces)

# Perform the statistical test for alpha diversity
alpha_diversity <- estimate_richness(all, measures = measures)
alpha_diversity_feces <- estimate_richness(feces, measures = measures)
# alpha_diversity[,3] <- pd[,1]

df_feces <- merge(alpha_diversity_feces, feces@sam_data, by = "row.names")
df <- merge(alpha_diversity, all@sam_data, by = "row.names")
names(df)[1] <- "Sample"
df$sequencing <- factor(df$sequencing, levels = c("16S","WMS","2bRAD"))
df$type <- factor(df$type, levels = c("feces","milk","meconium"))

# Box plots
plots <- list()
for (measure in measures) {
  # test statistics
  stat_test <- df_feces %>%
    wilcox_test(as.formula(paste(measure, "~ sequencing"))) %>%
    add_significance()  %>% add_xy_position(x = "sequencing")
  stat_test <- stat_test[-2, ]
  # box plots
  plot <- ggplot(df, aes_string(x = "interaction(sequencing, type)", y = measure)) + 
    geom_boxplot(aes(fill = type), outlier.shape = NA) +
    scale_fill_manual(values = c("chocolate3", "gold", "palegreen3")) +
    geom_beeswarm(alpha = 0.6, color = "black") +
    stat_pvalue_manual(stat_test, label = "p.adj.signif", tip.length = 0.02) +  
    labs(y = paste(measure, "diversity"), x = "") + 
    theme(legend.position = "none",
          panel.background = element_rect(fill = "white", color = NA),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(size = 0.5, color = "black"))+
    scale_x_discrete(labels = function(x) sapply(strsplit(x, "\\."), `[`, 1))
  # Add the plot to the list
  plots[[measure]] <- plot
}
alpha_diversity_plot <- ggarrange(plotlist = plots, ncol = 2, nrow = 1,
                                  common.legend = TRUE, legend = "right")
alpha_diversity_plot
ggsave("./plots/1C.png", alpha_diversity_plot, width = 8, height = 3)
