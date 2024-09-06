#######################################
# Monther-infant Study
#
# Figure 2e
#
# Author: HOU Shuwen
#######################################

# load RDS file
setwd("~/Downloads/data")
all <- readRDS("all_phyloseq.rds")

# load libraries
library(metagMisc)
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)
library(viridis)
library(openxlsx)
library(ggplot2)

# subset phyloseq object
RAD <- subset_samples(all, sequencing == "2bRAD")
OTU <- otu_table(RAD,taxa_are_rows = TRUE)

# TAX table
TAX <- read.table("phyloseq/2b-tax.txt", sep = '\t')
temp <- TAX$V7
temp <- sapply("s__", paste, temp, sep="")
TAX <- tax_table(TAX)
colnames(TAX) <- c("kindom", "phylum", "class", "order", "family", "genus", "species")
rownames(TAX) <- temp

#temp <- rownames(OTU)
#TAX <- data.frame(Species = rep("A", length(temp)))
#TAX[,1] <- temp
## genus
#split_strings <- strsplit(as.character(TAX$Species), "_")
#TAX$Genus <- sapply(split_strings, function(x) x[3])
#rm(split_strings)

# load tax table to phyloseq object
#TAX <- tax_table(TAX)
#row.names(TAX) <- temp
#colnames(TAX) <- c("Species","Genus")
RAD@tax_table <- TAX

# filter keystone taxa
# at least 0.05% abundance (5671->1256)
filter_1 <- phyloseq_filter_sample_wise_abund_trim(
  RAD,
  minabund = 5e-4)

# exist in at least 10% samples (5671->785)
filter_2 <- phyloseq_filter_prevalence(
  RAD,
  prev.trh = 0.1)

# combine the above two (5671->185)
filtered_all <- phyloseq_filter_prevalence(
  filter_1,
  prev.trh = 0.1)

# tax profile (class level)
tax <- data.frame(filtered_all@tax_table)
tip_class <- unique(tax$class)

# find parent node
a=as.numeric(gcc_plot$data[a,1])
class <- data.frame(id=c(191, 226,43,231,263,338,339,355), 
                    type=unique(tax$class))
class_color <- c("Gammaproteobacteria" = "#eacc76", 
                 "Alphaproteobacteria" = "#5c7272",
                 "Desulfovibrionia" = "#718c70", 
                 "Bacteroidia" = "#acab4b",
                 "Clostridia" = "#7c99bc", 
                 "Negativicutes" = "#a5b3c1",
                 "Bacilli" = "#c7a8a3", 
                 "Actinomycetia" = "#e9b962")
# ggtree skeleton
gcc_plot = ggtree(filtered_all@phy_tree, layout="fan", open.angle = 2, size = 0.15, alpha = 1) +
  geom_tiplab(linesize=.15, mapping = aes(label = NA)) +
  geom_highlight(data=class,mapping=aes(node=id, fill=type),alpha=0.8) +
  scale_fill_manual(values = class_color, guide = "none") + 
  coord_polar(theta = 'y', start = 0, direction = -1) +
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(0,0,0,0))
gcc_plot
ggsave("./plots/tree.png", gcc_plot, width = 4, height = 4)


# first layer: prevalence in different sample types
feces <- subset_samples(filtered_all, type == "feces")
milk <- subset_samples(filtered_all, type == "milk")
meconium <- subset_samples(filtered_all, type == "meconium")
t1 <- microbiome::prevalence(feces)
t2 <- microbiome::prevalence(milk)
t3 <- microbiome::prevalence(meconium)
prevalence <- data.frame(
  Feces = t1,
  Milk = t2,
  Meconium = t3)

prevalence_long <- prevalence %>%
  rownames_to_column("taxa") %>%
  gather(key = "type", value = "Prevalence", -taxa)

gcc_plot_data <- gcc_plot$data %>% filter(isTip)  %>%
  left_join(prevalence_long, by = c("label" = "taxa"))

gcc_plot_data <- gcc_plot_data%>% filter(Prevalence > 0)
gcc_plot_data$type_numeric <- as.numeric(
  factor(gcc_plot_data$type, levels = c("Feces", "Milk", "Meconium")))

color_type <- c("Feces" = "chocolate3", "Milk" = "gold", "Meconium" = "palegreen3")
sample_types <- data.frame(
  type = c("Feces", "Milk", "Meconium"),
  x = 1,
  y = c(3, 2, 1)
)

label_plot <- ggplot(sample_types, aes(x = x, y = y)) +
  geom_point(aes(color = type), size = 10) +
  scale_color_manual(values = color_type) +
  geom_text(aes(label = type), nudge_x = 0.1, hjust = 0, size = 5) +
  theme_void() +
  theme(legend.position = "none") +
  scale_y_continuous(limits = c(0.5, 3.5), breaks = c(1, 2, 3)) +
  coord_cartesian(clip = "off")
label_plot

layer2 <- ggplot() +
  geom_point(data = gcc_plot_data, 
             aes(x = y, y = type_numeric, size = Prevalence, 
                 fill = type, color = type), shape = 21) +
  scale_fill_manual(values = color_type) +
  scale_color_manual(values = color_type) +
  scale_size_continuous(range = c(0, 1)) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_line(size = 0.15),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0, 0, 0, 0),
        panel.border = element_blank()) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, max(gcc_plot_data$y, na.rm = TRUE) + 1) +
  scale_y_continuous(breaks = c(1, 2, 3), limits = c(-20, 3))
layer2
ggsave("./plots/layer1.png", layer2, width = 4, height = 4)

# second layer: taxonomic information
tax <- as.data.frame(filtered_all@tax_table)
tax$taxa <- rownames(tax)
gcc_plot_data <- gcc_plot$data %>% filter(isTip)%>%
  left_join(tax, by = c("label" = "taxa"))
genus_color <- viridis(82)
genus_color <- setNames(genus_color[1:82], unique(genus$Genus))

# Create a separation indicator
gcc_plot_data <- gcc_plot_data %>%
  arrange(y) %>%
  mutate(change = c(TRUE, diff(as.numeric(factor(Genus))) != 0))

layer1 <- ggplot() + 
  geom_tile(data = gcc_plot_data, aes(x = y, y = 2, fill = Genus), color = "transparent") +
  scale_fill_manual(values = genus_color) +
  # Add white lines between different genera to create separation
  geom_segment(data = gcc_plot_data %>% filter(change), 
               aes(x = y - 0.5, xend = y - 0.5, y = 1.5, yend = 2.5), 
               color = "white", size = 0.5) +
  theme_minimal() +
  theme(rect = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.grid = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        plot.margin = margin(0,0,0,0)) +
  coord_polar(start = 0, direction = -1) +
  xlim(0, max(gcc_plot_data$y, na.rm = TRUE) + 1) +
  scale_y_continuous(breaks = NULL, limits = c(-15,3))

# annotate genus name
gcc_plot_data_processed <- gcc_plot_data %>%
  mutate(previous_genus = lag(Genus)) %>% # Create a column for the previous genus
  filter(Genus != previous_genus | is.na(previous_genus)) %>% # Keep only rows where the genus is different from the previous row
  select(-previous_genus)

write.xlsx(gcc_plot_data_processed, "genus.xlsx")
gcc_plot_data_processed <- read.xlsx("genus.xlsx")

gcc_plot_data_processed <- gcc_plot_data_processed %>%
  mutate(angle = (y / max(gcc_plot_data$y, na.rm = TRUE) * 360) + 90)

layer1 <- layer1 +
  geom_text(data = gcc_plot_data_processed, 
            aes(x = y, y = 3, label = Genus, angle = angle),
            hjust = 0, vjust = 0.5, size = 3)
layer1
ggsave("./plots/layer2.png", layer1, width = 9, height = 9)
