#######################################
# Monther-infant Study
#
# Figure 3A
#
# Author: HOU Shuwen
#######################################

# Load necessary library
library(readxl)
library(reshape2)
library(ade4)
library(vegan)
library(dplyr)
library(tidyverse)

# data import
setwd("~/Downloads")

# Read Excel file
metadata <- read_excel("Clinical_Meta_all.xlsx", sheet = "Quantile")
host_ID <- read_excel("Clinical_Meta_all.xlsx", sheet = "ID")

# Perform inner join
metadata <- merge(metadata, host_ID, by = "ID", all = FALSE)

# takes a metadata dataframe as input and performs various checks and summaries on it. 
info <- check_metadata(metadata)
check_metadata <- function(metadata, more_missing_values=NULL, unique_rate_thres=0.2){
  if(is.null(more_missing_values)){
    missing_values<-c("not applicable", "Not applicable", "Missing:not collected",
                      "Not provided", "missing: not provided", "unknown",
                      "not provided", "not collected", "NA", NA, "")
  }
  
  check_integers <- function(vector, all=TRUE){
    checker <- grepl("^[0-9]+$", as.character(vector), perl = T)
    total_len <- length(vector)
    n_integers <- sum(checker)
    n_non_integers <- total_len - n_integers
    #cat("Number of integers: ", n_integers, "\n")
    #cat("Number of non-integers: ", n_non_integers, "\n")
    if(all) if_all_integers <- all(checker) ## if any TRUE in the checker
    out <- list(if_all_integers=if_all_integers,
                n_integers=n_integers,
                n_non_integers=n_non_integers)
    return(out)
  }
  
  
  
  n_unique_values <- sapply(metadata, function(x) nlevels(factor(x)))
  unique_values <- sapply(metadata, function(x) paste0(levels(factor(x)), collapse="|"))
  value_balances <- sapply(metadata, function(x) paste0(table(factor(x)), collapse = "|"))
  n_missing_values <- sapply(metadata, function(x) sum(x %in% missing_values))
  n_real_values <- sapply(metadata, function(x) sum(!x %in% missing_values))
  completeness <- n_real_values/nrow(metadata)
  unique_rate <- ifelse(n_real_values==0, 0, n_unique_values/n_real_values)
  #n_unique_real_values <- sapply(metadata, function(x) nlevels(as.factor(x[!x %in% missing_values])))
  n_integers <- sapply(metadata, function(x) check_integers(x[!x %in% missing_values])[[2]])
  all_values_identical <- apply(metadata, 2, function(x) length(unique(x))==1)
  metadata_summ<-data.frame(metadata=colnames(metadata),
                            all_values_identical,
                            total_n=nrow(metadata),
                            n_missing_values,
                            n_real_values,
                            completeness,
                            n_unique_values,
                            unique_values,
                            value_balances,
                            unique_rate)
  
  metadata_summ
}

# convert to numeric, filters the metadata dataframe based on completeness, uniqueness, and data type
metadata <- as.data.frame(lapply(metadata, convert_to_numeric), stringsAsFactors = FALSE)
info_all <- check_metadata(metadata)

convert_to_numeric <- function(x) {
  as.numeric_possible <- suppressWarnings(as.numeric(x))
  ifelse(is.na(as.numeric_possible), x, as.numeric_possible)
}

# subset phyloseq object by sample type
feces <- subset_samples(all, (sequencing == "2bRAD") & (type == "feces"))
meconium <- subset_samples(all, (sequencing == "2bRAD") & (type == "meconium"))
milk <- subset_samples(all, (sequencing == "2bRAD") & (type == "milk"))

# Choose weighted unifrac
distance_feces <-  phyloseq::distance(feces, method = "unifrac")
distance_meconium <-  phyloseq::distance(meconium, method = "unifrac")
distance_milk <-  phyloseq::distance(milk, method = "unifrac")

# Found significant variables to perform adonis
# feces
meta <- as(sample_data(feces), "data.frame") %>% select(-type) %>% select(-sequencing)
merged_data <- merge(meta, metadata, by = "host_ID", all = FALSE)

# reorder to match distance matrix
rownames(merged_data) <- merged_data[,1]
meta <- merged_data[meta$host_ID, ]
meta[is.na(meta)] <- 99
meta <- trim_metadata(meta)

# PERMANOVA
adonis2(distance_feces ~ TC_T2_quan, data = meta)
adonis2(distance_feces ~ CDC, data = meta)

# milk
meta <- as(sample_data(milk), "data.frame") %>% select(-type) %>% select(-sequencing)
merged_data <- merge(meta, metadata, by = "host_ID", all = FALSE)

# reorder to match distance matrix
rownames(merged_data) <- merged_data[,1]
meta <- merged_data[meta$host_ID, ]
meta[is.na(meta)] <- 99
meta <- trim_metadata(meta)

# PERMANOVA
adonis2(distance_milk ~ delivery + cornsBeans, data = meta)

# meconium
meta <- as(sample_data(meconium), "data.frame") %>% select(-type) %>% select(-sequencing)
merged_data <- merge(meta, metadata, by = "host_ID", all = FALSE)

# reorder to match distance matrix
rownames(merged_data) <- merged_data[,1]
meta <- merged_data[meta$host_ID, ]
meta[is.na(meta)] <- 99
meta <- trim_metadata(meta)

# PERMANOVA
# individual
adonis2(distance_meconium ~ deliveryWay, data = meta)
adonis2(distance_meconium ~ abortionNo, data = meta)
adonis2(distance_meconium ~ missed_A, data = meta)
adonis2(distance_meconium ~ pregnancyNo, data = meta)
adonis2(distance_meconium ~ education, data = meta)
adonis2(distance_meconium ~ IL_1b_pg_quan, data = meta)



# cumulative
adonis2(distance_meconium ~ deliveryWay + abortionNo, data = meta)
adonis2(distance_meconium ~ deliveryWay + abortionNo + missed_A, data = meta)
adonis2(distance_meconium ~ deliveryWay + abortionNo + missed_A + pregnancyNo, data = meta)
adonis2(distance_meconium ~ deliveryWay + abortionNo + missed_A + pregnancyNo + education, data = meta)
adonis2(distance_meconium ~ deliveryWay + abortionNo + missed_A + pregnancyNo + education + IL_1b_pg_quan, data = meta)

# adonis2(distance_meconium ~ thyroidfunction, data = meta)  
# adonis2(distance_meconium ~ Insulin_T2_quan, data = meta)       
# adonis2(distance_meconium ~ LCR_T2_quan, data = meta)         
# adonis2(distance_meconium ~ MA_T2_quan, data = meta)         
# adonis2(distance_meconium ~ TC_T2_quan, data = meta)          
        
# ggplot          
# Example data structure
Permanova <- read_excel("sig_variable_info.xlsx", sheet = "significant", skip = 4)
Permanova$variate <- factor(Permanova$variate, levels = unique(Permanova$variate))
Permanova_long <- Permanova %>%
  pivot_longer(cols = c(cumulative_R2, adonis_R2), names_to = "Type", values_to = "R2")
Permanova_long$Type <- factor(Permanova_long$Type, levels = c("cumulative_R2","adonis_R2"), 
                              labels = c("Cumulative","Individual"))

plot <- ggplot(Permanova_long, aes(x = variate, y = R2, fill = Type)) +
  geom_col(position = "dodge") +
  coord_flip() +  # This flips the axes to match your image
  labs(x = "Covariates", y = "Effect Size (R2)", fill = "Variable") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) +
  scale_fill_manual(values = c("Individual" = "gold","Cumulative" = "lightgrey"))
plot
ggsave("./data/plots/3A.png", plot, width = 4, height = 3)


# Try to do strain-level analysis on mother-infant pair
# but no shared species
phyloseq_extract_shared_otus(all, samp_names = c("MST.100.2b", "MEC.100.2b"))
phyloseq_extract_shared_otus(all, samp_names = c("MK.100.2b", "MEC.100.2b"))





