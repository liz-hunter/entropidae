library(tidyverse)
library(vegan)
library(readr)
library(taxonomizr)
library(gtools)

# something weird is going on here bc there are 6255 samples in virid but we sampled 6257, but it's a later problem
# maybe use the euk taxonomy or refseq virid to try and chase down the 2 entries
#acc_sp <- read_tsv("/projects/handy_grp/liz_working/entropy/data_summaries/db_accessions.tsv")
#sample <- unique(unlist(acc_sp))
#ref <- all_virid %>% select('Assembly Accession')
#ref <- unique(ref$`Assembly Accession`)
#missing_values <- setdiff(ref, sample)

# load in database tsv
taxid_sp <- read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/keys/db_comp/combined_taxid_symlinklist.tsv")
all_virid = read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/keys/viridiplantae_all_ncbi20250303.tsv", col_names=TRUE) 
sql_path <- "/projects/handy_grp/taxonomizr20250305/accessionTaxa.sql"

# grab higher level taxonomy
taxid_gen <- as.data.frame(lapply(taxid_sp, function(col) {
  taxonomy <- getTaxonomy(col, sql_path, desiredTaxa = "genus", getNames = FALSE)
  taxonomy[, "genus"]
}))

taxid_fam <- as.data.frame(lapply(taxid_sp, function(col) {
  taxonomy <- getTaxonomy(col, sql_path, desiredTaxa = "family", getNames = FALSE)
  taxonomy[, "family"]
}))

# create a list
taxid_levels <- list(
  species = taxid_sp,
  genus = taxid_gen,
  family = taxid_fam
)

# calculate basic metrics
calc_diversity <- function(df, level_name) {
  lapply(names(df), function(colname) {
    vec <- df[[colname]]
    freq_table <- table(vec)
    
    data.frame(
      Database = colname,
      Level = level_name,
      Shannon = diversity(freq_table, index = "shannon"),
      Simpson = diversity(freq_table, index = "simpson"),
      Richness = length(unique(vec))
    )
  }) %>% bind_rows()
}

# Combine metrics for all levels
diversity_df <- bind_rows(
  lapply(names(taxid_levels), function(level) {
    calc_diversity(taxid_levels[[level]], level)
  })
)

diversity_long <- diversity_df %>%
  pivot_longer(cols = c("Shannon", "Simpson", "Richness"),
               names_to = "Metric", values_to = "Value")

diversity_long$Database <- factor(diversity_long$Database,
                                  levels = mixedsort(unique(diversity_long$Database)))

diversity_long$Level <- factor(diversity_long$Level,
                               levels = c("species", "genus", "family"))

nice_colors <- c(
  species = "#E69F00",
  genus   = "#56B4E9",
  family  = "#009E73"
)

ggplot(diversity_long, aes(x = Database, y = Value, fill = Level)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = nice_colors) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.title.position = "plot"   
  ) +
  labs(
    title = "Diversity Metrics by Database & Taxonomic Level",
    y = "Diversity Value",
    fill = "Rank"
  )

# compare the databases to the overarching virid dataset
sample_stat <- function(master, samples_df) {
  # Master taxid frequencies (with duplicates)
  master_taxids <- master$`Organism Taxonomic ID`
  master_freq <- as.data.frame(table(master_taxids), stringsAsFactors = FALSE)
  colnames(master_freq) <- c("taxid", "Master_Freq")
  
  # Initialize storage
  taxid_stats_list <- list()
  summary_stats_list <- list()
  
  for (sample_name in colnames(samples_df)) {
    sample_taxids <- samples_df[[sample_name]]
    sample_taxids <- sample_taxids[!is.na(sample_taxids)]
    
    # Frequency of taxids in sample
    sample_freq <- as.data.frame(table(sample_taxids), stringsAsFactors = FALSE)
    colnames(sample_freq) <- c("taxid", "Sample_Freq")
    
    # Merge master and sample
    merged <- merge(master_freq, sample_freq, by = "taxid", all = TRUE)
    merged[is.na(merged)] <- 0
    
    # Add proportions
    merged$Proportion_In_Master <- merged$Master_Freq / sum(master_freq$Master_Freq)
    merged$Proportion_In_Sample <- merged$Sample_Freq / sum(sample_freq$Sample_Freq)
    
    # Store per-sample taxid stats
    taxid_stats_list[[sample_name]] <- merged
    
    # Summary stats for this sample
    sample_stats <- data.frame(
      sample = sample_name,
      total_taxids_in_sample = length(sample_taxids),
      unique_taxids_in_sample = length(unique(sample_taxids)),
      total_taxids_in_master = length(master_taxids),
      unique_taxids_in_master = length(unique(master_taxids)),
      richness_prop = length(unique(sample_taxids)) / length(unique(master_taxids)),
      percent_taxids_shared = mean(merged$Sample_Freq > 0)  # proportion of master taxids present in sample
    )
    
    summary_stats_list[[sample_name]] <- sample_stats
  }
  
  # Combine summary stats into one df
  summary_stats_df <- do.call(rbind, summary_stats_list)
  
  return(list(
    taxid_stats = taxid_stats_list,
    summary_stats = summary_stats_df
  ))
}

metrics <- sample_stat(all_virid, taxid_sp)
metrics$summary_stats

# samples reduc comes out of this function
process_samples <- function(samples_df, accession_df = NULL) {
  # List of taxids per sample
  samples_reduc <- lapply(samples_df, function(col) col[!is.na(col)])
  
  # Optional: list of accessions per sample
  samples_acc <- NULL
  if (!is.null(accession_df)) {
    samples_acc <- lapply(accession_df, function(col) col[!is.na(col)])
  }
  
  # Combine all taxids from all samples and compute frequency
  all_taxids <- unlist(samples_reduc)
  sample_freq <- as.data.frame(table(all_taxids), stringsAsFactors = FALSE)
  colnames(sample_freq) <- c("Value", "Sampled_Count")
  
  # Return result in the expected format
  return(list(
    samples_reduc = samples_reduc,
    samples_acc = samples_acc,
    freq_data = sample_freq
  ))
}

# function to grab unique taxid counts for each sample
get_histogram_data <- function(samples_reduc, category) {
  counts <- lapply(samples_reduc, function(x) length(unique(x)))
  counts_vector <- unlist(counts)
  histogram_data <- data.frame(UniqueCounts = counts_vector, Category = category)
  return(histogram_data)
}

# Stats 
freq <- virid_results$freq_data

cor_test <- cor.test(freq$Original_Count, freq$Sampled_Count)
cor_value <- cor_test$estimate
p_value <- cor_test$p.value

# Scatter plot of original count vs sampled count
ggplot(freq, aes(x = Original_Count, y = Sampled_Count)) +
  geom_jitter(color = "#FF3CC7", alpha = 0.6, width = 0.2, height = 500) +  
  geom_smooth(method = "lm", se = FALSE, color = "#4C1A57", linewidth = 1) +  
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 10000)) +  
  labs(
    title = "Original Frequency vs. Sampled Frequency",
    x = "Frequency in Original Data",
    y = "Frequency in Sampled Data"
  ) +
  annotate("text", x = 10, y = 9000, 
           label = paste("r =", round(cor_value, 4), "\np =", format.pval(p_value)), 
           color = "#4C1A57", size = 3.5, hjust = 0) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Optionally, grab accessions that are common to all samples 
common_accessions <- reduce(grab, intersect) %>% as.data.frame()
colnames(common_accessions) <- "accessions"

# Grab the count of unique taxids in each sample 
counts <- lapply(virid_results$samples_reduc, function(x) length(unique(x)))
counts_vector <- unlist(counts)

# Make a quick histogram to visualize the data distribution
ggplot(data.frame(UniqueCounts = counts_vector), aes(x = UniqueCounts)) +
  geom_histogram(binwidth = 10, fill = "#4C1A57", color = "black", alpha = 0.7) +
  labs(title = "Histogram of Unique Taxid Counts", x = "Unique Taxid Count", y = "Frequency") +
  theme_bw()

# Define colors for categories
category_colors <- c("Archaea" = "#FF3CC7", 
                     "Bacteria" = "#2C8CFF", 
                     "Viridiplantae" = "#57FF2C", 
                     "Virus" = "#4C1A57", 
                     "Eukaryota" = "#FFC300")

# Combine results into one data frame
combined_freq_data <- bind_rows(
  mutate(archea_results$freq_data, Category = "Archaea"),
  mutate(bacteria_results$freq_data, Category = "Bacteria"),
  mutate(virid_results$freq_data, Category = "Viridiplantae"),
  mutate(virus_results$freq_data, Category = "Virus"),
  mutate(euk_results$freq_data, Category = "Eukaryota")
)

# Scatter plot
ggplot(combined_freq_data, aes(x = Original_Count, y = Sampled_Count)) +
  geom_jitter(aes(color = Category), alpha = 0.6, width = 0.2, height = 500) +
  geom_smooth(aes(color="black"), method = "lm", se = FALSE, linewidth = 1) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 10000)) +
  labs(
    title = "Original Frequency vs. Sampled Frequency by Category",
    x = "Frequency in Original Data",
    y = "Frequency in Sampled Data"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_color_manual(values = category_colors) +
  facet_wrap(~ Category, scales = "free")  # Create separate panels for each category

# Define colors for categories
category_colors <- c("Archaea" = "#FF3CC7", 
                     "Bacteria" = "#2C8CFF", 
                     "Viridiplantae" = "#57FF2C", 
                     "Virus" = "#4C1A57", 
                     "Eukaryota" = "#FFC300")

# Collect histogram data for each dataset
histogram_archea <- get_histogram_data(archea_results$samples_reduc, "Archaea")
histogram_bacteria <- get_histogram_data(bacteria_results$samples_reduc, "Bacteria")
histogram_virid <- get_histogram_data(virid_results$samples_reduc, "Viridiplantae")
histogram_virus <- get_histogram_data(virus_results$samples_reduc, "Virus")
histogram_euk <- get_histogram_data(euk_results$samples_reduc, "Eukaryota")

# Combine all histogram data into one data frame
combined_histogram_data <- bind_rows(
  histogram_archea,
  histogram_bacteria,
  histogram_virid,
  histogram_virus,
  histogram_euk
)

# Calculate the mean number of unique taxa per sample
means <- combined_histogram_data %>% group_by(Category) %>%  
  summarise(mean = mean(UniqueCounts)) 

# Calculate the total number of taxids in each dataset, and the total number of unique taxids in each dataset
annotation_data <- data.frame(
  Category = c("Eukaryota", "Viridiplantae", "Bacteria", "Archaea", "Virus"),
  n = c(nrow(euk), nrow(virid), nrow(bac), nrow(archea), nrow(virus)),
  nuniq = c(length(unique(euk$`Organism Taxonomic ID`)), 
            length(unique(virid$`Organism Taxonomic ID`)),
            length(unique(bac$`Organism Taxonomic ID`)),
            length(unique(archea$`Organism Taxonomic ID`)),
            length(unique(virus$`Organism Taxonomic ID`))),
  x = c(16300, 1915, 99300, 2754, 59670)) %>% #define horizontal placement of labels
  mutate(label1 = paste0("n: ", n)) %>%
  mutate(label2 = paste0("uniq: ", nuniq))

# Calculate the overall proportion of the database that is unique (unique taxa/overall n taxa)
annotation_data <- merge(annotation_data, means, by="Category")
annotation_data$label3 <- paste0("prop: ", round(annotation_data$nuniq/annotation_data$n, digits=2))
annotation_data$label4 <- paste0("mean: ", annotation_data$mean)

# Plot the giant histogram with all the annotation data
ggplot(combined_histogram_data, aes(x = UniqueCounts, fill = Category)) +
  geom_histogram(binwidth = 20, color = "black", alpha = 0.7) +
  labs(title = "Histogram of Unique Taxid Counts by Category", 
       x = "Unique Taxid Count", 
       y = "Frequency") +
  scale_fill_manual(values = category_colors) +
  theme_bw()+
  theme(axis.text.x = element_text(angle=75, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.73, 0.3835),
        legend.background = element_rect(fill = "white", color = "darkgray")) +
  facet_wrap(~ Category, scales = "free") +
  ylim(0,41) +
  geom_vline(data = annotation_data, 
             aes(xintercept = mean), 
             color = "red", linetype = "dashed", size = 1) +
  geom_text(data = annotation_data, 
            aes(x=x, y=41, label = label1), inherit.aes = FALSE, size = 4) +
  geom_text(data = annotation_data, 
            aes(x=x, y=39, label = label2), inherit.aes = FALSE, size = 4) +
  geom_text(data = annotation_data, 
            aes(x=x, y=37, label = label3), inherit.aes = FALSE, size = 4) +
  geom_text(data = annotation_data, 
            aes(x = x, y = 35, label = label4),
            color = "red", size = 4)



