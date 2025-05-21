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
