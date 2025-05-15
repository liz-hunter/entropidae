#Setup and define functions
library(tidyverse, verbose=FALSE)

#Function for sampling with replacement ONLY
grabSR <- function(input, n, t) {
  lapply(1:t, function(i) {
    sampled_indices <- sample(seq_len(nrow(input)), n, replace = TRUE)
    input[sampled_indices, , drop = FALSE]
  })
}

# Function for sampling with replacement AND doing some basic calculations
# Outputs sample_reduc (just taxids), sample_acc (just accessions), and freq_data
# freq_data gives counts for each taxid, how many were present originally and how many were sampled in the bootstrap experiment
process_samples <- function(data, n_replicates = 100) {
  library(dplyr)
  
  # Select relevant columns
  data <- data %>% select(taxid = `Organism Taxonomic ID`, `Assembly Accession`)
  data_reduc <- data %>% select(taxid)
  
  # Sample with replacement
  samples <- grabSR(data, n = nrow(data), t = n_replicates)
  
  # Extract sampled columns
  samples_reduc <- lapply(samples, `[[`, "taxid")
  samples_acc <- lapply(samples, `[[`, "Assembly Accession")
  
  # Frequency calculations
  original_freq <- as.data.frame(table(data_reduc))
  colnames(original_freq) <- c("Value", "Original_Count")
  
  sample_freq <- as.data.frame(table(unlist(samples_reduc)))
  colnames(sample_freq) <- c("Value", "Sampled_Count")
  
  # Merge frequency counts
  freq_data <- merge(original_freq, sample_freq, by = "Value", all = TRUE)
  freq_data[is.na(freq_data)] <- 0
  
  # Return results as a list
  return(list(
    samples_reduc = samples_reduc,
    samples_acc = samples_acc,
    freq_data = freq_data
  ))
}

# Function to grab unique taxid counts for each sample
get_histogram_data <- function(samples_reduc, category) {
  counts <- lapply(samples_reduc, function(x) length(unique(x)))
  counts_vector <- unlist(counts)
  histogram_data <- data.frame(UniqueCounts = counts_vector, Category = category)
  return(histogram_data)
}

virid <- read_tsv("accessions/viridiplantae_all_ncbi20250303.tsv", show_col_types = FALSE)
bac <- read_tsv("accessions/bacteria_all_ncbi20250325.tsv", show_col_types = FALSE)
archea <- read_tsv("accessions/archea_all_ncbi20250325.tsv", show_col_types = FALSE)
virus <- read_tsv("accessions/viruses_all_ncbi20250325.tsv", show_col_types = FALSE)
euk <- read_tsv("accessions/euk_ncbi20250327.tsv", show_col_types = FALSE)

archea_results <- process_samples(archea, n=100)
bacteria_results <- process_samples(bac, n=100)
virid_results <- process_samples(virid, n=100)
virus_results <- process_samples(virus, n=100)
euk_results <- process_samples(euk, n=100)

# Grab the first 10 columns of accessions of the first 10 samples
grab <- virid_results$samples_acc[1:10] %>% data.frame()
colnames(grab) <- paste0("sample", 1:10)

write_tsv(grab, "samples/virid_samples.tsv")

# Sample a single accession from the dataset (in case of removal/suppression)
acc <- sample(virid$`Assembly Accession`, size=1)
print(acc)


