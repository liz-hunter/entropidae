#Setup and define functions
library(tidyverse, verbose=FALSE)

setwd("/Users/elizabeth.hunter/Library/CloudStorage/OneDrive-FDA/Documents/current_projects/entropy/accessions")

virid <- read_tsv("accessions/viridiplantae_all_ncbi20250303.tsv", show_col_types = FALSE)
bac <- read_tsv("accessions/bacteria_all_ncbi20250325.tsv", show_col_types = FALSE)
archea <- read_tsv("accessions/archea_all_ncbi20250325.tsv", show_col_types = FALSE)
virus <- read_tsv("accessions/viruses_all_ncbi20250325.tsv", show_col_types = FALSE)
euk <- read_tsv("accessions/euk_ncbi20250327.tsv", show_col_types = FALSE)

#Function for sampling with replacement
grabSR <- function(input, n, t) {
  lapply(1:t, function(i) {
    sampled_indices <- sample(seq_len(nrow(input)), n, replace = TRUE)
    input[sampled_indices, , drop = FALSE]
  })
}

# Function for doing some basic calculations on the sampled table
# Outputs sample_reduc (just taxids), sample_acc (just accessions), and freq_data
# freq_data gives counts for each taxid, how many were present originally and how many were sampled in the bootstrap experiment
process_samples <- function(data, sampled_data) {
  
  # Rename and select relevant columns
  data <- data %>%
    rename(taxid = `Organism Taxonomic ID`,
           accession = `Assembly Accession`) %>%
    select(taxid, accession)
  
  sampled_data <- sampled_data %>%
    rename(taxid = `Organism Taxonomic ID`,
           accession = `Assembly Accession`) %>%
    select(taxid, accession)
  
  # Frequency table for original data
  original_freq <- data %>%
    count(taxid, name = "Original_Count")
  
  # Frequency table for sampled data
  sampled_freq <- sampled_data %>%
    count(taxid, name = "Sampled_Count")
  
  # Merge frequencies, fill NAs with 0
  freq_data <- full_join(original_freq, sampled_freq, by = "taxid") %>%
    mutate(across(everything(), ~replace_na(.x, 0)))
  
  # Extract lists of taxid/accession from sampled data
  samples_taxid <- sampled_data$taxid
  samples_accession <- sampled_data$accession
  
  # Return results
  return(list(
    samples_taxid = samples_taxid,
    samples_accession = samples_accession,
    freq_data = freq_data
  ))
}

# need to fix this - pulled grabSR out of process_samples to improve downstream flexibility
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


