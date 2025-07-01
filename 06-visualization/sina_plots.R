library(tidyverse)
library(ggforce)

install.packages("viridis")

# Setdirectory
data_dir <- "$MYPATH/entropy/summaries_mixes/test/"  

# Grab a list of TSVs
files <- list.files(data_dir, pattern = "1M_taxdist_per_read.tsv", full.names = TRUE)

# Helper function to read and reshape one file
read_and_reshape <- function(file) {
  df <- read_tsv(file, show_col_types = FALSE)
  
  sample_id <- tools::file_path_sans_ext(basename(file))
  
  df_long <- df %>%
    select(starts_with("dist")) %>%
    pivot_longer(
      cols = everything(),
      names_to = "database",
      values_to = "entropy"
    ) %>%
    mutate(sample_id = sample_id)
  
  return(df_long)
}

# Apply to all files and combine
all_data <- map_dfr(files, read_and_reshape)

# Enforce law!! and order!
all_data <- all_data %>%
  mutate(
    sample_id = factor(sample_id),
    database = factor(database, levels = paste0("dist", 1:6))
  )

# Sina plot with violin outline
ggplot(all_data, aes(x=sample_id, y=entropy, color=database)) +
  scale_y_log10() +
  geom_violin() +
  geom_sina() + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("Tax Entropy Across Databases") +
  xlab("SampleID")
