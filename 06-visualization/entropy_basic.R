library(tidyverse)
library(viridis)

# grab the count distribution for samples containing these taxids 
myfiles = list.files(path="/hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_counts", pattern="*_entropy_counts.tsv", full.names=TRUE)

# read in and bind all files
output = lapply(myfiles, function(x) {
  dat = read_tsv(x, col_names = TRUE)
  return(dat)
})

entropy_df <- bind_rows(output)

# alternatively, for a single file:
#entropy_df <- read_tsv("analysis/entropy/lod-extract-12-01-miseq1M_taxdist_per_read.tsv")
#entropy_df <- entropy_df %>% filter(true_taxid == 38859, preserve=TRUE)

# format as long data
entropy_long <- entropy_df %>%
  pivot_longer(cols = starts_with("dist"),
               names_to = "database",
               values_to = "entropy") %>%
  mutate(
    database = gsub("dist", "virid", database),
    entropy_bin = round(entropy)
  ) %>%
  rename(taxid = true_taxid) 

# get counts
entropy_counts <- entropy_long %>%
  group_by(taxid, database, entropy_bin) %>%
  summarise(count = n(), .groups = "drop")

# drop NAs (unclassified)
entropy_counts_clean <- entropy_counts %>%
  filter(!is.na(entropy_bin), !is.na(count))

# Join the db_score and cv_group info
entropy_counts_with_score <- entropy_counts_clean %>%
  left_join(cv_long, by = c("taxid", "database"))

summary(entropy_counts_clean)
midpoint <- median(entropy_counts_with_score$db_score, na.rm = TRUE)

# setup plot names
my_name <- unique(entropy_counts_with_score$name)
my_taxid <- unique(entropy_counts_with_score$taxid)
plot_title <- paste0("Classifications for ", my_name, " (", my_taxid, ")")

# plot single plot
ggplot() +
  geom_area(
    data = entropy_counts_with_score,
    aes(x = entropy_bin, y = count, group = database, fill = db_score),
    alpha = 0.4, position = "identity"
  ) +
  geom_line(
    data = entropy_counts_with_score,
    aes(x = entropy_bin, y = count, group = database, color = db_score),
    size = 1
  ) +
  scale_color_viridis_c(option = "plasma", name = "% Representation") +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(NA, 10), ylim = c(0, NA)) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(title = plot_title, x = "Entropy Bin", y = "Count")

