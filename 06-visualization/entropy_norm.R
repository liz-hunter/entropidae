library(tidyverse)
library(viridis)
library(ggforce) 

# this is for manually processing this in R 
# recommend subsets for larger analyses, and/or using the python instead

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

# join the db_score and cv_group info
entropy_counts <- entropy_df %>%
  left_join(cv_long, by = c("taxid", "database", "name"))

# create title column & order by CV
entropy_counts <- entropy_counts %>%
  mutate(label = paste0(name, " (", taxid, ")")) %>%
  mutate(label = fct_reorder(label, -cv))

# format CoV text for labels
cv_labels <- entropy_counts %>%
  select(label, cv) %>%
  distinct() %>%
  mutate(cv_label = paste0("Coefficient of Variation: ", round(cv, 3))) %>%
  mutate(x = 10, y = 1500000)

# normalize - your choice
entropy_counts_norm <- entropy_counts %>%
  group_by(taxid, database) %>%
  mutate(
    total_count = sum(count, na.rm = TRUE),
    #count_fraction = count / total_count,
    count_percent = (count / total_count) * 100 # pick one 
  ) %>%
  ungroup()

# make a faceted line plot
ggplot(entropy_counts_norm) +
  geom_area(
    aes(x = entropy_bin, y = count_fraction, group = database, fill = db_score),
    alpha = 0.2, position = "identity") +
  geom_line(
    aes(x = entropy_bin, y = count_fraction, group = database, color = db_score),
    size = 1) +
  geom_text(
    data = cv_labels,
    aes(x = x, y = y, label = cv_label),
    vjust = 1.2, hjust = 1.1,
    inherit.aes = FALSE,
    size = 3) +
  scale_color_viridis_c(option = "plasma", name = "% Database Composition") +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  coord_cartesian(xlim = c(NA, 10), ylim = c(0, NA)) +
  facet_wrap(~ label) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 9)
  ) +
  labs(
    title = "Entropy in Classification by Taxon",
    x = "Entropy Score",
    y = "Count"
  )

# make a faceted normalized log scaled line plot
ggplot(entropy_counts_norm) +
  geom_ribbon(
    aes(x = entropy_bin, ymin = 1e-8, ymax = count_fraction, group = database, fill = db_score),
    alpha = 0.1
  ) +
  geom_line(
    aes(x = entropy_bin, y = count_fraction, group = database, color = db_score),
    size = 1
  ) +
  scale_y_log10(
    labels = scales::label_scientific(digits = 2),
    breaks = scales::log_breaks(n = 6)
  ) +
  scale_color_viridis_c(option = "plasma", name = "% Database Composition") +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  facet_wrap(~ label) +
  coord_cartesian(xlim = c(NA, 20)) +  # leave y scaling to scale_y_log10
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 9)
  ) +
  labs(
    title = "Normalized Entropy in Classification by Taxon (Log Scale)",
    x = "Entropy Score",
    y = "Fraction of Reads (log10)"
  )

# make a normalized proportional chart
ggplot(entropy_counts_norm) +
  geom_area(
    aes(x = entropy_bin, y = count_percent, group = database, fill = db_score),
    alpha = 0.2, position = "identity"
  ) +
  geom_line(
    aes(x = entropy_bin, y = count_percent, group = database, color = db_score),
    size = 1
  ) +
  scale_color_viridis_c(option = "plasma", name = "% Database Composition") +
  scale_fill_viridis_c(option = "plasma", guide = "none") +
  coord_cartesian(xlim = c(NA, 10), ylim = c(0, 30)) +  # fixed percent range
  facet_wrap(~ label) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    strip.text = element_text(size = 9)
  ) +
  labs(
    title = "Normalized Entropy in Classification by Taxon",
    x = "Entropy Score",
    y = "Percent of Reads"
  )
