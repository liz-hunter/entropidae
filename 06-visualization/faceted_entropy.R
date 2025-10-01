library(tidyverse)
library(viridis)
library(ggforce) 

# set variables & paths
report_head = c("perc_comp", "incl_min_count", "excl_min_count", "level", "taxid", "name")
true_taxids = read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/all_true_taxids.txt", col_names="taxid")
direc = "/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/"
myfiles = list.files(path="/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect", pattern="*_inspect.txt", full.names=TRUE)

# define extremes function
get_extremes <- function(df, sum_stat, n = 10, which = c("high", "low")) {
  ToB <- match.arg(which)
  
  df_filtered <- df %>%
    filter(!is.na(.data[[sum_stat]])) %>%
    arrange(.data[[sum_stat]])
  
  extreme_rows <- if (ToB == "low") {
    slice_head(df_filtered, n = n)
  } else {
    slice_tail(df_filtered, n = n)
  }
  
  extreme_rows %>%
    mutate(stat_group = ToB)
}

# read in the full kraken inspect reports for all databases
inspect = lapply(myfiles, function(x) {
  dat = read_tsv(x, col_names = report_head)
  dat$db = x  # add a column with the db ID
  dat$db = gsub(direc, "", dat$db) # cleanup
  dat$db = gsub("_inspect.txt", "", dat$db)
  return(dat)
})

# extract inspection data for taxids represented in the test data
inspect_min = lapply(inspect, function(df) {
  df = inner_join(df, true_taxids, by = "taxid")
  df = df %>% select(taxid, perc_comp, name, db)
  return(df)
})

# make database ID the column header
taxid_stats = lapply(inspect_min, function(df) { 
  new_colname <- unique(na.omit(df$db)) 
  df %>%
    select(taxid, name, !!new_colname := perc_comp)
})

# merge into a single df, and replace NAs with 0 (0% composition)
merged_df <- reduce(taxid_stats, full_join, by = c("taxid", "name"))
merged_df[is.na(merged_df)] <- 0

# add statistical summaries for each taxon (max, min, mean, range, std deviation, & coefficient of variation)
virid_cols <- grep("^virid_db", names(merged_df), value = TRUE)

inspect_summary <- merged_df %>%
  rowwise() %>%
  mutate(
    perc_values = list(c_across(all_of(virid_cols))),
    max_val = max(perc_values),
    min_val = min(perc_values),
    value_range = max_val - min_val,
    std_dev = round(sd(perc_values), 3),
    mean_val = round(mean(perc_values), 3),
    cv = ifelse(mean_val == 0, NA, round(std_dev / mean_val, 3))
  ) %>%
  ungroup() %>%
  select(-perc_values)

# grab the 10 top and bottom results for the value CV
top_cv <- get_extremes(inspect_summary, "cv", n=10, which="high") #higher the cv, the more variation
bottom_cv <- get_extremes(inspect_summary, "cv", n=10, which="low") #lower the cv, the lower the variation
cv_df <- rbind(bottom_cv, top_cv)

# grab the names of the taxids you want to plot
taxids_of_interest <- cv_df %>% select(taxid, name, cv, stat_group) # grab the taxids you want to graph

# make list of taxids represented in dataset (parsed from filenames)
present_taxids <- c(13894, 171251, 261450, 38859, 51239,
                    22663, 3755, 44588, 4577, 65409)

# Filter cv_df to only taxa present in your dataset
cv_present <- cv_df %>%
  filter(taxid %in% present_taxids)

# sort each group by cv ascending and get names in order
low_order <- cv_present %>%
  filter(stat_group == "low") %>%
  arrange(cv) %>%
  pull(name)

high_order <- cv_present %>%
  filter(stat_group == "high") %>%
  arrange(cv) %>%
  pull(name)

# Combine low then high
taxa_order <- c(low_order, high_order)

# pivot the coefficient of variation info to plot 
cv_long <- cv_df %>%
  pivot_longer(
    cols = starts_with("virid_db"),
    names_to = "database",
    values_to = "db_score"
  ) %>%
  mutate(
    database = str_replace(database, "^virid_db", "DB-")
  ) %>%
  select(taxid, name, stat_group, cv, database, db_score)

# grab pre-formatted and filtered files
csv_files <- list.files("/hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_filtered/subsamples/",
                        pattern = "*.tsv",
                        full.names = TRUE)

# coerce into a single df, long format 
df_all <- bind_rows(lapply(csv_files, function(file_path) {
  read_tsv(file_path, show_col_types = FALSE)
}))

df_all_long <- df_all %>%
  rename(taxid = true_taxid) %>%
  pivot_longer(
    cols = starts_with("dist"),
    names_to = "database",
    values_to = "entropy"
  ) %>%
  mutate(
    database = str_replace(database, "dist", "DB-"))  

# join with the other database summary data
df_joined <- df_all_long %>%
  inner_join(cv_long %>% select(taxid, database, db_score, cv, stat_group, name), 
             by = c("taxid", "database"))

df_joined <- df_joined %>%
  mutate(
    stat_group = factor(stat_group, levels = c("high", "low")),
    name = factor(name, levels = taxa_order)
  )

facet_settings <- facet_grid(stat_group ~ name, scales = "fixed")

# set a base theme for consistent formatting
base_theme <- theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(size = 9)
  )

ggplot(df_joined, aes(x = database, y = entropy, fill = db_score)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, color = NA, alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "% DB Composition", limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(0, 10), oob = scales::oob_keep) +
  facet_grid(stat_group ~ name) +
  base_theme +
  labs(title = "Violin Plot of Entropy Distributions", x = "Database", y = "Entropy Score")

# og violin plot
ggplot(df_joined, aes(x = database, y = entropy, fill = db_score)) +
  geom_violin(scale = "width", adjust = 1, trim = FALSE, color = NA, alpha = 0.7) +
  scale_fill_viridis_c(option = "plasma", name = "% DB Composition", limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(0, 10), oob = scales::oob_keep) +
  facet_settings +
  base_theme +
  labs(title = "Violin Plot of Entropy Distributions", x = "Database", y = "Entropy Score")

# box plot
ggplot(df_all, aes(x = database, y = entropy, fill = db_score)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.6) +
  scale_fill_viridis_c(option = "plasma", name = "% DB Composition", limits = c(0, 0.5)) +
  scale_y_continuous(limits = c(0, 10), oob = scales::oob_keep) +
  facet_settings +
  base_theme +
  labs(title = "Boxplot of Entropy Distributions", x = "Database", y = "Entropy Score")
