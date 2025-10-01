library(tidyverse)
library(vegan)
library(readr)
library(forcats)

# set variables & paths
report_head = c("perc_comp", "incl_min_count", "excl_min_count", "level", "taxid", "name")
true_taxids = read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/all_true_taxids.txt", col_names="taxid")
direc = "/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/"
myfiles = list.files(path="/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect", pattern="*_inspect.txt", full.names=TRUE)

# read in the full kraken inspect reports for all databases
inspect = lapply(myfiles, function(x) {
  dat = read_tsv(x, col_names = report_head)
  dat$db = x  # add a column with the db ID
  dat$db = gsub(direc, "", dat$db) # cleanup
  dat$db = gsub("_inspect.txt", "", dat$db)
  return(dat)
})

# optionally extract inspection data for taxids represented in test data ONLY 
inspect_min = lapply(inspect, function(df) {
  df = inner_join(df, true_taxids, by = "taxid")
  df = df %>% select(taxid, perc_comp, name, db)
  return(df)
})

# make database ID the column header
taxid_stats = lapply(inspect_min, function(df) { # use original 'inspect' object to show all or 'inspect_min' for sampled taxa only
  #df <- filter(df, level == "S") # if you use the 'inspect' object, you need to filter by desired taxonomic level
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

# calculate boundaries for the coefficient of variation 
# may want to calculate this for the full dataset even if you only visualize taxa that are represented
Q1 <- summary(inspect_summary$cv)[['1st Qu.']]
Q3 <- summary(inspect_summary$cv)[['3rd Qu.']]
min_cv <- min(inspect_summary$cv, na.rm = TRUE)
max_cv <- max(inspect_summary$cv, na.rm = TRUE)

# function to grab extremes for any summary statistic
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

# grab the 10 top and bottom results for the value CV
top_cv <- get_extremes(inspect_summary, "cv", n=10, which="high") #higher the cv, the more variation
bottom_cv <- get_extremes(inspect_summary, "cv", n=10, which="low") #lower the cv, the lower the variation
cv_df <- rbind(bottom_cv, top_cv)

# optionally, grab the 10 top and bottom results for the range 
# top_range <- get_extremes(inspect_summary, "value_range", n=10, which="high") #higher the cv, the more variation
# bottom_range <- get_extremes(inspect_summary, "value_range", n=10, which="low") #lower the cv, the lower the variation
# range_df <- rbind(bottom_range, top_range)

# reformat for plotting
cv_plot <- inspect_summary %>%
  mutate(
    cv_label = ifelse(is.na(cv), paste0(name, " *"), name),
    cv_stat = case_when(
      is.na(cv) ~ "NA",
      cv <= Q1 ~ "Lower Quartile (Q1)",
      cv >= Q3 ~ "Upper Quartile (Q3)",
      TRUE ~ "Middle Quartile (Q2)"
    )
  ) %>%
  arrange(is.na(cv), desc(cv)) %>%
  mutate(cv_label = factor(cv_label, levels = rev(cv_label)))

# coefficient of variation bar plot
ggplot(cv_plot, aes(x = cv_label, y = cv)) +
  geom_col(aes(fill = cv_stat)) +  # map fill to cv_stat
  geom_text(aes(label = ifelse(is.na(cv), "NA", round(cv, 2)),
                color = cv_stat),
            hjust = -0.1, size = 2.5) +
  scale_fill_manual(values = c(
    "Upper Quartile (Q3)" = "#E5B611",
    "Lower Quartile (Q1)" = "#007CBA",
    "Middle Quartile (Q2)" = "#ADADAD",
    "NA" = "#2E2925"
  )) +
  scale_color_manual(values = c(
    "Upper Quartile (Q3)" = "black",
    "Lower Quartile (Q1)" = "black",
    "Middle Quartile (Q2)" = "black",
    "NA" = "#2E2925"
  )) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  coord_flip() +
  labs(title = "Coefficient of Variation Across Databases for a Given Taxa", 
       x = "", y = "Coefficient of Variation (CV)",
       color = "CV Category") +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5), 
    panel.background = element_rect(fill = "transparent", color = NA),
  )

# plot a heatmap of % database composition for each taxa
perc_plot <- merged_df %>%
  select(taxid, name, starts_with("virid_db")) %>%
  pivot_longer(cols = starts_with("virid_db"),
               names_to = "database", values_to = "percentage")

# add the CV values to the df for ordering of taxa
perc_merg <- inner_join(perc_plot, inspect_summary, by = c("taxid", "name")) %>%
  select(taxid, name, database, percentage, cv)

# fix the names
perc_merg$database <- str_replace(perc_merg$database, "virid_db", "DB-")

# fix the order of the databases 1-10
perc_merg$database <- factor(perc_merg$database,
                             levels = perc_merg$database %>%
                               unique() %>%
                               str_extract("\\d+") %>%           # Extract numeric part
                               as.integer() %>%                  # Convert to number
                               order() %>%
                               {perc_merg$database[.]}          # Reorder original levels
)

# organize the taxons by CV instead of alphabetically
cv_order <- perc_merg %>%
  distinct(taxid, name, cv) %>%
  arrange(desc(cv)) %>%
  pull(name)

perc_merg$name <- factor(perc_merg$name, levels = rev(cv_order))  # rev for highest at top

# heatmap of % database composition
ggplot(perc_merg, aes(x = database, y = name, fill = percentage)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#222C67", "#4F96C4", "#F3F9FC"),
    na.value = "grey90"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 70, hjust=1),
        plot.title = element_text(hjust = 0.5)) +
  labs(title = "Taxon Abundance in Bootstrapped Databases",
       x = "", y = "",
       fill = "% of Database")
