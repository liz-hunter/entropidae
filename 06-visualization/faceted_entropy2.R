library(tidyverse)

# Load variables
report_head = c("perc_comp", "incl_min_count", "excl_min_count", "level", "taxid", "name")
true_taxids = read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/all_true_taxids.txt", col_names="taxid")
direc = "/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/"
myfiles = list.files(path="/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect", pattern="*_inspect.txt", full.names=TRUE)
csv_files <- list.files("/hpc/scratch/Elizabeth.Hunter/entropy/analysis/entropy_filtered/subsamples/",
                        pattern = "*.tsv",
                        full.names = TRUE)
present_taxids <- c(13894, 171251, 261450, 38859, 51239,
                    22663, 3755, 44588, 4577, 65409)

## ------- LOAD FUNCTIONS

# LOAD + CLEAN INSPECT DATA
load_inspect <- function(files, direc, report_head) {
  lapply(files, function(x) {
    dat <- read_tsv(x, col_names = report_head)
    dat$db <- gsub("_inspect.txt", "", gsub(direc, "", x))
    dat
  })
}


# FILTER TO TAXIDS OF INTEREST
filter_taxids <- function(inspect_list, taxid_df) {
  lapply(inspect_list, function(df) {
    df %>%
      inner_join(taxid_df, by = "taxid") %>%
      select(taxid, perc_comp, name, db)
  })
}

# PIVOT TO WIDE + MERGE
pivot_to_wide <- function(df_list) {
  taxid_stats <- lapply(df_list, function(df) { 
    colname <- unique(na.omit(df$db))
    if (length(colname) != 1) {
      warning("Multiple DB names found in one dataframe: ", paste(colname, collapse = ", "))
      colname <- colname[1] # choose first to avoid errors
    }
    df %>% select(taxid, name, !!colname := perc_comp)
  })
  merged_df <- reduce(taxid_stats, full_join, by = c("taxid", "name"))
  merged_df[is.na(merged_df)] <- 0
  merged_df
}

pivot_to_wide <- function(df_list) {
  taxid_stats <- lapply(df_list, function(df) { 
    colname <- unique(na.omit(df$db))
    if (length(colname) != 1) {
      warning("Multiple DB names found in one dataframe: ", paste(colname, collapse = ", "))
      colname <- colname[1] # choose first to avoid errors
    }
    
    # Fix DB naming quirks
    if (colname == "dist0") {
      colname <- "virid_db10"
    }
    
    df %>% select(taxid, name, !!colname := perc_comp)
  })
  
  merged_df <- reduce(taxid_stats, full_join, by = c("taxid", "name"))
  merged_df[is.na(merged_df)] <- 0
  merged_df
}

# SUMMARISE TAXON STATS
summarise_taxa <- function(merged_df, pattern = "^virid_db") {
  virid_cols <- grep(pattern, names(merged_df), value = TRUE)
  
  merged_df %>%
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
}

# GRAB EXTREMES 
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


# LOAD ENTROPY DATA
load_entropy_data <- function(file_paths) {
  bind_rows(lapply(file_paths, function(file_path) {
    read_tsv(file_path, show_col_types = FALSE)
  })) %>%
    rename(taxid = true_taxid) %>%
    pivot_longer(
      cols = starts_with("dist"),
      names_to = "database",
      values_to = "entropy"
    ) %>%
    mutate(
      database = str_replace(database, "dist", "DB-")
    ) %>%
    mutate(
      database = str_replace(database, "DB-0", "DB-10")
    )
}


# JOIN ENTROPY WITH CV DATA
join_entropy_with_cv <- function(entropy_df, cv_df, taxa_order, filter_taxids = NULL) {
  
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
  
  df <- entropy_df %>%
    inner_join(cv_long, by = c("taxid", "database")) %>%
    mutate(
      stat_group = factor(stat_group, levels = c("high", "low")),
      name = factor(name, levels = taxa_order)
    )
  
  if (!is.null(filter_taxids)) {
    df <- df %>% filter(taxid %in% filter_taxids)
  }
  
  db_levels <- unique(df$database) # set DB order numerically
  db_levels <- db_levels[order(as.numeric(str_remove(db_levels, "DB-")))]
  df <- df %>%
    mutate(database = factor(database, levels = db_levels))
  
  if (!is.null(filter_taxids)) {
    df <- df %>% filter(taxid %in% filter_taxids)
  }
  
  df
}


# PLOT IN BATCHES
plot_batches <- function(df, taxa_order, batch_size = 5, out_dir = "plots") {
  dir.create(out_dir, showWarnings = FALSE)
  
  batches <- split(taxa_order, ceiling(seq_along(taxa_order) / batch_size))
  
  cv_colors <- colorRampPalette(c("#A3CAE1", "#222C67"))(100)
  
  for (i in seq_along(batches)) {
    this_taxa <- batches[[i]]
    df_batch <- df %>% filter(name %in% this_taxa)
    
    dropped <- df_batch %>%
      filter(!is.finite(entropy) | entropy < 0 | entropy > 20)
    
    if (nrow(dropped) > 0) {
      message("Batch ", i, ": ", nrow(dropped), 
              " rows dropped due to non-finite or out-of-range entropy values.")
      print(
        dropped %>%
          group_by(stat_group, name, database) %>%
          summarise(
            min_entropy = min(entropy, na.rm = TRUE),
            max_entropy = max(entropy, na.rm = TRUE),
            count = n(),
            .groups = "drop"
          )
      )
    }
    
    # Plot
    p <- ggplot(df_batch, aes(x = database, y = entropy, fill = db_score)) +
      geom_violin(scale = "width", adjust = 1, trim = FALSE, color = NA, alpha = 0.7) +
      scale_fill_gradientn(colors = cv_colors, name = "% DB Composition", limits = c(0, 0.5)) +
      scale_y_continuous(limits = c(0, 20), oob = scales::oob_keep) +
      facet_grid(stat_group ~ name, scales = "fixed") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(size = 9)
      ) +
      labs(
        title = "Entropy Distributions",
        x = "Database", y = "Entropy Score"
      )
    
    ggsave(file.path(out_dir, paste0("entropy_batch_", i, ".png")),
           plot = p, width = 12, height = 8)
  }
}

# PERMUTATION TEST 
perm_test_global <- function(df, n_perm = 10) {
  groups <- unique(df$stat_group)
  if(length(groups) != 2) stop(paste("perm_test requires exactly 2 groups, but found", length(groups)))
  
  observed_diff <- df %>%
    group_by(stat_group) %>%
    summarise(mean_entropy = mean(entropy, na.rm = TRUE)) %>%
    summarise(diff = diff(mean_entropy)) %>%
    pull(diff)
  
  perm_diffs <- replicate(n_perm, {
    shuffled_groups <- sample(df$stat_group)
    perm_mean_diff <- df %>%
      mutate(stat_group = shuffled_groups) %>%
      group_by(stat_group) %>%
      summarise(mean_entropy = mean(entropy, na.rm = TRUE)) %>%
      summarise(diff = diff(mean_entropy)) %>%
      pull(diff)
    perm_mean_diff
  })
  
  p_value <- mean(abs(perm_diffs) >= abs(observed_diff))
  
  tibble(observed_diff = observed_diff, p_value = p_value)
}

### ------ RUN IT 

# Load and clean inspection reports
inspect <- load_inspect(myfiles, direc, report_head)

# Filter to relevant taxids
inspect_min <- filter_taxids(inspect, true_taxids)

# Pivot + merge to wide format
merged_df <- pivot_to_wide(inspect_min)

# Summarize statistics
inspect_summary <- summarise_taxa(merged_df)

# Get top/bottom CV taxa
top_cv <- get_extremes(inspect_summary, "cv", n = 10, which = "high")
bottom_cv <- get_extremes(inspect_summary, "cv", n = 10, which = "low")
cv_df <- rbind(bottom_cv, top_cv)

# Create taxa order for plotting
low_order <- cv_df %>%
  filter(stat_group == "low", taxid %in% present_taxids) %>%
  arrange(cv) %>%
  pull(name)

high_order <- cv_df %>%
  filter(stat_group == "high", taxid %in% present_taxids) %>%
  arrange(cv) %>%
  pull(name)

taxa_order <- c(low_order, high_order)

# Load entropy data
entropy_df <- load_entropy_data(csv_files)

# Join entropy with CV + filter to present taxa
df_joined <- join_entropy_with_cv(entropy_df, cv_df, taxa_order,
                                  filter_taxids = present_taxids)

# Plot in batches
plot_batches(df_joined, taxa_order, batch_size = 5, out_dir = "plots_new")

# Permutation test
perm_results <- perm_test_global(df_joined)
print(perm_results)
