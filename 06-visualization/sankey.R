library(tidyverse)
library(vegan)
library(readr)
library(forcats)
library(taxonomizr)
library(plotly)

# set variables & paths
report_head = c("perc_comp", "incl_min_count", "excl_min_count", "level", "taxid", "name")
true_taxids = read_tsv("/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/all_true_taxids.txt", col_names="taxid")
direc = "/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect/"
myfiles = list.files(path="/flash/storage/scratch/Elizabeth.Hunter/entropy/analysis/inspect", pattern="*_inspect.txt", full.names=TRUE)
sqlpath="/projects/handy_grp/taxonomizr20250305/accessionTaxa.sql"

# read in the kraken inspect reports
inspect = lapply(myfiles, function(x) {
  dat = read_tsv(x, col_names = report_head)
  dat$db = x  # add a column with the db ID
  dat$db = gsub(direc, "", dat$db) # cleanup
  dat$db = gsub("_inspect.txt", "", dat$db)
  return(dat)
})

# make database ID the column header
taxid_stats = lapply(inspect, function(df) { # use original 'inspect' object to show all or 'inspect_min' for sampled taxa only
  new_colname <- unique(na.omit(df$db)) 
  df %>%
    select(taxid, name, level, !!new_colname := perc_comp)
})

# merge into a single df, and replace NAs with 0 (0% composition)
merged_df <- reduce(taxid_stats, full_join, by = c("taxid", "name", "level"))
merged_df[is.na(merged_df)] <- 0

# add statistical summaries for each taxon (max, min, mean, range, std deviation, & coefficient of variation)
virid_cols <- grep("^virid_db", names(merged_df), value = TRUE)

inspect_summary <- merged_df %>%
  rowwise() %>%
  mutate(
    perc_values = list(c_across(all_of(virid_cols))),
    std_dev = round(sd(perc_values), 3),
    mean_val = round(mean(perc_values), 3),
    cv = ifelse(mean_val == 0, NA, round(std_dev / mean_val, 3))
  ) %>%
  ungroup() %>%
  select(-any_of(virid_cols)) %>%
  select(-perc_values)

desired_taxa <- c("D", "P", "C", "O", "F")
desired_taxa_full <- c("domain", "phylum", "class", "order", "family")

# grab the desired major levels 
major_only <- inspect_summary %>% 
  filter(!str_detect(level, "[0-9]$")) %>%
  filter(level %in% desired_taxa) %>% 
  filter(!is.na(cv))

# get the full taxonomy for everything the remains
taxid_vec <- major_only$taxid

taxonomy_matrix <- getTaxonomy(taxid_vec, sqlFile = sqlpath)

lineage_df <- taxonomy_matrix %>%
  as.data.frame() %>%
  rownames_to_column("taxid") %>%
  mutate(taxid = gsub("[^0-9]", "", taxid),  # Remove ALL non-digit characters
         taxid = as.numeric(taxid)) %>%
  select(all_of(desired_taxa_full), taxid)

# merge in the abundance data
all_info <- major_only %>%
  select(taxid, mean_val, cv) %>%
  left_join(lineage_df, by = "taxid")

# gather up the info 
df_tax <- all_info %>%
  mutate(across(c(desired_taxa_full), as.character))

# make the link df
links <- df_tax %>%
  rowwise() %>%
  mutate(
    taxa = list(na.omit(c_across(all_of(desired_taxa_full)))),
    valid_idx = list(seq_along(taxa)),
    n_valid = length(taxa),
    source = if(n_valid >= 2) taxa[n_valid - 1] else NA_character_,
    target = if(n_valid >= 2) taxa[n_valid] else NA_character_
  ) %>%
  ungroup() %>%
  filter(!is.na(source) & !is.na(target)) %>%
  select(source, target, mean_val, cv) %>%
  rename(value = mean_val)

nodes <- data.frame(name = unique(c(links$source, links$target)))

links <- links %>%
  mutate(
    source_id = match(source, nodes$name) - 1,
    target_id = match(target, nodes$name) - 1
  )

# Create a color scale function for your cv values, e.g. blue (low) to red (high)
cv_colors <- colorRampPalette(c("#A3CAE1", "#E5B611"))(100)

# Map cv values to 1-100 range for indexing colors
cv_scaled <- as.integer(cut(links$cv, breaks = 100, labels = FALSE))
link_colors <- cv_colors[cv_scaled]

plot_ly(
  type = "sankey",
  arrangement = "snap",
  node = list(
    #label = nodes$name,
    pad = 50,
    thickness = 20,
    line = list(color = "black", width = 0.5),
    color = "#007CBA"  # single color for all nodes
  ),
  link = list(
    source = links$source_id,
    target = links$target_id,
    value = links$value,
    color = link_colors  # keep the colored links by cv
  )
)

plot_ly(
  type = "sankey",
  arrangement = "snap",
  node = list(
    pad = 20,
    thickness = 30,
    line = list(color = "black", width = 0.5),
    color = "gray"
  ),
  link = list(
    source = links$source_id,
    target = links$target_id,
    value = links$value,
    color = link_colors
  )
) %>%
  layout(margin = list(r = 150))

# add a custom color gradient legend

# Convert to Plotly colorscale format
plotly_colorscale <- lapply(seq_along(cv_colors), function(i) {
  list((i - 1) / (length(cv_colors) - 1), cv_colors[i])
})

# Dummy scatter trace to create just the colorbar
plot_ly(
  type = "scatter",
  x = c(0), y = c(0),  # dummy point
  mode = "markers",
  marker = list(
    size = 0.0001, # make the point invisible
    colorscale = plotly_colorscale,
    cmin = 0,
    cmax = 3,
    colorbar = list(
      title = "CV",
      tickvals = seq(0, 3, by = 0.5),
      ticks = "outside",
      len = 1 # full height
    ),
    color = c(0) # just one value to attach the scale to
  ),
  showlegend = FALSE,
  hoverinfo = "none"
)

