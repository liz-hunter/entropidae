suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
})

sanitize_label <- function(label) {
  label %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_-]", "_")
}

# Diversity functions
shannon_H <- function(counts) {
  counts <- counts[counts > 0]
  p <- counts / sum(counts)
  -sum(p * log(p))
}

simpson_D <- function(counts) {
  counts <- counts[counts > 0]
  p <- counts / sum(counts)
  sum(p^2)  # dominance
}

gini_simpson <- function(counts) 1 - simpson_D(counts)

inverse_simpson <- function(counts) {
  D <- simpson_D(counts)
  if (D == 0) NA_real_ else 1 / D
}

pielou_evenness <- function(counts) {
  S <- sum(counts > 0)
  if (S <= 1) return(NA_real_)
  H <- shannon_H(counts)
  H / log(S)
}

read_taxid_composite <- function(taxids_path) {
  if (!file.exists(taxids_path)) stop("Taxids composite file not found: ", taxids_path, call. = FALSE)

  tax_mat <- read_tsv(taxids_path, show_col_types = FALSE)
  if (ncol(tax_mat) < 1) stop("Taxids composite file has no columns.", call. = FALSE)

  tax_mat %>%
    mutate(.row = row_number()) %>%
    pivot_longer(cols = - .row, names_to = "sample_id", values_to = "taxid") %>%
    select(sample_id, taxid) %>%
    filter(!is.na(taxid), str_trim(as.character(taxid)) != "") %>%
    count(sample_id, taxid, name = "count")
}

# Calculate and compile diversity metrics
compute_diversity_summary <- function(tax_counts) {
  tax_counts %>%
    group_by(sample_id) %>%
    summarise(
      n_draws = sum(count),
      n_unique_taxa = n_distinct(taxid),
      shannon = shannon_H(count),
      simpson = simpson_D(count),
      gini_simpson = gini_simpson(count),
      inverse_simpson = inverse_simpson(count),
      pielou_evenness = pielou_evenness(count),
      .groups = "drop"
    ) %>%
    arrange(sample_id)
}

# Write a nice summary in a new directory
write_diversity_summary <- function(summary_tbl, out_path) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  write_tsv(summary_tbl, out_path)
  out_path
}

# Plotting functions
metrics_to_long <- function(summary_tbl,
                            metrics = c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness")) {
  keep <- intersect(metrics, names(summary_tbl))
  if (length(keep) == 0) stop("No requested metrics found in summary table.", call. = FALSE)

  summary_tbl %>%
    select(sample_id, all_of(keep)) %>%
    pivot_longer(cols = -sample_id, names_to = "metric", values_to = "value")
}

plot_metric_distributions <- function(summary_tbl,
                                      metrics = c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness"),
                                      facet_scales = "free_y") {
  long <- metrics_to_long(summary_tbl, metrics = metrics)

  ggplot(long, aes(x = value)) +
    geom_histogram(bins = 30) +
    facet_wrap(~ metric, scales = facet_scales) +
    labs(x = NULL, y = "Count")
}

plot_metric_boxplots <- function(summary_tbl,
                                 metrics = c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness"),
                                 facet_scales = "free_y") {
  long <- metrics_to_long(summary_tbl, metrics = metrics)

  ggplot(long, aes(x = metric, y = value)) +
    geom_boxplot(outlier.alpha = 0.4) +
    facet_wrap(~ metric, scales = facet_scales) +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
}

plot_metrics_facet <- function(summary_tbl,
                               metrics = c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness")) {
  keep <- intersect(metrics, names(summary_tbl))
  if (length(keep) == 0) stop("No requested metrics found in summary table.", call. = FALSE)

  int_metrics <- c("n_draws", "n_unique_taxa")

  fmt_value <- function(metric, value) {
    if (is.na(value)) return(NA_character_)
    if (metric %in% int_metrics) {
      return(formatC(round(value), format = "d"))
    }
    # default: 2 decimals
    formatC(value, format = "f", digits = 2)
  }

  long <- summary_tbl %>%
    select(sample_id, all_of(keep)) %>%
    pivot_longer(cols = -sample_id, names_to = "metric", values_to = "value") %>%
    mutate(
      db_index = as.integer(stringr::str_match(sample_id, "_db(\\d+)$")[, 2])
    )

  if (all(is.na(long$db_index))) {
    long <- long %>%
      group_by(metric) %>%
      mutate(db_index = row_number()) %>%
      ungroup()
  }

  long <- long %>%
    mutate(
      db_index = factor(db_index, levels = sort(unique(db_index))),
      value_label = mapply(fmt_value, metric, value, USE.NAMES = FALSE)
    )

  ggplot(long, aes(x = db_index, y = value, fill = metric)) +
    geom_col(show.legend = FALSE) +
    geom_text(
      aes(label = value_label),
      vjust = -0.25,
      size = 3,
      na.rm = TRUE
    ) +
    facet_wrap(~ metric, scales = "free_y") +
    labs(x = "Database / bootstrap replicate", y = NULL) +
    expand_limits(y = 0) +
    theme(
      strip.text = element_text(face = "bold")
    )
}

save_plot <- function(p, path, width = 10, height = 6, dpi = 300) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggsave(filename = path, plot = p, width = width, height = height, dpi = dpi)
  path
}

# Wrapper
run_diversity_metrics <- function(label,
                                  bootdir = NULL,
                                  taxids_path = NULL,
                                  out_path = NULL,
                                  make_plots = FALSE,
                                  plots_outdir = NULL,
                                  plot_metrics = c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness"),
                                  verbose = TRUE) {
  if (is.null(label) || stringr::str_trim(label) == "") {
    stop("label must be non-empty.", call. = FALSE)
  }
  label_safe <- sanitize_label(label)

  # Defaults (still here for direct library use)
  if (is.null(bootdir) || stringr::str_trim(bootdir) == "") bootdir <- paste0(label_safe, "_boot")
  if (!dir.exists(bootdir)) stop("bootdir does not exist: ", bootdir, call. = FALSE)

  if (is.null(taxids_path) || stringr::str_trim(taxids_path) == "") {
    taxids_path <- file.path(bootdir, paste0("composite_", label_safe, "_taxids.tsv"))
  }

  # Default output base directory
  metrics_dir <- paste0(label_safe, "_metrics")

  # Default summary TSV goes into metrics_dir
  if (is.null(out_path) || stringr::str_trim(out_path) == "") {
    out_path <- file.path(metrics_dir, paste0("diversity_", label_safe, "_summary.tsv"))
  }

  tax_counts <- read_taxid_composite(taxids_path)
  summary_tbl <- compute_diversity_summary(tax_counts)
  write_diversity_summary(summary_tbl, out_path)

  if (verbose) {
    cat(
      "Wrote diversity summary (taxid-based): ",
      out_path, " (", nrow(summary_tbl), " samples)\n",
      sep = ""
    )
  }

  plot_paths <- list()

  if (isTRUE(make_plots)) {
    # Default plots directory: directly into metrics_dir
    if (is.null(plots_outdir) || stringr::str_trim(plots_outdir) == "") {
      plots_outdir <- metrics_dir
    }

    p_facet <- plot_metrics_facet(summary_tbl, metrics = plot_metrics)

    plot_paths$facet_by_database <- save_plot(
      p_facet,
      file.path(plots_outdir, paste0("diversity_", label_safe, "_metrics_by_database.png")),
      width = 12, height = 7
    )

    if (verbose) {
      cat("Wrote plots to: ", normalizePath(plots_outdir, winslash = "/", mustWork = FALSE), "\n", sep = "")
    }
  }

  invisible(list(
    label = label,
    label_safe = label_safe,
    bootdir = normalizePath(bootdir, winslash = "/", mustWork = FALSE),
    taxids_path = taxids_path,
    out_path = out_path,
    summary = summary_tbl,
    plots = plot_paths
  ))
}
