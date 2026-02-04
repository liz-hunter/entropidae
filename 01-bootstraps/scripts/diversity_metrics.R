#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript scripts/diversity_metrics.R --label <string> [--bootdir DIR] [--taxids FILE] [--out FILE]",
  "                                      [--plots] [--plots_outdir DIR] [--plot_metrics m1,m2,...] [--quiet]",
  "",
  "Defaults (if not provided):",
  "  --bootdir        {label}_boot",
  "  --taxids         {bootdir}/composite_{label}_taxids.tsv",
  "  --out            {label}_metrics/diversity_{label}_summary.tsv",
  "  --plots_outdir   {label}_metrics",
  "",
  "Plot options:",
  "  --plots          Write basic plots (facet hist + facet boxplot).",
  "  --plot_metrics   Comma-separated metrics to include in facet plots.",
  "                  Default: n_unique_taxa,shannon,gini_simpson,inverse_simpson,pielou_evenness",
  "",
  "Example:",
  "  Rscript scripts/diversity_metrics.R --label plants --plots",
  sep = "\n"
)

if (length(args) == 0 || any(args %in% c("-h", "--help"))) {
  cat(usage, "\n")
  quit(status = 0)
}

get_flag_value <- function(flag, args) {
  i <- which(args == flag)
  if (length(i) == 0) return(NULL)
  if (i[length(i)] == length(args)) stop("Flag ", flag, " requires a value.\n\n", usage, call. = FALSE)
  args[i[length(i)] + 1]
}

label <- get_flag_value("--label", args)
if (is.null(label) || str_trim(label) == "") {
  cat("ERROR: --label is required and must be non-empty.\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

label_safe <- label %>%
  str_trim() %>%
  str_replace_all("\\s+", "_") %>%
  str_replace_all("[^A-Za-z0-9_-]", "_")

bootdir <- get_flag_value("--bootdir", args)
taxids_path <- get_flag_value("--taxids", args)
out_path <- get_flag_value("--out", args)

make_plots <- "--plots" %in% args
plots_outdir <- get_flag_value("--plots_outdir", args)

plot_metrics_raw <- get_flag_value("--plot_metrics", args)
if (is.null(plot_metrics_raw) || str_trim(plot_metrics_raw) == "") {
  plot_metrics <- c("n_unique_taxa", "shannon", "gini_simpson", "inverse_simpson", "pielou_evenness")
} else {
  plot_metrics <- str_split(plot_metrics_raw, ",")[[1]]
  plot_metrics <- str_trim(plot_metrics)
  plot_metrics <- plot_metrics[plot_metrics != ""]
}

quiet <- "--quiet" %in% args

# validate unknown flags
known_flags <- c("--label", "--bootdir", "--taxids", "--out", "--plots", "--plots_outdir", "--plot_metrics", "--quiet", "-h", "--help")
is_flag <- str_starts(args, "-")
unknown <- setdiff(args[is_flag], known_flags)
if (length(unknown) > 0) {
  cat("ERROR: Unknown option(s): ", paste(unknown, collapse = ", "),
      "\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

# ---- restore defaults here (so behavior is predictable in wrapper) ----
if (is.null(bootdir) || str_trim(bootdir) == "") {
  bootdir <- paste0(label_safe, "_boot")
}

if (is.null(taxids_path) || str_trim(taxids_path) == "") {
  taxids_path <- file.path(bootdir, paste0("composite_", label_safe, "_taxids.tsv"))
}

if (is.null(out_path) || str_trim(out_path) == "") {
  out_path <- file.path(paste0(label_safe, "_metrics"), paste0("diversity_", label_safe, "_summary.tsv"))
}

# Warning (stdout) so this doesn't look "silent" if stderr is swallowed
if (!dir.exists(bootdir)) {
  cat("WARNING: Default bootdir not found: ", bootdir, "\n", sep = "")
  cat("         If your bootstraps directory is elsewhere, pass --bootdir DIR\n", sep = "")
  quit(status = 2)
}

if (!file.exists(taxids_path)) {
  cat("WARNING: Taxids composite file not found: ", taxids_path, "\n", sep = "")
  cat("         Pass --taxids FILE (or confirm bootstraps were generated for this label)\n", sep = "")
  quit(status = 2)
}

# source lib relative to this script
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  match <- grep("--file=", cmd_args, value = TRUE)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub("--file=", "", match[1]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
stage_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
lib_path <- file.path(stage_dir, "lib", "diversity_metrics_lib.R")

if (!file.exists(lib_path)) {
  cat("ERROR: Could not find library file: ", lib_path, "\n", sep = "")
  quit(status = 2)
}

source(lib_path)

run_diversity_metrics(
  label = label,
  bootdir = bootdir,
  taxids_path = taxids_path,
  out_path = out_path,
  make_plots = make_plots,
  plots_outdir = plots_outdir,
  plot_metrics = plot_metrics,
  verbose = !quiet
)
