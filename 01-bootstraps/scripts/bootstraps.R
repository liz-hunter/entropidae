#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- paste(
  "Usage:",
  "  Rscript scripts/bootstraps.R --in <cleaned.tsv> --n_boot <int> --label <string>",
  "                      [--seed <int>] [--outdir DIR]",
  "                      [--outputs composite,files,unique_files,accession_counts,taxid_counts]",
  "                      [--quiet]",
  "",
  "Required:",
  "  --in        TSV with columns: accession, name, taxid",
  "  --n_boot    Number of bootstrap samples (with replacement)",
  "  --label     String used for naming outputs (e.g., plants)",
  "",
  "Optional:",
  "  --seed      Integer seed for reproducible sampling",
  "  --outdir    Output directory (default: {label}_boot)",
  "  --outputs   Comma-separated outputs to produce (reduces from default).",
  "             Allowed: composite, files, unique_files, accession_counts, taxid_counts",
  "             Default (if omitted): ALL outputs",
  "  --quiet     Suppress progress messages (still errors on invalid input).",
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

# parse required flags
infile     <- get_flag_value("--in", args)
n_boot_raw <- get_flag_value("--n_boot", args)
label      <- get_flag_value("--label", args)

missing_required <- character(0)
if (is.null(infile)) missing_required <- c(missing_required, "--in")
if (is.null(n_boot_raw)) missing_required <- c(missing_required, "--n_boot")
if (is.null(label)) missing_required <- c(missing_required, "--label")

if (length(missing_required) > 0) {
  cat("ERROR: Missing required flag(s): ", paste(missing_required, collapse = ", "),
      "\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

# optional flags
seed_raw    <- get_flag_value("--seed", args)
outdir      <- get_flag_value("--outdir", args)
outputs_raw <- get_flag_value("--outputs", args)
quiet <- "--quiet" %in% args

# validate unknown flags early
known_flags <- c("--in", "--n_boot", "--label", "--seed", "--outdir", "--outputs", "--quiet", "-h", "--help")
is_flag <- str_starts(args, "-")
unknown <- setdiff(args[is_flag], known_flags)
if (length(unknown) > 0) {
  cat("ERROR: Unknown option(s): ", paste(unknown, collapse = ", "),
      "\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

# validate values
n_boot <- suppressWarnings(as.integer(n_boot_raw))
if (is.na(n_boot) || n_boot <= 0) {
  cat("ERROR: --n_boot must be a positive integer.\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

seed <- NULL
if (!is.null(seed_raw)) {
  seed <- suppressWarnings(as.integer(seed_raw))
  if (is.na(seed)) {
    cat("ERROR: --seed must be an integer.\n\n", usage, "\n", sep = "")
    quit(status = 2)
  }
}

if (str_trim(label) == "") {
  cat("ERROR: --label must be non-empty.\n\n", usage, "\n", sep = "")
  quit(status = 2)
}

if (!is.null(outdir) && str_trim(outdir) == "") outdir <- NULL

# source lib relative to this script's location (needed for allowed_outputs)
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  match <- grep(file_arg, cmd_args, value = TRUE)
  if (length(match) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", match[1]), winslash = "/", mustWork = FALSE)))
  }
  normalizePath(getwd(), winslash = "/", mustWork = FALSE)
}

script_dir <- get_script_dir()
stage_dir <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
lib_path <- file.path(stage_dir, "lib", "bootstraps_lib.R")

if (!file.exists(lib_path)) {
  cat("ERROR: Could not find library file: ", lib_path, "\n", sep = "")
  quit(status = 2)
}

source(lib_path)

# outputs: default = everything, unless user specifies --outputs
if (is.null(outputs_raw) || str_trim(outputs_raw) == "") {
  outputs <- allowed_outputs
} else {
  outputs <- str_split(outputs_raw, ",")[[1]]
  outputs <- str_trim(outputs)
  outputs <- outputs[outputs != ""]
}

# run
run_bootstraps(
  infile = infile,
  n_boot = n_boot,
  label = label,
  seed = seed,
  outdir = outdir,
  outputs = outputs,
  verbose = !quiet
)
