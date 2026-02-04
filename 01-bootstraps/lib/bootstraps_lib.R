suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
})

allowed_outputs <- c("composite", "files", "unique_files", "accession_counts", "taxid_counts")

sanitize_label <- function(label) {
  label %>%
    str_trim() %>%
    str_replace_all("\\s+", "_") %>%
    str_replace_all("[^A-Za-z0-9_-]", "_")
}

validate_outputs <- function(outputs) {
  bad <- setdiff(outputs, allowed_outputs)
  if (length(bad) > 0) {
    stop(
      "Unknown output(s): ", paste(bad, collapse = ", "),
      "\nAllowed: ", paste(allowed_outputs, collapse = ", "),
      call. = FALSE
    )
  }
  outputs
}

read_cleaned_input <- function(infile) {
  if (!file.exists(infile)) stop("Input file does not exist: ", infile, call. = FALSE)

  df <- read_tsv(infile, show_col_types = FALSE)
  needed <- c("accession", "name", "taxid")
  missing_cols <- setdiff(needed, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required column(s) in cleaned input: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  lookup <- df %>%
    select(accession, name, taxid) %>%
    distinct(accession, .keep_all = TRUE)

  n <- nrow(lookup)
  if (n == 0) stop("No accessions found in input.", call. = FALSE)

  lookup
}

make_bootstrap_matrices <- function(lookup, n_boot, seed = NULL) {
  n_boot <- suppressWarnings(as.integer(n_boot))
  if (is.na(n_boot) || n_boot <= 0) stop("n_boot must be a positive integer.", call. = FALSE)

  accessions <- lookup$accession
  taxids     <- lookup$taxid

  n <- length(accessions)

  if (!is.null(seed)) {
    seed <- suppressWarnings(as.integer(seed))
    if (is.na(seed)) stop("seed must be an integer.", call. = FALSE)
    set.seed(seed)
  }

  # sample indices for all bootstraps at once
  idx <- sample.int(n, size = n * n_boot, replace = TRUE)

  mat_acc <- matrix(accessions[idx], nrow = n, ncol = n_boot)
  mat_tax <- matrix(taxids[idx],     nrow = n, ncol = n_boot)

  list(mat_acc = mat_acc, mat_tax = mat_tax, idx = idx)
}

set_bootstrap_colnames <- function(mat, label_safe) {
  n_boot <- ncol(mat)
  colnames(mat) <- paste0(label_safe, "_db", seq_len(n_boot))
  mat
}

ensure_outdir <- function(outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  normalizePath(outdir, winslash = "/", mustWork = FALSE)
}

write_composites <- function(mat_acc, mat_tax, outdir, label_safe) {
  composite_acc_path <- file.path(outdir, paste0("composite_", label_safe, "_samples.tsv"))
  composite_tax_path <- file.path(outdir, paste0("composite_", label_safe, "_taxids.tsv"))

  write_tsv(as.data.frame(mat_acc), composite_acc_path)
  write_tsv(as.data.frame(mat_tax), composite_tax_path)

  list(accessions = composite_acc_path, taxids = composite_tax_path)
}

write_bootstrap_files <- function(mat_acc, outdir, label_safe) {
  n_boot <- ncol(mat_acc)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {
    fpath <- file.path(outdir, paste0(label_safe, "_db", i, ".txt"))
    writeLines(mat_acc[, i], con = fpath, sep = "\n")
    paths[i] <- fpath
  }

  paths
}

write_unique_bootstrap_files <- function(mat_acc, outdir, label_safe) {
  n_boot <- ncol(mat_acc)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {
    uniq <- unique(mat_acc[, i])
    fpath <- file.path(outdir, paste0(label_safe, "_db", i, "_unique.txt"))
    writeLines(uniq, con = fpath, sep = "\n")
    paths[i] <- fpath
  }

  paths
}

write_accession_count_tables <- function(lookup, idx, n_boot, outdir, label_safe) {
  accessions <- lookup$accession
  taxids     <- lookup$taxid
  names_vec  <- lookup$name
  n <- length(accessions)

  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {
    inds <- idx[((i - 1) * n + 1):(i * n)]
    t <- table(inds)
    ids <- as.integer(names(t))

    out <- tibble(
      accession = accessions[ids],
      name = names_vec[ids],
      taxid = taxids[ids],
      count = as.integer(t)
    ) %>%
      arrange(desc(count), accession)

    out_path <- file.path(outdir, paste0(label_safe, "_db", i, "_accession_counts.tsv"))
    write_tsv(out, out_path)
    paths[i] <- out_path
  }

  paths
}

write_taxid_count_tables <- function(mat_tax, outdir, label_safe) {
  n_boot <- ncol(mat_tax)
  paths <- character(n_boot)

  for (i in seq_len(n_boot)) {
    out <- tibble(taxid = mat_tax[, i]) %>%
      count(taxid, name = "count") %>%
      arrange(desc(count), taxid)

    out_path <- file.path(outdir, paste0(label_safe, "_db", i, "_taxid_counts.tsv"))
    write_tsv(out, out_path)
    paths[i] <- out_path
  }

  paths
}

run_bootstraps <- function(infile,
                           n_boot,
                           label,
                           seed = NULL,
                           outdir = NULL,
                           outputs = c("composite", "files"),
                           verbose = TRUE) {

  if (is.null(label) || str_trim(label) == "") stop("label must be non-empty.", call. = FALSE)

  label_safe <- sanitize_label(label)
  outputs <- validate_outputs(outputs)

  if (is.null(outdir) || str_trim(outdir) == "") outdir <- paste0(label_safe, "_boot")
  outdir_norm <- ensure_outdir(outdir)

  lookup <- read_cleaned_input(infile)
  mats <- make_bootstrap_matrices(lookup, n_boot = n_boot, seed = seed)

  mat_acc <- set_bootstrap_colnames(mats$mat_acc, label_safe)
  mat_tax <- set_bootstrap_colnames(mats$mat_tax, label_safe)

  results <- list(
    label = label,
    label_safe = label_safe,
    outdir = outdir_norm,
    n = nrow(lookup),
    n_boot = as.integer(n_boot),
    outputs = outputs
  )

  if ("composite" %in% outputs) {
    paths <- write_composites(mat_acc, mat_tax, outdir_norm, label_safe)
    results$composite_accessions <- paths$accessions
    results$composite_taxids <- paths$taxids
    if (verbose) {
      cat("Wrote composite accessions: ", paths$accessions, " (", nrow(lookup), " rows x ", ncol(mat_acc), " cols)\n", sep = "")
      cat("Wrote composite taxids:     ", paths$taxids,     " (", nrow(lookup), " rows x ", ncol(mat_tax), " cols)\n", sep = "")
    }
  }

  if ("files" %in% outputs) {
    results$bootstrap_files <- write_bootstrap_files(mat_acc, outdir_norm, label_safe)
    if (verbose) cat("Wrote ", length(results$bootstrap_files), " bootstrap accession list files to: ", outdir_norm, "\n", sep = "")
  }

  if ("unique_files" %in% outputs) {
    results$unique_bootstrap_files <- write_unique_bootstrap_files(mat_acc, outdir_norm, label_safe)
    if (verbose) cat("Wrote ", length(results$unique_bootstrap_files), " unique-accession list files to: ", outdir_norm, "\n", sep = "")
  }

  if ("accession_counts" %in% outputs) {
    results$accession_count_tables <- write_accession_count_tables(
      lookup = lookup, idx = mats$idx, n_boot = as.integer(n_boot),
      outdir = outdir_norm, label_safe = label_safe
    )
    if (verbose) cat("Wrote ", length(results$accession_count_tables), " per-bootstrap accession count tables to: ", outdir_norm, "\n", sep = "")
  }

  if ("taxid_counts" %in% outputs) {
    results$taxid_count_tables <- write_taxid_count_tables(mat_tax, outdir_norm, label_safe)
    if (verbose) cat("Wrote ", length(results$taxid_count_tables), " per-bootstrap taxid count tables to: ", outdir_norm, "\n", sep = "")
  }

  if (verbose) cat("Output directory: ", outdir_norm, "\n", sep = "")

  invisible(results)
}
