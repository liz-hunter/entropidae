library(taxonomizr)
library(tidyverse)
library(data.table)
library(data.tree)

# Set working directory
sql <- "/Users/elizabeth.hunter/taxonomizr_20250305/accessionTaxa.sql"
local <- "/Users/elizabeth.hunter/Library/CloudStorage/OneDrive-FDA/Documents/current_projects/entropy/"
remote <- "/Volumes/handy_group/liz_working/entropy/"

# Read in the raw taxonomy node data
nodes <- fread("/Users/elizabeth.hunter/Library/CloudStorage/OneDrive-FDA/Documents/current_projects/entropy/entropidae/06-entropy/taxdump/nodes.dmp", sep = "|", header = FALSE, quote = "", stringsAsFactors = FALSE)

# Extract only the columns we care about: tax_id, parent_tax_id, rank
tax_nodes <- nodes[, .(tax_id = as.integer(V1),
                       parent_id = as.integer(V2),
                       rank = trimws(V3))]

# Build a parent lookup table: named vector where names = tax_id, values = parent_id
parent_map <- setNames(tax_nodes$parent_id, tax_nodes$tax_id)

# Function to trace a taxid's lineage up to the root
get_lineage <- function(taxid) {
  path <- c()
  while (!is.na(taxid) && taxid != 1) {
    path <- c(taxid, path)
    taxid <- parent_map[as.character(taxid)]
  }
  return(c(1, path))  # include root node (1)
}

# Grab a report + corresponding key
base <- "citrus-rl-hiseq1M"
db <- "virid2"

# Read in the report and rename
report <- fread(paste0(remote, "kraken_reports/", base, "_", db, ".txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
report_names <- c("kraken_abundance", "reads_inclusive", "reads_direct", 
                  "minimizer_count", "minimizers_distinct", "ID_level", "taxID", "ID")
colnames(report) <- report_names

# Read in the out file, sort, filter, and rename
reads <- fread(paste0(remote, "kraken_reports/", base, "_", db, ".out"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
read_names <- c("CoU", "SequenceID", "taxid", "read_length", "lca2kmer")
colnames(reads) <- read_names

reads <- reads %>%
  separate(taxid, into = c("taxa", "taxid"), sep = " \\(taxid ", fill = "right")
reads$taxid <- gsub(")", "", reads$taxid)

reads_cut <- reads %>% select(SequenceID, taxid)

# Read in key file and reformat
key <- fread(paste0(remote, "keys/", base, "_key_final.txt"), sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(key) <- c("SequenceID", "accession", "taxid")
key$SequenceID <- gsub("@", "", key$SequenceID)
key$SequenceID <- gsub("/1", "", key$SequenceID)

# Calculate the distance between 2 different classifications
get_taxonomic_distance <- function(taxid1, taxid2) {
  lineage1 <- get_lineage(as.integer(taxid1))
  lineage2 <- get_lineage(as.integer(taxid2))
  
  # Find the lowest common ancestor (LCA)
  common <- intersect(lineage1, lineage2)
  if (length(common) == 0) return(NA)  # no common ancestor, rare
  
  lca <- common[which.max(match(common, lineage1))]  # deepest shared node
  
  # Compute distances from each taxid to LCA
  d1 <- length(lineage1) - match(lca, lineage1)
  d2 <- length(lineage2) - match(lca, lineage2)
  
  return(d1 + d2)
}

reads_merged <- merge(reads_cut, key[, .(SequenceID, true_taxid = taxid)], by = "SequenceID", all.x = TRUE)
reads_merged$tax_dist <- mapply(get_taxonomic_distance, reads_merged$taxid, reads_merged$true_taxid)
