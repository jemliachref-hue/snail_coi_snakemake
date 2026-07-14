## Rule: download
## Sections 1-6 of original script

library(rentrez)
library(Biostrings)
library(dplyr)
library(stringr)
library(tibble)
library(readr)

# --- snakemake I/O -----------------------------------------------------
out_fasta_raw <- snakemake@output[["fasta_raw"]]
out_metadata  <- snakemake@output[["metadata"]]

families    <- snakemake@params[["families"]]
gene        <- snakemake@params[["gene"]]
batch_size  <- snakemake@params[["batch_size"]]
max_records <- snakemake@params[["max_records"]]
sleep_time  <- snakemake@params[["sleep_time"]]

# --- Build NCBI queries -------------------------------------------------
european_countries <- c(
  "Albania", "Andorra", "Austria", "Belarus", "Belgium",
  "Bosnia and Herzegovina", "Bulgaria", "Croatia", "Cyprus",
  "Czech Republic", "Denmark", "Estonia", "Finland", "France",
  "Germany", "Greece", "Hungary", "Iceland", "Ireland", "Italy",
  "Kosovo", "Latvia", "Liechtenstein", "Lithuania", "Luxembourg",
  "Malta", "Moldova", "Monaco", "Montenegro", "Netherlands",
  "North Macedonia", "Norway", "Poland", "Portugal", "Romania",
  "Russia", "San Marino", "Serbia", "Slovakia", "Slovenia",
  "Spain", "Sweden", "Switzerland", "Turkey", "Ukraine",
  "United Kingdom"
)

country_query <- paste0(
  "(",
  paste(paste0('"', european_countries, '"[Country]'), collapse = " OR "),
  ")"
)

build_query <- function(family) {
  paste0(
    '"', family, '"[Organism] AND ',
    '(', gene, ') AND ',
    '200[SLEN]:2000[SLEN] AND ',
    country_query
  )
}

queries <- setNames(lapply(families, build_query), families)
cat("Queries built for:", paste(families, collapse = ", "), "\n")

# --- Search & count ------------------------------------------------------
search_results <- list()
for (fam in families) {
  cat("\nSearching:", fam, "...\n")
  res <- entrez_search(
    db          = "nucleotide",
    term        = queries[[fam]],
    retmax      = max_records,
    use_history = TRUE
  )
  search_results[[fam]] <- res
  cat("  ->", res$count, "records found (fetching up to", max_records, ")\n")
  Sys.sleep(sleep_time)
}

# --- Fetch sequences in batches -------------------------------------
all_sequences <- list()
all_metadata  <- list()

for (fam in families) {
  res   <- search_results[[fam]]
  ids   <- res$ids
  n_ids <- length(ids)

  if (n_ids == 0) {
    cat("No sequences found for", fam, "-- skipping.\n")
    next
  }

  cat("\nFetching", n_ids, "sequences for", fam, "...\n")

  seq_list  <- character(0)
  meta_list <- list()
  starts    <- seq(1, n_ids, by = batch_size)

  for (s in starts) {
    end       <- min(s + batch_size - 1, n_ids)
    batch_ids <- ids[s:end]

    fasta_raw <- entrez_fetch(
      db      = "nucleotide",
      id      = batch_ids,
      rettype = "fasta",
      retmode = "text"
    )
    seq_list <- c(seq_list, fasta_raw)

    summ <- entrez_summary(db = "nucleotide", id = batch_ids)

    for (id in names(summ)) {
      s_item <- summ[[id]]
      meta_list[[id]] <- tibble(
        accession   = s_item$caption,
        gi          = id,
        title       = s_item$title,
        length_bp   = s_item$slen,
        organism    = s_item$organism,
        family      = fam,
        taxid       = s_item$taxid,
        update_date = s_item$updatedate
      )
    }

    cat("  Fetched batch", s, "-", end, "\n")
    Sys.sleep(sleep_time)
  }

  fasta_combined <- paste(seq_list, collapse = "\n")
  tmp_fasta <- tempfile(fileext = ".fasta")
  writeLines(fasta_combined, tmp_fasta)
  dna <- readDNAStringSet(tmp_fasta)

  all_sequences[[fam]] <- dna
  all_metadata[[fam]]  <- bind_rows(meta_list)
  cat("  ->", length(dna), "sequences parsed for", fam, "\n")
}

# --- Export ------------------------------------------------------------
all_dna <- Reduce(c, all_sequences)
cat("Total sequences:", length(all_dna), "\n")

writeXStringSet(all_dna, filepath = out_fasta_raw)
cat("Raw FASTA written to:", out_fasta_raw, "\n")

metadata_df <- bind_rows(all_metadata)
write_csv(metadata_df, out_metadata)
cat("Metadata written to:", out_metadata, "\n")

