## Rule: clean_fasta
## Section 8 of original script.

library(Biostrings)
library(dplyr)
library(stringr)
library(tibble)
library(readr)

in_fasta_filt <- snakemake@input[["fasta_filtered"]]
in_metadata   <- snakemake@input[["metadata"]]
out_clean     <- snakemake@output[["fasta_clean"]]

all_dna_filtered <- readDNAStringSet(in_fasta_filt)
metadata_df      <- read_csv(in_metadata, show_col_types = FALSE)

# 8.1 Standardise headers
extract_accession <- function(header) {
  str_extract(header, "^[^\\s.]+")
}

raw_accessions <- extract_accession(names(all_dna_filtered))

header_df <- tibble(
  accession_full = names(all_dna_filtered),
  accession      = str_remove(raw_accessions, "\\.\\d+$")
) |>
  left_join(
    metadata_df |> select(accession, organism, family),
    by = "accession"
  ) |>
  mutate(
    organism_clean = str_replace_all(organism, " ", "_"),
    new_header     = paste0(accession, "|", organism_clean, "|", family)
  )

clean_dna <- all_dna_filtered
names(clean_dna) <- header_df$new_header
cat("Clean header example:\n ", names(clean_dna)[1], "\n")

# 8.2 Remove ambiguous sequences
count_ambiguous <- function(seq) {
  bases <- strsplit(as.character(seq), "")[[1]]
  valid <- c("A", "T", "G", "C", "-")
  sum(!bases %in% valid) / length(bases)
}

ambig_prop <- sapply(clean_dna, count_ambiguous)
keep       <- ambig_prop <= 0.05

cat("Sequences removed (>5% ambiguous):", sum(!keep), "\n")
cat("Sequences kept                    :", sum(keep), "\n")

clean_dna <- clean_dna[keep]

# 8.3 Remove duplicates
clean_dna <- clean_dna[!duplicated(names(clean_dna))]
cat("Sequences after deduplication:", length(clean_dna), "\n")

# 8.4 Export clean FASTA (single-line)
writeXStringSet(clean_dna, filepath = out_clean, width = 200000L)
cat("Clean FASTA written to:", out_clean, "\n")
cat("Total sequences        :", length(clean_dna), "\n")