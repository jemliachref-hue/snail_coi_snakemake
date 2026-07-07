## Rule: filter_length
## Section 7 of original script.

library(Biostrings)

in_fasta_raw   <- snakemake@input[["fasta_raw"]]
out_fasta_filt <- snakemake@output[["fasta_filtered"]]
min_length     <- snakemake@params[["min_length"]]

all_dna <- readDNAStringSet(in_fasta_raw)

cat("Before filtering:", length(all_dna), "sequences\n")
cat("Length range    :", min(width(all_dna)), "-", max(width(all_dna)), "bp\n")

all_dna_filtered <- all_dna[width(all_dna) >= min_length]
cat("After filtering :", length(all_dna_filtered), "sequences (>=", min_length, "bp)\n")

writeXStringSet(all_dna_filtered, filepath = out_fasta_filt)
cat("Filtered FASTA written to:", out_fasta_filt, "\n")