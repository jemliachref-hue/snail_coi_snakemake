## Rule: align
## Section 9 of original script.

library(Biostrings)
library(DECIPHER)

in_clean    <- snakemake@input[["fasta_clean"]]
out_aligned <- snakemake@output[["fasta_aligned"]]

clean_dna <- readDNAStringSet(in_clean)

cat("Running alignment...\n")
cat("(This may take a few minutes)\n\n")

aligned_dna <- AlignSeqs(
  DNAStringSet(clean_dna),
  anchor  = NA,
  verbose = TRUE
)

cat("\nAlignment done!\n")
cat("Alignment width:", unique(width(aligned_dna)), "columns\n")

writeXStringSet(aligned_dna, filepath = out_aligned, width = 200000L)
cat("Aligned FASTA written to:", out_aligned, "\n")