## Rule: download_outgroup
## Section 13.1 of original script.

library(rentrez)
library(Biostrings)

out_fasta <- snakemake@output[["outgroup_fasta"]]

cat("Downloading Helix pomatia COI sequences...\n")

outgroup_search <- entrez_search(
  db     = "nucleotide",
  term   = '"Helix pomatia"[Organism] AND (COI[Gene] OR COX1[Gene]) AND 500[SLEN]:700[SLEN]',
  retmax = 5
)

outgroup_fasta <- entrez_fetch(
  db      = "nucleotide",
  id      = outgroup_search$ids,
  rettype = "fasta",
  retmode = "text"
)

tmp_out <- tempfile(fileext = ".fasta")
writeLines(outgroup_fasta, tmp_out)
outgroup_dna <- readDNAStringSet(tmp_out)

cat("Outgroup sequences downloaded:", length(outgroup_dna), "\n")

outgroup_dna <- outgroup_dna[which.max(width(outgroup_dna))]
names(outgroup_dna) <- "Outgroup|Helix_pomatia|Helicidae"
cat("Outgroup selected:", names(outgroup_dna), "\n")

writeXStringSet(outgroup_dna, filepath = out_fasta)
cat("Outgroup FASTA written to:", out_fasta, "\n")
