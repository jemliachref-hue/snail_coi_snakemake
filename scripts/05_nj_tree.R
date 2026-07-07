## Rule: nj_tree
## Section 10 of original script.

library(Biostrings)
library(ape)
library(phangorn)
options(ignore.negative.edge = TRUE)

in_aligned <- snakemake@input[["fasta_aligned"]]
out_tree   <- snakemake@output[["tree"]]

aligned_dna <- readDNAStringSet(in_aligned)
dna_bin     <- as.DNAbin(aligned_dna)

cat("Computing K80 distance matrix...\n")
dist_matrix <- dist.dna(
  dna_bin,
  model             = "K80",
  pairwise.deletion = TRUE
)
cat("Distance range:", round(min(dist_matrix, na.rm = TRUE), 4),
    "-", round(max(dist_matrix, na.rm = TRUE), 4), "\n")

nj_tree        <- nj(dist_matrix)
nj_tree_rooted <- midpoint(nj_tree)

cat("NJ tree built!\n")
cat("Tips  :", Ntip(nj_tree_rooted), "\n")
cat("Nodes :", Nnode(nj_tree_rooted), "\n")

write.tree(nj_tree_rooted, file = out_tree)
cat("Tree saved to:", out_tree, "\n")