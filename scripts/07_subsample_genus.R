## Rule: subsample_genus
## Section 12 of original script.

library(Biostrings)
library(DECIPHER)
library(ape)
library(phangorn)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
options(ignore.negative.edge = TRUE)

in_clean      <- snakemake@input[["fasta_clean"]]
out_pdf       <- snakemake@output[["pdf"]]
out_fasta_sub <- snakemake@output[["fasta_sub"]]

clean_dna <- readDNAStringSet(in_clean)

family_colours <- c(
  "Limacidae"      = "#E63946",
  "Milacidae"      = "#2A9D8F",
  "Agriolimacidae" = "#E9C46A",
  "Arionidae"      = "#457B9D"
)

species_df2 <- tibble(
  name       = names(clean_dna),
  organism   = str_extract(names(clean_dna), "(?<=\\|)[^|]+(?=\\|)"),
  family     = str_extract(names(clean_dna), "[^|]+$"),
  seq_length = width(clean_dna)
) |>
  mutate(genus = str_extract(organism, "^[^_]+")) |>
  group_by(genus) |>
  slice_max(seq_length, n = 1, with_ties = FALSE) |>
  ungroup()

clean_dna_sub <- clean_dna[species_df2$name]
cat("Sequences after subsampling (1 per genus):", length(clean_dna_sub), "\n")

writeXStringSet(clean_dna_sub, filepath = out_fasta_sub, width = 200000L)
cat("Subsampled FASTA written to:", out_fasta_sub, "\n")

aligned_sub <- AlignSeqs(DNAStringSet(clean_dna_sub), anchor = NA, verbose = FALSE)

dna_bin_sub <- as.DNAbin(aligned_sub)
dist_sub    <- dist.dna(dna_bin_sub, model = "K80", pairwise.deletion = TRUE)
nj_sub      <- nj(dist_sub)
nj_sub_rooted <- midpoint(nj_sub)
cat("Tips:", Ntip(nj_sub_rooted), "\n")

tip_df_sub <- tibble(label = nj_sub_rooted$tip.label) |>
  mutate(
    family   = str_extract(label, "[^|]+$"),
    organism = str_extract(label, "(?<=\\|)[^|]+(?=\\|)") |>
      str_replace_all("_", " ")
  )

p_sub <- ggtree(nj_sub_rooted, layout = "rectangular", linewidth = 0.5) %<+% tip_df_sub +
  geom_tippoint(aes(color = family), size = 3) +
  geom_tiplab(
    aes(label = organism, color = family),
    size     = 3.5,
    offset   = 0.002,
    fontface = "italic"
  ) +
  scale_color_manual(values = family_colours, name = "Family") +
  theme_tree2() +
  xlim(0, 0.22) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Neighbor-Joining Tree — Slug COI (Europe)",
    subtitle = "K80 model | Midpoint rooted | 1 sequence per genus"
  )

ggsave(filename = out_pdf, plot = p_sub, width = 10, height = 10, device = "pdf")
cat("Subsampled tree saved to:", out_pdf, "\n")
