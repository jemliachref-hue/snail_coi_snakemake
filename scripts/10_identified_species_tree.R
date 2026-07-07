## Rule: identified_species_tree
## Section 14 of original script.

library(Biostrings)
library(DECIPHER)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
options(ignore.negative.edge = TRUE)

in_clean    <- snakemake@input[["fasta_clean"]]
in_outgroup <- snakemake@input[["outgroup_fasta"]]
out_tree    <- snakemake@output[["tree"]]
out_pdf     <- snakemake@output[["pdf"]]

clean_dna    <- readDNAStringSet(in_clean)
outgroup_dna <- readDNAStringSet(in_outgroup)

family_colours_out <- c(
  "Limacidae"      = "#E63946",
  "Milacidae"      = "#2A9D8F",
  "Agriolimacidae" = "#E9C46A",
  "Arionidae"      = "#457B9D",
  "Helicidae"      = "#6A0572"
)

identified <- clean_dna[!grepl("_sp\\.", names(clean_dna))]
cat("After removing sp. sequences:", length(identified), "\n")

species_df3 <- tibble(
  name       = names(identified),
  organism   = str_extract(names(identified), "(?<=\\|)[^|]+(?=\\|)"),
  family     = str_extract(names(identified), "[^|]+$"),
  accession  = str_extract(names(identified), "^[^|]+"),
  seq_length = width(identified)
) |>
  group_by(organism) |>
  slice_max(seq_length, n = 1, with_ties = FALSE) |>
  ungroup()

identified_sub <- identified[species_df3$name]
cat("Unique identified species:", length(identified_sub), "\n")

identified_out <- c(identified_sub, outgroup_dna)
cat("Total with outgroup:", length(identified_out), "\n")

aligned_id <-
  AlignSeqs(DNAStringSet(identified_out), anchor = NA, verbose = FALSE)
cat("Alignment width:", unique(width(aligned_id)), "columns\n")

dna_bin_id   <- as.DNAbin(aligned_id)
dist_id      <- dist.dna(dna_bin_id, model = "K80", pairwise.deletion = TRUE)
nj_id        <- nj(dist_id)
nj_id_rooted <- root(nj_id, outgroup = "Outgroup|Helix_pomatia|Helicidae", resolve.root = TRUE)
cat("Tips:", Ntip(nj_id_rooted), "\n")
write.tree(nj_id_rooted, file = out_tree)
cat("Tree saved.\n")

tip_df_id <- tibble(label = nj_id_rooted$tip.label) |>
  mutate(
    accession = str_extract(label, "^[^|]+"),
    family    = str_extract(label, "[^|]+$"),
    organism  = str_extract(label, "(?<=\\|)[^|]+(?=\\|)")
  ) |>
  mutate(
    organism  = str_replace_all(organism, "_", " "),
    tip_label = ifelse(
      family == "Helicidae",
      paste0("Helix pomatia [outgroup]"),
      paste0(organism, "  [", accession, "]")
    )
  )

p_id <- ggtree(nj_id_rooted, layout = "rectangular", linewidth = 0.5) %<+% tip_df_id +
  geom_tippoint(aes(color = family), size = 2.5) +
  geom_tiplab(
    aes(label = tip_label, color = family),
    size     = 2.8,
    offset   = 0.001,
    fontface = "italic"
  ) +
  scale_color_manual(values = family_colours_out, name = "Family") +
  theme_tree2() +
  xlim(0, 0.35) +
  theme(
    legend.position  = "bottom",
    legend.title     = element_text(size = 12, face = "bold"),
    legend.text      = element_text(size = 11),
    plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Neighbor-Joining Tree — Slug COI (Europe)",
    subtitle = "K80 model | Helix pomatia outgroup | Identified species only"
  )

ggsave(filename = out_pdf, plot = p_id, width = 14, height = 26, device = "pdf")
cat("Tree saved to:", out_pdf, "\n")
