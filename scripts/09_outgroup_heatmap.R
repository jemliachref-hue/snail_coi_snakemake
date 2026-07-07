## Rule: outgroup_heatmap
## Sections 13.2-13.7 of original script.

library(Biostrings)
library(DECIPHER)
library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
options(ignore.negative.edge = TRUE)

in_fasta_sub  <- snakemake@input[["fasta_sub"]]
in_outgroup   <- snakemake@input[["outgroup_fasta"]]
out_tree      <- snakemake@output[["tree"]]
out_pdf       <- snakemake@output[["pdf"]]

clean_dna_sub <- readDNAStringSet(in_fasta_sub)
outgroup_dna  <- readDNAStringSet(in_outgroup)

clean_dna_out <- c(clean_dna_sub, outgroup_dna)
cat("Total sequences with outgroup:", length(clean_dna_out), "\n")

aligned_out <- AlignSeqs(DNAStringSet(clean_dna_out), anchor = NA, verbose = FALSE)
cat("Alignment done — width:", unique(width(aligned_out)), "columns\n")

dna_bin_out   <- as.DNAbin(aligned_out)
dist_out      <- dist.dna(dna_bin_out, model = "K80", pairwise.deletion = TRUE)
nj_out        <- nj(dist_out)
nj_out_rooted <- root(nj_out, outgroup = "Outgroup|Helix_pomatia|Helicidae", resolve.root = TRUE)

cat("Tips  :", Ntip(nj_out_rooted), "\n")
cat("Nodes :", Nnode(nj_out_rooted), "\n")

write.tree(nj_out_rooted, file = out_tree)
cat("Tree saved to:", out_tree, "\n")

tip_df_out <- tibble(label = nj_out_rooted$tip.label) |>
  mutate(family = str_extract(label, "[^|]+$")) |>
  mutate(organism = str_extract(label, "(?<=\\|)[^|]+(?=\\|)")) |>
  mutate(organism = str_replace_all(organism, "_", " "))

family_colours_out <- c(
  "Limacidae"      = "#E63946",
  "Milacidae"      = "#2A9D8F",
  "Agriolimacidae" = "#E9C46A",
  "Arionidae"      = "#457B9D",
  "Helicidae"      = "#6A0572"
)

dist_df <- as.matrix(dist_out)

p_out2 <- ggtree(nj_out_rooted, layout = "rectangular", linewidth = 0.6) %<+% tip_df_out +
  geom_tippoint(aes(color = family), size = 3) +
  scale_color_manual(values = family_colours_out, name = "Family") +
  theme_tree2() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 9, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Neighbor-Joining Tree — Slug COI (Europe)",
    subtitle = "K80 model | Helix pomatia outgroup | 1 sequence per genus"
  )

p_final <- gheatmap(
  p_out2,
  as.data.frame(dist_df),
  offset       = 0.03,
  width        = 0.6,
  colnames     = FALSE,
  legend_title = "K80 distance"
) +
  scale_fill_gradient(low = "#FFFFFF", high = "#E63946", name = "K80 distance") +
  geom_tiplab(
    aes(label = organism, color = family),
    size     = 3,
    offset   = 0.25,
    fontface = "italic"
  ) +
  hexpand(0.6)

ggsave(filename = out_pdf, plot = p_final, width = 14, height = 10, device = "pdf")
cat("Final plot saved to:", out_pdf, "\n")
