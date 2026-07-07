## Rule: plot_all_sequences
## Section 15 of original script.

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)
library(stringr)
library(tibble)
options(ignore.negative.edge = TRUE)

in_tree <- snakemake@input[["tree"]]
out_pdf <- snakemake@output[["pdf"]]

nj_tree_rooted <- read.tree(in_tree)

family_colours <- c(
  "Limacidae"      = "#E63946",
  "Milacidae"      = "#2A9D8F",
  "Agriolimacidae" = "#E9C46A",
  "Arionidae"      = "#457B9D"
)

tip_df_all <- tibble(label = nj_tree_rooted$tip.label) |>
  mutate(
    accession = str_extract(label, "^[^|]+"),
    family    = str_extract(label, "[^|]+$"),
    organism  = str_extract(label, "(?<=\\|)[^|]+(?=\\|)")
  ) |>
  mutate(
    organism  = str_replace_all(organism, "_", " "),
    tip_label = paste0(organism, "  [", accession, "]")
  )

p_all <- ggtree(nj_tree_rooted, layout = "rectangular", linewidth = 0.2) %<+% tip_df_all +
  geom_tiplab(
    aes(label = tip_label, color = family),
    size     = 1.2,
    offset   = 0.001,
    fontface = "italic",
    align    = TRUE,
    linetype = "solid",
    linesize = 0.2
  ) +
  scale_color_manual(values = family_colours, name = "Family") +
  theme_tree2() +
  xlim(0, 0.35) +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 11),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Neighbor-Joining Tree — Slug COI (Europe) — All Sequences",
    subtitle = "K80 model | Midpoint rooted | DECIPHER alignment | 1625 sequences"
  )

ggsave(
  filename  = out_pdf,
  plot      = p_all,
  width     = 20,
  height    = 120,
  device    = "pdf",
  limitsize = FALSE
)
cat("Full tree saved to:", out_pdf, "\n")
