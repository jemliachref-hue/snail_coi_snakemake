## Rule: plot_full_tree
## Section 11 of original script.

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

tip_df <- tibble(label = nj_tree_rooted$tip.label) |>
  mutate(
    family   = str_extract(label, "[^|]+$"),
    organism = str_extract(label, "(?<=\\|)[^|]+(?=\\|)")
  )

family_colours <- c(
  "Limacidae"      = "#E63946",
  "Milacidae"      = "#2A9D8F",
  "Agriolimacidae" = "#E9C46A",
  "Arionidae"      = "#457B9D"
)

p <- ggtree(nj_tree_rooted, layout = "rectangular", linewidth = 0.3) %<+% tip_df +
  geom_tippoint(aes(color = family), size = 1.2) +
  geom_tiplab(
    aes(label = organism, color = family),
    size     = 2,
    offset   = 0.001,
    fontface = "italic"
  ) +
  scale_color_manual(values = family_colours, name = "Family") +
  theme_tree2() +
  theme(
    legend.position = "bottom",
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10),
    plot.title      = element_text(size = 14, face = "bold", hjust = 0.5),
    plot.subtitle   = element_text(size = 10, hjust = 0.5, color = "grey40")
  ) +
  labs(
    title    = "Neighbor-Joining Tree — Slug COI (Europe)",
    subtitle = "K80 model | Midpoint rooted | DECIPHER alignment"
  )

ggsave(filename = out_pdf, plot = p, width = 14, height = 22, device = "pdf")
cat("Tree plot saved to:", out_pdf, "\n")
