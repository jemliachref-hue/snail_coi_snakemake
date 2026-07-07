# Snakemake Pipeline — COI Phylogenetics of European Slug Families

A reproducible Snakemake workflow for downloading, filtering, aligning, and visualising COI (cytochrome c oxidase subunit I) sequences of European slug families from NCBI GenBank.

Built as part of a bioinformatics retraining internship at the **Museum für Naturkunde Berlin** (April–June 2026), under the supervision of Dr. Jörg Freyhof.

---

## Biological Context

This pipeline targets four European slug families:

| Family | Colour in trees |
|---|---|
| Limacidae | Red `#E63946` |
| Milacidae | Teal `#2A9D8F` |
| Agriolimacidae | Yellow `#E9C46A` |
| Arionidae | Blue `#457B9D` |

An outgroup (*Helix pomatia*, Helicidae) is used for rooting in two of the five output trees.

---

## Pipeline Overview

The workflow consists of **11 Snakemake rules** chained in a DAG:

```
download ──► filter_length ──► clean_fasta ──┬──► align ──► nj_tree ──┬──► plot_full_tree
                                              │                        └──► plot_all_sequences
                                              ├──► subsample_genus ──► outgroup_heatmap ◄── download_outgroup
                                              └──► identified_species_tree ◄────────────── download_outgroup
```

| Rule | Script | Description |
|---|---|---|
| `download` | `01_download.R` | Fetch COI sequences from NCBI via `rentrez` |
| `filter_length` | `02_filter_length.R` | Remove sequences shorter than 400 bp |
| `clean_fasta` | `03_clean_fasta.R` | Standardise headers, remove ambiguous/duplicate sequences |
| `align` | `04_align.R` | Multiple sequence alignment with DECIPHER |
| `nj_tree` | `05_nj_tree.R` | K80 distance matrix + Neighbor-Joining tree |
| `plot_full_tree` | `06_plot_full_tree.R` | ggtree plot of all sequences |
| `subsample_genus` | `07_subsample_genus.R` | 1 sequence per genus + subtree |
| `download_outgroup` | `08_download_outgroup.R` | Fetch *Helix pomatia* from NCBI |
| `outgroup_heatmap` | `09_outgroup_heatmap.R` | Outgroup-rooted tree + K80 distance heatmap |
| `identified_species_tree` | `10_identified_species_tree.R` | Identified species only + GenBank accession labels |
| `plot_all_sequences` | `11_plot_all_sequences.R` | Zoomable PDF of all 1625 sequences |

---

## Project Structure

```
snail_coi_snakemake/
├── Snakefile                        # pipeline orchestration
├── scripts/
│   ├── 01_download.R
│   ├── 02_filter_length.R
│   ├── 03_clean_fasta.R
│   ├── 04_align.R
│   ├── 05_nj_tree.R
│   ├── 06_plot_full_tree.R
│   ├── 07_subsample_genus.R
│   ├── 08_download_outgroup.R
│   ├── 09_outgroup_heatmap.R
│   ├── 10_identified_species_tree.R
│   └── 11_plot_all_sequences.R
├── dag.png                          # DAG visualisation
├── .gitignore
└── Data/
    └── output_COI/                  # all outputs (gitignored)
```

---

## Requirements

### System
- Linux / WSL2 (Ubuntu)
- [Miniforge](https://github.com/conda-forge/miniforge) (conda/mamba)

### Installation

**1. Clone the repository:**
```bash
git clone https://github.com/jemliachref-hue/snail_coi_snakemake.git
cd snail_coi_snakemake
```

**2. Install Snakemake:**
```bash
conda install -c conda-forge -c bioconda snakemake-minimal
```

**3. Install R and required packages:**
```bash
# R itself
conda install -c conda-forge r-base

# CRAN packages
conda install -c conda-forge r-rentrez r-dplyr r-stringr r-tibble r-readr \
  r-ape r-phangorn r-ggplot2

# Bioconductor packages
conda install -c bioconda bioconductor-biostrings bioconductor-decipher \
  bioconductor-ggtree
```

> **Note:** Installing `ggtree` via conda (Bioconda) is strongly recommended over `BiocManager::install()` to avoid version mismatch errors with `treeio`.

---

## Running the Pipeline

```bash
# Dry run — see what will execute without running anything
snakemake -n

# Full run
snakemake --cores 1

# Run targeting a specific output file
snakemake --cores 1 Data/output_COI/snails_COI_NJ.nwk

# Force re-run of all rules
snakemake --cores 1 --forceall

# Visualise the DAG (requires graphviz)
snakemake --dag | dot -Tpng > dag.png
```

> Keep `--cores 1` — the NCBI download steps are rate-limited and should not run in parallel.

---

## Outputs

All outputs land in `Data/output_COI/` (excluded from git, reproducible by running the pipeline):

| File | Description |
|---|---|
| `snails_COI_Europe.fasta` | Raw sequences from NCBI (1648 sequences) |
| `snails_COI_metadata.csv` | Accession metadata (organism, family, length, taxid) |
| `snails_COI_Europe_filtered.fasta` | After length filtering (≥400 bp) |
| `snails_COI_clean.fasta` | After header standardisation + deduplication |
| `snails_COI_aligned.fasta` | DECIPHER multiple sequence alignment (709 columns) |
| `snails_COI_NJ.nwk` | Midpoint-rooted NJ tree (Newick format) |
| `snails_COI_NJ_tree.pdf` | Full tree — all sequences |
| `snails_COI_clean_subsampled.fasta` | 1 sequence per genus (17 genera) |
| `snails_COI_NJ_tree_genus.pdf` | Subtree — 1 sequence per genus |
| `outgroup_helix_pomatia.fasta` | *Helix pomatia* outgroup sequence |
| `snails_COI_NJ_outgroup.nwk` | Outgroup-rooted NJ tree |
| `snails_COI_NJ_outgroup_heatmap.pdf` | Outgroup tree + K80 distance heatmap |
| `snails_COI_NJ_identified.nwk` | Identified-species-only tree |
| `snails_COI_NJ_identified_species.pdf` | Identified species tree + GenBank accessions |
| `snails_COI_NJ_ALL_sequences.pdf` | Zoomable PDF of all 1625 sequences |

---

## Parameters

All tunable parameters are defined at the top of the `Snakefile`:

```python
FAMILIES    = ["Limacidae", "Milacidae", "Agriolimacidae", "Arionidae"]
GENE        = "COI[Gene] OR COX1[Gene] OR ..."
BATCH_SIZE  = 100      # sequences per NCBI fetch batch
MAX_RECORDS = 500      # max sequences per family
SLEEP_TIME  = 0.4      # seconds between NCBI requests
MIN_LENGTH  = 400      # minimum sequence length in bp
```

To change families or thresholds, edit these values and re-run `snakemake --cores 1` — only affected rules will re-execute.

---

## Author

**Dr. Achref Jemli**
PhD in Biological Sciences and Biotechnology
Retraining in Applied Bioinformatics and Biostatistics
Internship — Museum für Naturkunde Berlin, 2026
