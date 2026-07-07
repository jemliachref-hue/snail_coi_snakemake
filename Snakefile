OUT = "Data/output_COI"

FAMILIES    = ["Limacidae", "Milacidae", "Agriolimacidae", "Arionidae"]
GENE        = "COI[Gene] OR COX1[Gene] OR cytochrome oxidase subunit I[Title] OR cytochrome c oxidase subunit I[Title]"
BATCH_SIZE  = 100
MAX_RECORDS = 500
SLEEP_TIME  = 0.4

rule all:
    input:
        f"{OUT}/snails_COI_NJ_tree.pdf",
        f"{OUT}/snails_COI_NJ_tree_genus.pdf",
        f"{OUT}/snails_COI_NJ_outgroup_heatmap.pdf",
        f"{OUT}/snails_COI_NJ_identified_species.pdf",
        f"{OUT}/snails_COI_NJ_ALL_sequences.pdf",


rule download:
    output:
        fasta_raw = f"{OUT}/snails_COI_Europe.fasta",
        metadata  = f"{OUT}/snails_COI_metadata.csv",
    params:
        families    = FAMILIES,
        gene        = GENE,
        batch_size  = BATCH_SIZE,
        max_records = MAX_RECORDS,
        sleep_time  = SLEEP_TIME,
    script:
        "scripts/01_download.R"

MIN_LENGTH = 400
rule filter_length:
    input:
        fasta_raw = rules.download.output.fasta_raw,
    output:
        fasta_filtered = f"{OUT}/snails_COI_Europe_filtered.fasta",
    params:
        min_length = MIN_LENGTH,
    script:
        "scripts/02_filter_length.R"

rule clean_fasta:
    input:
        fasta_filtered = rules.filter_length.output.fasta_filtered,
        metadata       = rules.download.output.metadata,
    output:
        fasta_clean = f"{OUT}/snails_COI_clean.fasta",
    script:
        "scripts/03_clean_fasta.R"

rule align:
    input:
        fasta_clean = rules.clean_fasta.output.fasta_clean,
    output:
        fasta_aligned = f"{OUT}/snails_COI_aligned.fasta",
    script:
        "scripts/04_align.R"

rule nj_tree:
    input:
        fasta_aligned = rules.align.output.fasta_aligned,
    output:
        tree = f"{OUT}/snails_COI_NJ.nwk",
    script:
        "scripts/05_nj_tree.R"

rule plot_full_tree:
    input:
        tree = rules.nj_tree.output.tree,
    output:
        pdf = f"{OUT}/snails_COI_NJ_tree.pdf",
    script:
        "scripts/06_plot_full_tree.R"

rule subsample_genus:
    input:
        fasta_clean = rules.clean_fasta.output.fasta_clean,
    output:
        pdf       = f"{OUT}/snails_COI_NJ_tree_genus.pdf",
        fasta_sub = f"{OUT}/snails_COI_clean_subsampled.fasta",
    script:
        "scripts/07_subsample_genus.R"

rule download_outgroup:
    output:
        outgroup_fasta = f"{OUT}/outgroup_helix_pomatia.fasta",
    script:
        "scripts/08_download_outgroup.R"

rule outgroup_heatmap:
    input:
        fasta_sub      = rules.subsample_genus.output.fasta_sub,
        outgroup_fasta = rules.download_outgroup.output.outgroup_fasta,
    output:
        tree = f"{OUT}/snails_COI_NJ_outgroup.nwk",
        pdf  = f"{OUT}/snails_COI_NJ_outgroup_heatmap.pdf",
    script:
        "scripts/09_outgroup_heatmap.R"

rule identified_species_tree:
    input:
        fasta_clean    = rules.clean_fasta.output.fasta_clean,
        outgroup_fasta = rules.download_outgroup.output.outgroup_fasta,
    output:
        tree = f"{OUT}/snails_COI_NJ_identified.nwk",
        pdf  = f"{OUT}/snails_COI_NJ_identified_species.pdf",
    script:
        "scripts/10_identified_species_tree.R"

rule plot_all_sequences:
    input:
        tree = rules.nj_tree.output.tree,
    output:
        pdf = f"{OUT}/snails_COI_NJ_ALL_sequences.pdf",
    script:
        "scripts/11_plot_all_sequences.R"
