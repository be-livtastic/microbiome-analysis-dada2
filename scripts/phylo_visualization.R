############################################################
# MICROBIOME PHYLOGENETIC TREE VISUALIZATION WORKFLOW
# - Big-picture circular tree
# - Taxonomy-collapsed tree
# - Top abundant ASV subtree
# - Phylogenetic heatmap
############################################################

############################
# LOAD LIBRARIES
############################

library(phyloseq)
library(ggtree)
library(ggplot2)
library(ape)
library(phangorn)
library(dplyr)
library(pheatmap)
library(RColorBrewer)

############################
# LOAD DATA
############################

set.seed(1)

choose_file <- function(paths) {
    existing <- paths[file.exists(paths)]
    if (length(existing) == 0) {
        return(NA_character_)
    }
    existing[1]
}

output_dir <- here::here("outputs", "phylogeny")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Directory for figures (separate from tree objects/newick)
figures_dir <- here::here("outputs", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
# Tree
ml_tree <- readRDS(file.path(output_dir, "asv_tree_ml.rds"))

# ASV abundance table
asv_table_path <- choose_file(c(
    here::here("data", "processed", "dada2", "ASV_paired_end_table.rds"),
    here::here("outputs", "ASV_paired_end_table.rds")
))

if (is.na(asv_table_path)) {
    stop("No paired-end ASV table was found. Looked in data/processed/dada2 and outputs/.")
}

asv_table <- readRDS(asv_table_path)

# Taxonomy table
taxonomy_table_path <- choose_file(c(
    here::here("data", "processed", "dada2", "taxonomy.rds"),
    here::here("outputs", "taxonomy.rds"),
    here::here("outputs", "single_end", "taxonomy_SE.rds"),
    here::here("outputs", "single_end", "taxonomy_SE_genuslevel.rds")
))

taxonomy_table <- NULL
if (!is.na(taxonomy_table_path)) {
    taxonomy_table <- readRDS(taxonomy_table_path)
}
############################
# CREATE PHYLOSEQ OBJECT (with robust taxonomy handling)
############################

# Basic otu and tree objects (ASV table assumed samples x ASVs)
otu <- otu_table(as.matrix(asv_table), taxa_are_rows = FALSE)
tree <- phy_tree(ml_tree)

# Decide whether to use single-end taxonomy: only accept single-end taxonomy
# when the ASV table came from a single-end run. If a single-end taxonomy file
# is present but the ASV table is paired-end, ignore the single-end taxonomy.
if (!is.null(taxonomy_table) && grepl("single_end", taxonomy_table_path, ignore.case = TRUE)) {
    if (!grepl("single_end", asv_table_path, ignore.case = TRUE) && !grepl("_SE", basename(asv_table_path), ignore.case = TRUE)) {
        cat("Found single-end taxonomy but ASV table appears paired-end; ignoring single-end taxonomy.\n")
        taxonomy_table <- NULL
    }
}

# If taxonomy is present, validate that ASV identifiers match the ASV table.
if (!is.null(taxonomy_table)) {
    # ASV identifiers in the pipeline are stored as column names in the ASV table
    asv_ids <- colnames(asv_table)
    tax_ids <- rownames(taxonomy_table)

    common_ids <- intersect(asv_ids, tax_ids)

    if (length(common_ids) == 0) {
        cat("Taxonomy rownames do not match ASV identifiers; ignoring taxonomy to avoid invalid phyloseq object.\n")
        taxonomy_table <- NULL
        ps <- phyloseq(otu, tree)
    } else {
        if (length(common_ids) < length(asv_ids)) {
            cat(sprintf("Partial match: keeping %d/%d ASVs that have taxonomy assignments.\n", length(common_ids), length(asv_ids)))
            # Subset the ASV table to the intersection so components align
            asv_table_sub <- asv_table[, common_ids, drop = FALSE]
            otu <- otu_table(as.matrix(asv_table_sub), taxa_are_rows = FALSE)

            # Subset taxonomy to the matching ASVs and reorder to match otu taxa order
            taxonomy_table <- taxonomy_table[common_ids, , drop = FALSE]

            # Prune the tree to the matching tips as well (avoid mismatched tip labels)
            tips_to_keep <- common_ids
            # Some trees use slightly different tip label formatting; ensure intersection
            tree_tips <- ml_tree$tip.label
            tips_present <- intersect(tree_tips, tips_to_keep)
            if (length(tips_present) < length(tips_to_keep)) {
                cat("Warning: some ASVs with taxonomy are not present in the tree; pruning to the shared set.\n")
                tips_to_keep <- tips_present
                # also subset otu and taxonomy accordingly
                otu <- otu_table(as.matrix(asv_table_sub[, tips_to_keep, drop = FALSE]), taxa_are_rows = FALSE)
                taxonomy_table <- taxonomy_table[tips_to_keep, , drop = FALSE]
            }

            if (length(tips_to_keep) == 0) {
                stop("No overlapping ASV identifiers between ASV table, taxonomy and tree. Cannot build phyloseq object.")
            }

            pruned_tree <- drop.tip(ml_tree, setdiff(ml_tree$tip.label, tips_to_keep))
            tree <- phy_tree(pruned_tree)

            tax <- tax_table(as.matrix(taxonomy_table))
            ps <- phyloseq(otu, tax, tree)
        } else {
            # Full match — simply create phyloseq object with taxonomy
            tax <- tax_table(as.matrix(taxonomy_table))
            ps <- phyloseq(otu, tax, tree)
        }
    }
} else {
    ps <- phyloseq(otu, tree)
    cat("No taxonomy table was found; family and genus collapsed trees will be skipped.\n")
}

############################################################
# 1. BIG-PICTURE TREE
# Circular tree WITHOUT labels
############################################################

p_big <- ggtree(
    phy_tree(ps),
    layout = "circular",
    linewidth = 0.2
) +
    geom_tippoint(size = 0.3, alpha = 0.6) +
    ggtitle("Full ASV Phylogenetic Tree") +
    theme_void()

print(p_big)

ggsave(
    filename = file.path(
        figures_dir,
        "full_asv_tree_circular.png"
    ),
    plot = p_big,
    width = 12,
    height = 12,
    dpi = 300
)

############################################################
# 2. COLLAPSE TREE BY FAMILY
# Biological interpretation becomes easier
############################################################

if (!is.null(taxonomy_table)) {
    # Merge ASVs belonging to same Family
    ps_family <- tax_glom(ps, taxrank = "Family")

    # Remove taxa without family assignment
    ps_family <- prune_taxa(
        !is.na(as.vector(tax_table(ps_family)[, "Family"])),
        ps_family
    )

    family_tree <- phy_tree(ps_family)

    p_family <- ggtree(
        family_tree,
        layout = "circular",
        linewidth = 0.4
    ) +
        geom_tiplab(size = 2) +
        ggtitle("Family-Level Collapsed Tree") +
        theme_tree()

    print(p_family)

    ggsave(
        filename = file.path(
            figures_dir,
            "family_collapsed_tree.png"
        ),
        plot = p_family,
        width = 14,
        height = 14,
        dpi = 300
    )
}

############################################################
# 3. COLLAPSE TREE BY GENUS
############################################################

if (!is.null(taxonomy_table)) {
    ps_genus <- tax_glom(ps, taxrank = "Genus")

    ps_genus <- prune_taxa(
        !is.na(as.vector(tax_table(ps_genus)[, "Genus"])),
        ps_genus
    )

    genus_tree <- phy_tree(ps_genus)

    p_genus <- ggtree(
        genus_tree,
        layout = "circular",
        linewidth = 0.4
    ) +
        geom_tiplab(size = 1.8) +
        ggtitle("Genus-Level Collapsed Tree") +
        theme_tree()

    print(p_genus)

    ggsave(
        filename = file.path(
            figures_dir,
            "genus_collapsed_tree.png"
        ),
        plot = p_genus,
        width = 14,
        height = 14,
        dpi = 300
    )
}

############################################################
# 4. TOP 50 MOST ABUNDANT ASVs
# Best for thesis/manuscript figures
############################################################

# Total abundance per ASV
asv_abundance <- taxa_sums(ps)

# Top 50 ASVs
top50_taxa <- names(sort(
    asv_abundance,
    decreasing = TRUE
))[seq_len(min(50, length(asv_abundance)))]

# Prune phyloseq object
ps_top50 <- prune_taxa(top50_taxa, ps)

top50_tree <- phy_tree(ps_top50)

p_top50 <- ggtree(
    top50_tree,
    layout = "rectangular",
    linewidth = 0.5
) +
    geom_tiplab(size = 2) +
    ggtitle("Top 50 Most Abundant ASVs") +
    theme_tree2()

print(p_top50)

ggsave(
    filename = file.path(
        figures_dir,
        "top50_asv_tree.png"
    ),
    plot = p_top50,
    width = 14,
    height = 10,
    dpi = 300
)

############################################################
# 5. PHYLOGENETIC HEATMAP
# Much better for comparing samples
############################################################

top_n_asvs_heatmap <- 30

# Keep the display focused on the most common ASVs so the heatmap stays readable.
asv_abundance <- taxa_sums(ps)
top_asvs_heatmap <- names(
    sort(
        asv_abundance,
        decreasing = TRUE
    )
)[seq_len(min(top_n_asvs_heatmap, length(asv_abundance)))]

ps_top_heatmap <- prune_taxa(
    top_asvs_heatmap,
    ps
)

# Convert counts to within-sample relative abundance so samples can be compared on the same scale.
# A sample with more sequencing depth should not automatically look more intense.
ps_top_heatmap_rel <- transform_sample_counts(
    ps_top_heatmap,
    function(x) {
        total <- sum(x)
        if (total == 0) {
            return(x)
        }
        x / total
    }
)

# pheatmap expects a plain matrix; this extracts the ASV-by-sample abundance table.
heatmap_matrix <- as.matrix(
    otu_table(ps_top_heatmap_rel)
)

# Phyloseq can store taxa in rows or columns depending on the input object, so make the layout explicit.
if (!taxa_are_rows(ps_top_heatmap_rel)) {
    heatmap_matrix <- t(heatmap_matrix)
}

# Log scaling compresses large values and makes low-abundance ASVs visible without dominating the color range.
# The small pseudocount prevents log10(0).
heatmap_matrix <- log10(
    heatmap_matrix + 1e-6
)

# Use short labels in the plot because the original ASV sequences are too long to read on a figure.
short_asv_names <- paste0(
    "ASV_",
    seq_len(nrow(heatmap_matrix))
)

rownames(heatmap_matrix) <- short_asv_names

annotation_rows <- NULL
if (!is.null(taxonomy_table)) {
    taxonomy_subset <- as.data.frame(
        taxonomy_table
    )[top_asvs_heatmap, , drop = FALSE]

    # Row annotations show the Family and Genus for each ASV.
    # Unassigned values are labeled explicitly so missing taxonomy is easy to spot.
    annotation_rows <- data.frame(
        Family = ifelse(
            is.na(taxonomy_subset$Family) | taxonomy_subset$Family == "",
            "Unassigned",
            as.character(taxonomy_subset$Family)
        ),
        Genus = ifelse(
            is.na(taxonomy_subset$Genus) | taxonomy_subset$Genus == "",
            "Unassigned",
            as.character(taxonomy_subset$Genus)
        )
    )

    rownames(annotation_rows) <- short_asv_names
}

# Column annotations show library size after pruning to the selected ASVs.
# This is useful for spotting samples with unusually low remaining signal.
sample_depths <- sample_sums(ps_top_heatmap)
annotation_cols <- data.frame(
    Read_depth = sample_depths
)
rownames(annotation_cols) <- names(sample_depths)

pheatmap::pheatmap(
    heatmap_matrix,

    # Cluster rows to group ASVs with similar abundance profiles across samples.
    cluster_rows = TRUE,
    # Cluster columns to group samples with similar ASV composition.
    cluster_cols = TRUE,
    # Remove cell borders so the color pattern is easier to read.
    border_color = NA,
    # Keep the ASV IDs and sample names visible; the short labels make this practical.
    show_rownames = TRUE,
    show_colnames = TRUE,
    # Rotate sample labels so they remain legible when there are many samples.
    angle_col = 45,
    # Slightly smaller font for dense figures.
    fontsize_row = 7,
    fontsize_col = 8,
    # Light gray for any missing values after processing.
    na_col = "#f2f2f2",
    # Scale each row so the plot emphasizes relative patterns across samples.
    # This makes it easier to see which samples are enriched for a given ASV.
    scale = "row",
    # Diverging palette highlights enrichment and depletion around the row mean.
    color = colorRampPalette(
        rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu"))
    )(100),
    # Family and Genus help interpret whether clustered ASVs belong to related taxa.
    annotation_row = annotation_rows,
    # Read-depth annotation helps explain whether sample-level clustering reflects biology or sequencing depth.
    annotation_col = annotation_cols,
    # Title communicates that this is a normalized, log-scaled summary of the top ASVs.
    main = paste("Improved Phylogenetic Heatmap of Top", top_n_asvs_heatmap, "ASVs"),
    # Save directly to file so the figure is reproducible and does not depend on the active plotting device.
    filename = file.path(
        figures_dir,
        "asv_tree_ml_phylogenetic_heatmap_top_asvs_improved.png"
    ),
    width = 14,
    height = 12
)

############################################################
# 6. OPTIONAL:
# Midpoint-root the tree for cleaner visualization
############################################################

rooted_tree <- midpoint(ml_tree)

p_rooted <- ggtree(
    rooted_tree,
    layout = "fan"
) +
    geom_tippoint(size = 0.3) +
    ggtitle("Midpoint-Rooted Fan Tree") +
    theme_tree()

print(p_rooted)

ggsave(
    filename = file.path(
        figures_dir,
        "midpoint_rooted_fan_tree.png"
    ),
    plot = p_rooted,
    width = 12,
    height = 12,
    dpi = 300
)

############################################################
# OUTPUT SUMMARY
############################################################

cat("\nVisualization workflow complete.\n")

cat("\nSaved figures:\n")

cat("- full_asv_tree_circular.png\n")
cat("- family_collapsed_tree.png\n")
cat("- genus_collapsed_tree.png\n")
cat("- top50_asv_tree.png\n")
cat("- asv_tree_ml_phylogenetic_heatmap_top_asvs_improved.png\n")
cat("- midpoint_rooted_fan_tree.png\n")
