# Construct a phylogenetic tree for the ASVs in this dataset using DECIPHER and phangorn.

set.seed(1)

# Load data from the ASV table output by DADA2.
# This should be an RDS file containing a matrix of ASV counts with ASV sequences as column names.
asv_table_candidates <- c(
    here::here("data", "processed", "dada2", "ASV_paired_end_table.rds"),
    here::here("outputs", "ASV_paired_end_table.rds"),
    here::here("data", "processed", "dada2", "single_end", "ASV_table_SE.rds"),
    here::here("outputs", "single_end", "ASV_table_SE.rds")
)

asv_table_path <- asv_table_candidates[file.exists(asv_table_candidates)][1]
if (is.na(asv_table_path)) {
    stop("No ASV table was found. Looked for paired-end and single-end RDS outputs in the processed and outputs directories.")
}

output_dir <- here::here("outputs", "phylogeny")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Directory for figures (separate from tree objects/newick)
figures_dir <- here::here("outputs", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
asv_table <- readRDS(asv_table_path)
asv_table <- as.matrix(asv_table)

if (is.null(colnames(asv_table))) {
    stop("The ASV table does not have column names, so the ASV sequences cannot be recovered.")
}

asv_sequences <- gsub("\\s+", "", colnames(asv_table))
if (anyDuplicated(asv_sequences) > 0L) {
    stop("Duplicate ASV sequences were found in the ASV table column names.")
}

# Convert ASV sequences to DNAStringSet for alignment and tree building.
asv_dna <- Biostrings::DNAStringSet(asv_sequences)
names(asv_dna) <- asv_sequences

cat("Loaded", nrow(asv_table), "samples and", ncol(asv_table), "ASVs from:", asv_table_path, "\n")

# Align the ASV sequences before distance estimation and tree inference.
aligned_asvs <- DECIPHER::AlignSeqs(asv_dna, processors = 1L)
phang.alignment <- phangorn::phyDat(as(aligned_asvs, "matrix"), type = "DNA")

# Use a model-based distance matrix for the NJ starting tree.
# The JC69 model is a simple substitution model that assumes equal base frequencies and
# equal substitution rates.
asv_distance_matrix <- phangorn::dist.ml(phang.alignment, model = "JC69")
saveRDS(asv_distance_matrix, file.path(output_dir, "asv_distance_matrix.rds"))

# Build a neighbor-joining tree as a starting point for maximum likelihood optimization.
nj_tree <- ape::nj(asv_distance_matrix)
nj_tree$edge.length[nj_tree$edge.length < 0] <- 0 # NJ can produce negative branch lengths,
# which are not biologically meaningful. Set them to zero.
saveRDS(nj_tree, file.path(output_dir, "asv_tree_nj.rds"))
ape::write.tree(nj_tree, file = file.path(output_dir, "asv_tree_nj.newick"))
# Save the NJ tree in Newick format for compatibility with other software.
cat("Neighbor-joining tree constructed and saved.\n")

# Refine the NJ tree with maximum likelihood using the aligned ASV sequences.
cat("Starting maximum likelihood optimization of the tree. This may take a while...\n")
fit_nj <- phangorn::pml(nj_tree, data = phang.alignment)
fit_ml <- phangorn::optim.pml(
    fit_nj,
    model = "GTR",
    optInv = TRUE,
    optGamma = TRUE,
    rearrangement = "stochastic",
    control = phangorn::pml.control(trace = 0) # Suppress optimization output for cleaner logs.
)
cat("Maximum likelihood tree optimization complete.\n")

# This step will take a while for large ASV datasets, especially with complex models.
# The GTR model is a general time-reversible model that allows for different substitution rates between all pairs of nucleotides,
# which can provide a better fit to the data than simpler models like JC69.
cat("Optimized tree log-likelihood:", fit_ml$logLik, "\n")
ml_tree <- fit_ml$tree
ml_tree$edge.length[ml_tree$edge.length < 0] <- 0
cat("Negative branch lengths set to zero.\n")

# Save the maximum likelihood tree and the fitted model for downstream analysis and visualization.
saveRDS(fit_ml, file.path(output_dir, "asv_tree_ml_fit.rds"))
saveRDS(ml_tree, file.path(output_dir, "asv_tree_ml.rds"))
ape::write.tree(ml_tree, file = file.path(output_dir, "asv_tree_ml.newick"))
cat("Maximum likelihood tree saved in RDS and Newick formats.\n")

# Optional branch support. Increase this value when you want bootstrap support values.
cat("Calculating bootstrap support values for the tree (this can be very time-consuming)...\n")
bootstrap_replicates <- 0L
if (bootstrap_replicates > 0L) {
    bootstrap_support <- phangorn::bootstrap.pml(
        fit_ml,
        bs = bootstrap_replicates,
        optNni = TRUE,
        multicore = TRUE,
        control = phangorn::pml.control(trace = 0)
    )
    saveRDS(bootstrap_support, file.path(output_dir, "asv_tree_bootstrap_support.rds"))
}

# The resulting tree can be used for phylogenetically informed analyses, such as UniFrac distance calculations,
# phylogenetic diversity metrics, and visualizations that incorporate evolutionary relationships among ASVs.

# Visualize the tree with ggtree (optional, requires ggtree package).
cat("Visualizing the tree with ggtree and optionally exporting an interactive Plotly HTML (if installed)...\n")
if (requireNamespace("ggtree", quietly = TRUE)) {
    library(ggtree)
    library(ggplot2)
    tree_plot <- ggtree(ml_tree) +
        geom_tiplab(size = 2) +
        theme_tree2() +
        ggtitle("Maximum Likelihood Phylogenetic Tree of ASVs")

    # Save static PNG to the figures directory
    ggsave(filename = file.path(figures_dir, "asv_tree_ml_ggtree.png"), plot = tree_plot, width = 8, height = 10)

    # If plotly is available, export an interactive HTML version
    if (requireNamespace("plotly", quietly = TRUE)) {
        if (exists("tree_plot") && inherits(tree_plot, "ggplot")) {
            tryCatch(
                {
                    interactive_plot <- plotly::ggplotly(tree_plot)
                    if (requireNamespace("htmlwidgets", quietly = TRUE)) {
                        htmlwidgets::saveWidget(htmlwidgets::as_widget(interactive_plot),
                            file = file.path(figures_dir, "asv_tree_ml_ggtree_interactive.html"),
                            selfcontained = FALSE
                        )
                        cat("Interactive Plotly HTML saved to figures directory.\n")
                    } else {
                        warning("htmlwidgets package not available; cannot save interactive HTML.")
                    }
                },
                error = function(e) {
                    warning("Failed to create or save interactive Plotly plot: ", conditionMessage(e))
                }
            )
        } else {
            warning("`tree_plot` not created or not a ggplot object; skipping interactive export.")
        }
    } else {
        cat("Plotly not installed; skipping interactive export.\n")
    }
}

cat("Phylogeny outputs saved to:", output_dir, "\n")

# Additional analyses could include:
# - Calculating UniFrac distances using the phylogenetic tree and ASV abundances.
# - Computing Faith's Phylogenetic Diversity for each sample.
# - Visualizing the tree with metadata annotations (e.g., sample groups, ASV abundances).

# =============================================================================
# PIPELINE COMPLETE
# =============================================================================
cat("\n")
cat("============================================\n")
cat("   PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("============================================\n")
cat("\n")

# =============================================================================
# END OF PIPELINE
# =============================================================================
