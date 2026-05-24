# Construct a phylogenetic tree for the ASVs in this dataset using DECIPHER and phangorn.

set.seed(1)

cat("\n============================================================\n")
cat("   PHYLOGENETIC TREE CONSTRUCTION\n")
cat("   DECIPHER alignment · NJ tree · ML optimisation (GTR+G)\n")
cat("============================================================\n\n")

output_dir <- here::here("outputs", "phylogeny")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Directory for figures (separate from tree objects/newick)
figures_dir <- here::here("outputs", "figures")
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Load an ASV table and a matching taxonomy table from the same analysis mode.
# This lets the script create both a full circular ASV tree and reduced genus/family views.
analysis_candidates <- list(
    list(
        label = "paired-end",
        asv = c(
            here::here("data", "processed", "dada2", "ASV_paired_end_table.rds"),
            here::here("outputs", "ASV_paired_end_table.rds")
        ),
        tax = c(
            here::here("data", "processed", "dada2", "taxonomy.rds"),
            here::here("outputs", "taxonomy.rds")
        )
    ),
    list(
        label = "single-end",
        asv = c(
            here::here("data", "processed", "dada2", "single_end", "ASV_table_SE.rds"),
            here::here("outputs", "single_end", "ASV_table_SE.rds")
        ),
        tax = c(
            here::here("data", "processed", "dada2", "single_end", "taxonomy_SE.rds"),
            here::here("outputs", "single_end", "taxonomy_SE.rds"),
            here::here("data", "processed", "dada2", "single_end", "taxonomy_SE_genuslevel.rds"),
            here::here("outputs", "single_end", "taxonomy_SE_genuslevel.rds")
        )
    )
)

# Minimum shared taxa required to accept a taxonomy table.
# If this value is between 0 and 1 it is interpreted as a fraction of total ASVs
# (e.g. 0.1 = 10% overlap). If >= 1 it is treated as an absolute count.
# Default: 10% overlap required for reasonably confident taxonomy collapse.
min_shared_taxa_for_taxonomy <- 0.1

standardize_taxonomy <- function(taxonomy_object) {
    taxonomy_matrix <- as.matrix(taxonomy_object)
    if (is.null(rownames(taxonomy_matrix))) {
        stop("The taxonomy table does not have row names, so it cannot be matched to ASVs.")
    }

    standard_names <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    lower_lookup <- setNames(standard_names, tolower(standard_names))
    matched_names <- unname(lower_lookup[tolower(colnames(taxonomy_matrix))])
    colnames(taxonomy_matrix) <- ifelse(is.na(matched_names), colnames(taxonomy_matrix), matched_names)
    taxonomy_matrix
}

load_analysis_inputs <- function(candidate_sets) {
    fallback_result <- NULL

    for (candidate_set in candidate_sets) {
        for (asv_path in candidate_set$asv[file.exists(candidate_set$asv)]) {
            asv_table <- as.matrix(readRDS(asv_path))
            if (is.null(colnames(asv_table))) {
                next
            }

            asv_sequences <- gsub("\\s+", "", colnames(asv_table))
            if (anyDuplicated(asv_sequences) > 0L) {
                stop("Duplicate ASV sequences were found in the ASV table column names.")
            }

            matching_taxonomy_path <- NA_character_
            taxonomy_table <- NULL
            for (tax_path in candidate_set$tax[file.exists(candidate_set$tax)]) {
                taxonomy_candidate <- tryCatch(readRDS(tax_path), error = function(e) NULL)
                if (is.null(taxonomy_candidate)) {
                    next
                }

                taxonomy_candidate <- standardize_taxonomy(taxonomy_candidate)
                shared_taxa <- intersect(asv_sequences, rownames(taxonomy_candidate))
                # Compute required overlap based on configured threshold
                required_shared <- if (min_shared_taxa_for_taxonomy > 0 && min_shared_taxa_for_taxonomy < 1) {
                    ceiling(length(asv_sequences) * min_shared_taxa_for_taxonomy)
                } else {
                    as.integer(min_shared_taxa_for_taxonomy)
                }
                if (length(shared_taxa) < required_shared) {
                    next
                }

                taxonomy_table <- matrix(
                    NA_character_,
                    nrow = length(asv_sequences),
                    ncol = ncol(taxonomy_candidate),
                    dimnames = list(asv_sequences, colnames(taxonomy_candidate))
                )
                matched_rows <- match(asv_sequences, rownames(taxonomy_candidate))
                present_rows <- !is.na(matched_rows)
                taxonomy_table[present_rows, ] <- taxonomy_candidate[matched_rows[present_rows], , drop = FALSE]
                matching_taxonomy_path <- tax_path
                break
            }

            if (is.null(fallback_result)) {
                fallback_result <- list(
                    mode = candidate_set$label,
                    asv_path = asv_path,
                    asv_table = asv_table,
                    asv_sequences = asv_sequences,
                    taxonomy_path = NA_character_,
                    taxonomy_table = NULL
                )
            }

            if (!is.null(taxonomy_table)) {
                return(list(
                    mode = candidate_set$label,
                    asv_path = asv_path,
                    asv_table = asv_table,
                    asv_sequences = asv_sequences,
                    taxonomy_path = matching_taxonomy_path,
                    taxonomy_table = taxonomy_table
                ))
            }
        }
    }

    if (!is.null(fallback_result)) {
        return(fallback_result)
    }

    stop("No ASV table was found. Looked for paired-end and single-end RDS outputs in the processed and outputs directories.")
}

analysis_input <- load_analysis_inputs(analysis_candidates)
asv_table <- analysis_input$asv_table
asv_sequences <- analysis_input$asv_sequences
taxonomy_table <- analysis_input$taxonomy_table
asv_table_path <- analysis_input$asv_path

# Convert ASV sequences to DNAStringSet for alignment and tree building.
asv_dna <- Biostrings::DNAStringSet(asv_sequences)
names(asv_dna) <- asv_sequences

cat("Loaded", nrow(asv_table), "samples and", ncol(asv_table), "ASVs from:", asv_table_path, " (analysis mode:", analysis_input$mode, ")\n")
if (!is.na(analysis_input$taxonomy_path)) {
    cat("Matched taxonomy table:", analysis_input$taxonomy_path, "\n")
} else {
    cat("No matching taxonomy table was found; only the full circular ASV tree will be created.\n")
}

cat("   Step: Aligning ASV sequences (DECIPHER)...\n")
# Align the ASV sequences before distance estimation and tree inference.
aligned_asvs <- DECIPHER::AlignSeqs(asv_dna, processors = 1L)
phang.alignment <- phangorn::phyDat(as(aligned_asvs, "matrix"), type = "DNA")

# Use a model-based distance matrix for the NJ starting tree.
# The JC69 model is a simple substitution model that assumes equal base frequencies and
# equal substitution rates.
asv_distance_matrix <- phangorn::dist.ml(phang.alignment, model = "JC69")
saveRDS(asv_distance_matrix, file.path(output_dir, "asv_distance_matrix.rds"))

cat("   Step: Building neighbor-joining starting tree...\n")
# Build a neighbor-joining tree as a starting point for maximum likelihood optimization.
nj_tree <- ape::nj(asv_distance_matrix)
nj_tree$edge.length[nj_tree$edge.length < 0] <- 0 # NJ can produce negative branch lengths,
# which are not biologically meaningful. Set them to zero.
saveRDS(nj_tree, file.path(output_dir, "asv_tree_nj.rds"))
ape::write.tree(nj_tree, file = file.path(output_dir, "asv_tree_nj.newick"))
# Save the NJ tree in Newick format for compatibility with other software.
cat("Neighbor-joining tree constructed and saved.\n")

cat("   Step: Optimising tree with maximum likelihood (GTR+G+I) — this may take a while...\n")
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
