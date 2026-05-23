#!/usr/bin/env Rscript
# Alpha and Beta diversity analysis pipeline
# Uses outputs from DADA2/previous pipeline (ASV table, taxonomy, phylogeny, sample metadata)
# Provides logging via cat(), error handling, and commented explanations of the statistics used.

options(warn = 1)
set.seed(1)

cat("Starting alpha and beta diversity analysis pipeline...\n")

# Helper: choose the first existing file from candidates
choose_file <- function(cands) {
    existing <- cands[file.exists(cands)]
    if (length(existing) == 0) {
        return(NA_character_)
    }
    existing[1]
}

library_here <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(sprintf("Required package '%s' is not installed. Install it with install.packages('%s') or BiocManager::install('%s') if appropriate.", pkg, pkg, pkg))
    }
}

# Check key package
required_pkgs <- c("phyloseq", "vegan", "ggplot2", "dplyr", "tidyr", "ape")
for (p in required_pkgs) library_here(p)

library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ape)

# File candidates (match what earlier scripts may have written)
asv_candidates <- c(
    file.path("data", "processed", "dada2", "ASV_paired_end_table.rds"),
    file.path("outputs", "ASV_paired_end_table.rds"),
    file.path("data", "processed", "dada2", "single_end", "ASV_table_SE.rds"),
    file.path("outputs", "single_end", "ASV_table_SE.rds")
)

tax_candidates <- c(
    file.path("outputs", "single_end", "taxonomy_SE.rds"),
    file.path("outputs", "single_end", "taxonomy_SE_genuslevel.rds"),
    file.path("data", "processed", "dada2", "taxonomy.rds")
)

tree_candidates <- c(
    file.path("outputs", "phylogeny", "asv_tree_ml.rds"),
    file.path("outputs", "asv_tree_ml.rds"),
    file.path("phylogeny", "asv_tree_ml.rds")
)

meta_candidates <- c(
    file.path("data", "processed", "dada2", "sample_metadata.csv"),
    file.path("data", "processed", "dada2", "metadata.csv"),
    file.path("outputs", "sample_metadata.csv"),
    file.path("outputs", "metadata.csv")
)

asv_path <- choose_file(asv_candidates)
if (is.na(asv_path)) stop("No ASV table found. Looked in standard locations.")
cat("Loading ASV table from:", asv_path, "\n")

asv_obj <- tryCatch(
    {
        readRDS(asv_path)
    },
    error = function(e) stop("Failed to read ASV RDS: ", conditionMessage(e))
)

# Ensure ASV table is a matrix (samples x taxa)
asv_mat <- as.matrix(asv_obj)
if (is.null(colnames(asv_mat))) stop("ASV table has no column names (ASV sequences).")

# Load taxonomy if available (optional but recommended for summaries)
tax_path <- choose_file(tax_candidates)
tax_table_obj <- NULL
if (!is.na(tax_path)) {
    cat("Loading taxonomy from:", tax_path, "\n")
    tax_table_obj <- tryCatch(readRDS(tax_path), error = function(e) {
        warning("Failed to read taxonomy RDS: ", conditionMessage(e))
        NULL
    })
}

# Load phylogenetic tree if available (for Faith's PD and UniFrac)
tree_path <- choose_file(tree_candidates)
phy_tree_obj <- NULL
if (!is.na(tree_path)) {
    cat("Loading phylogenetic tree from:", tree_path, "\n")
    phy_tree_obj <- tryCatch(readRDS(tree_path), error = function(e) {
        warning("Failed to read tree RDS: ", conditionMessage(e))
        NULL
    })
}

# Load sample metadata if available; if not, proceed but many tests require grouping variables
meta_path <- choose_file(meta_candidates)
sample_metadata <- NULL
if (!is.na(meta_path)) {
    cat("Loading sample metadata from:", meta_path, "\n")
    sample_metadata <- tryCatch(read.csv(meta_path, stringsAsFactors = FALSE, row.names = 1), error = function(e) {
        warning("Failed to read metadata CSV: ", conditionMessage(e))
        NULL
    })
} else {
    cat("No metadata CSV found in standard locations; some tests (group comparisons) will be skipped.\n")
}

# Construct a phyloseq object when possible
ps <- NULL
try(
    {
        otu <- otu_table(asv_mat, taxa_are_rows = FALSE)
        components <- list(otu = otu)
        if (!is.null(tax_table_obj)) {
            # If tax is a matrix/data.frame convert to tax_table
            if (is.matrix(tax_table_obj) || is.data.frame(tax_table_obj)) {
                tt <- tax_table(as.matrix(tax_table_obj))
                components$tax <- tt
            }
        }
        if (!is.null(sample_metadata)) {
            sd <- sample_data(sample_metadata)
            components$sd <- sd
        }
        if (!is.null(phy_tree_obj)) {
            components$tree <- phy_tree_obj
        }
        ps <- do.call(phyloseq, components)
    },
    silent = TRUE
)

if (is.null(ps)) {
    cat("Could not fully build a phyloseq object; proceeding with matrix-based analyses where possible.\n")
}

# Create output directories
out_dir <- file.path("outputs", "alpha_beta")
fig_dir <- file.path("outputs", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

#############################
# ALPHA DIVERSITY
#############################
cat("Calculating alpha diversity metrics...\n")

# Explanation: Alpha diversity measures within-sample diversity. Common metrics:
# - Observed (richness): count of ASVs observed in a sample.
# - Shannon index: accounts for both richness and evenness; sensitive to rare taxa. # nolint
# - Simpson index: probability that two randomly selected individuals belong to the same species.
#   For interpretability we report the Simpson diversity as 1 - D (the complement), so larger values mean
#   greater diversity (this transformation is applied consistently below).
# - Faith's PD: phylogenetic diversity summing branch lengths represented in a sample (requires a phylogenetic tree).

alpha_df <- NULL
try(
    {
        if (!is.null(ps)) {
            alpha_df <- estimate_richness(ps, measures = c("Observed", "Shannon", "Simpson"))
            # estimate_richness returns Simpson's D; convert to Simpson diversity (1 - D) for consistency
            if ("Simpson" %in% colnames(alpha_df)) {
                alpha_df$Simpson <- 1 - alpha_df$Simpson
            }
        } else {
            # Fallback: compute from matrix
            obs <- rowSums(asv_mat > 0)
            shannon <- diversity(asv_mat, index = "shannon")
            simpson <- diversity(asv_mat, index = "simpson")
            alpha_df <- data.frame(Observed = obs, Shannon = shannon, Simpson = 1 - simpson, row.names = rownames(asv_mat))
        }
        # If phylogeny present, compute Faith's PD using picante::pd
        if (!is.null(phy_tree_obj) && requireNamespace("picante", quietly = TRUE)) {
            # Ensure tree tip labels match taxa names and subset tree properly using ape::keep.tip
            comm <- as.data.frame(asv_mat)
            if (!is.null(phy_tree_obj$tip.label)) {
                common_taxa <- intersect(colnames(comm), phy_tree_obj$tip.label)
                if (length(common_taxa) > 2) {
                    pruned_tree <- tryCatch(ape::keep.tip(phy_tree_obj, common_taxa), error = function(e) NULL)
                    if (!is.null(pruned_tree)) {
                        pd_res <- picante::pd(comm[, common_taxa, drop = FALSE], pruned_tree, include.root = TRUE)
                        alpha_df$Faith_PD <- pd_res$PD
                    } else {
                        warning("Failed to prune tree for Faith's PD calculation; skipping Faith's PD.")
                    }
                } else {
                    warning("Not enough overlapping taxa between ASV table and tree tips to compute Faith's PD.")
                }
            }
        }
        # Save alpha metrics
        write.csv(alpha_df, file = file.path(out_dir, "alpha_metrics.csv"), row.names = TRUE)
        cat("Alpha diversity metrics saved to:", file.path(out_dir, "alpha_metrics.csv"), "\n")

        # Compute Hill numbers (orders q = 0, 1, 2). Hill numbers are a family of diversity metrics:
        # - q = 0: species richness (counts all species equally)
        # - q = 1: exponential of Shannon entropy (gives weight proportional to abundance; sensitive to common and rare taxa)
        # - q = 2: inverse Simpson (gives more weight to common taxa)
        hill_number <- function(x, q) {
            if (sum(x) == 0) {
                return(NA_real_)
            }
            p <- x / sum(x)
            p <- p[p > 0]
            if (q == 1) {
                return(exp(-sum(p * log(p))))
            }
            return((sum(p^q))^(1 / (1 - q)))
        }

        q_vals <- c(0, 1, 2)
        hill_mat <- t(apply(asv_mat, 1, function(x) sapply(q_vals, function(q) hill_number(x, q))))
        colnames(hill_mat) <- paste0("Hill_q", q_vals)
        hill_df <- as.data.frame(hill_mat, stringsAsFactors = FALSE)
        write.csv(hill_df, file = file.path(out_dir, "hill_numbers_q0_q1_q2.csv"), row.names = TRUE)
        cat("Hill numbers (q=0,1,2) saved to:", file.path(out_dir, "hill_numbers_q0_q1_q2.csv"), "\n")

        # Also compute Hill numbers using vegan::renyi with hill=TRUE for cross-check (if available).
        if (requireNamespace("vegan", quietly = TRUE)) {
            renyi_hill <- tryCatch(
                vegan::renyi(asv_mat, scales = q_vals, hill = TRUE),
                error = function(e) {
                    warning("vegan::renyi failed: ", conditionMessage(e))
                    NULL
                }
            )
            if (!is.null(renyi_hill)) {
                # renyi returns samples as rows; ensure rownames
                write.csv(renyi_hill, file = file.path(out_dir, "renyi_hill_vegan.csv"), row.names = TRUE)
                cat("Vegan renyi (Hill numbers) saved to:", file.path(out_dir, "renyi_hill_vegan.csv"), "\n")
            }
        }
    },
    silent = FALSE
)

#############################
# RAREFACTION CURVES
#############################
cat("Plotting rarefaction curves...\n")
try(
    {
        asv_int <- round(asv_mat)
        rarefaction_path <- file.path(fig_dir, "rarefaction_curve.png")
        png(rarefaction_path, width = 10, height = 7, units = "in", res = 150)
        rarecurve(asv_int,
            step = 50,
            xlab = "Reads sampled",
            ylab = "ASVs observed",
            main = "Rarefaction Curves",
            label = FALSE
        )
        abline(v = min(rowSums(asv_int)), lty = 2, col = "red")
        dev.off()
        cat("Rarefaction curve saved to:", rarefaction_path, "\n")
    },
    silent = FALSE
)

#############################
# STATISTICAL TESTS ON ALPHA
#############################
cat("Performing alpha diversity group comparisons (if metadata available)...\n")

if (!is.null(sample_metadata)) {
    # Choose a grouping variable: prefer common names; otherwise pick the first categorical column.
    candidate_names <- c("Group", "group", "Treatment", "treatment", "SampleType", "sample_type")
    group_var <- intersect(candidate_names, colnames(sample_metadata))
    if (length(group_var) == 0) {
        # Helper to detect categorical-ish columns (not all unique and not constant)
        is_categorical_col <- function(vec, n) {
            u <- unique(vec)
            lu <- length(u)
            # treat as categorical if there are between 2 and n-1 unique values
            lu > 1 && lu < n
        }
        possible <- colnames(sample_metadata)[sapply(sample_metadata, is_categorical_col, n = nrow(sample_metadata))]
        group_var <- if (length(possible) > 0) possible[1] else NA_character_
    } else {
        group_var <- group_var[1]
    }

    if (!is.na(group_var)) {
        cat("Using grouping variable:", group_var, "for tests.\n")
        meta <- sample_metadata
        # Keep sample IDs as a column for safe merging
        meta$.sample_id <- rownames(meta)
        alpha_long <- data.frame(Sample = rownames(alpha_df), alpha_df, stringsAsFactors = FALSE)
        # Merge on the explicit Sample column and metadata row names
        merged <- merge(alpha_long, meta, by.x = "Sample", by.y = "row.names", all.x = TRUE)
        rownames(merged) <- merged$Sample
        # For each alpha metric run Kruskal-Wallis (non-parametric alternative to ANOVA):
        alpha_stats <- lapply(c("Observed", "Shannon", "Simpson", "Faith_PD"), function(metric) {
            if (!metric %in% colnames(merged)) {
                return(NULL)
            }
            # Kruskal-Wallis tests whether medians differ between groups (non-parametric)
            # Build a clean subset data.frame for the test to avoid inline subsetting in formulas
            test_df <- merged[!is.na(merged[[group_var]]), c(metric, group_var)]
            colnames(test_df) <- c("metric_val", "grouping")
            formula <- metric_val ~ grouping
            # Kruskal-Wallis tests whether medians differ between groups (non-parametric)
            kw <- kruskal.test(as.formula(formula), data = test_df)
            # Pairwise Wilcoxon for post-hoc (with BH p-value correction)
            pairwise <- pairwise.wilcox.test(test_df$metric_val, test_df$grouping, p.adjust.method = "BH")
            list(metric = metric, kruskal = kw, pairwise = pairwise)
        })
        saveRDS(alpha_stats, file = file.path(out_dir, "alpha_stats.rds"))
        cat("Alpha diversity statistics saved to:", file.path(out_dir, "alpha_stats.rds"), "\n")
    } else {
        cat("No suitable grouping variable found in metadata; skipping group tests for alpha diversity.\n")
    }
} else {
    cat("No metadata available; skipping alpha group statistics.\n")
}

#############################
# BETA DIVERSITY
#############################
cat("Calculating beta diversity (dissimilarity matrices) and running PERMANOVA...\n")

beta_results <- list()
try(
    {
        # Use relative abundances for Bray-Curtis
        comm <- asv_mat
        rel_comm <- decostand(comm, method = "total")
        # Bray-Curtis distance (community composition)
        bray <- vegdist(rel_comm, method = "bray")
        beta_results$bray <- bray
        saveRDS(bray, file = file.path(out_dir, "bray_distance.rds"))
        cat("Bray-Curtis distance saved to:", file.path(out_dir, "bray_distance.rds"), "\n")

        # UniFrac distances if phylogeny available and phyloseq built
        if (!is.null(ps) && !is.null(phy_tree_obj)) {
            cat("Computing UniFrac distances (requires phyloseq).\n")
            wunifrac <- phyloseq::distance(ps, method = "unifrac", weighted = TRUE)
            runifrac <- phyloseq::distance(ps, method = "unifrac", weighted = FALSE)
            beta_results$wunifrac <- wunifrac
            beta_results$runifrac <- runifrac
            saveRDS(wunifrac, file = file.path(out_dir, "unifrac_weighted.rds"))
            saveRDS(runifrac, file = file.path(out_dir, "unifrac_unweighted.rds"))
        }
    },
    silent = FALSE
)

# Runs PERMANOVA (adonis2) + betadisper as a paired unit for one distance matrix.
# betadisper tests homogeneity of dispersions; it must accompany every PERMANOVA
# because adonis2 is sensitive to differences in dispersion as well as centroids.
run_permanova_betadisper <- function(dist_obj, meta_df, group_col, label, out_dir) {
    if (is.null(dist_obj)) {
        cat("Skipping PERMANOVA/betadisper for", label, ": distance object is NULL.\n")
        return(invisible(NULL))
    }
    common_samples <- intersect(rownames(meta_df), attr(dist_obj, "Labels"))
    if (length(common_samples) < 3) {
        warning(sprintf("Not enough overlapping samples for PERMANOVA on %s (%d found).", label, length(common_samples)))
        return(invisible(NULL))
    }
    sub_meta <- meta_df[common_samples, , drop = FALSE]
    sub_meta$grouping <- sub_meta[[group_col]]
    ad <- adonis2(dist_obj ~ grouping, data = sub_meta, permutations = 999)
    saveRDS(ad, file = file.path(out_dir, paste0("permanova_", label, ".rds")))
    cat("PERMANOVA results saved to:", file.path(out_dir, paste0("permanova_", label, ".rds")), "\n")
    bd <- betadisper(dist_obj, sub_meta$grouping)
    bd_anova <- anova(bd)
    saveRDS(list(betadisper = bd, anova = bd_anova), file = file.path(out_dir, paste0("betadisper_", label, ".rds")))
    cat("Betadisper results saved to:", file.path(out_dir, paste0("betadisper_", label, ".rds")), "\n")
    invisible(list(permanova = ad, betadisper = bd, betadisper_anova = bd_anova))
}

#############################
# PERMANOVA and dispersion (betadisper)
#############################
if (!is.null(sample_metadata)) {
    # reuse group_var from alpha section
    if (!exists("group_var") || is.na(group_var)) {
        cat("No grouping variable found earlier; skipping PERMANOVA.\n")
    } else {
        cat("Running PERMANOVA (adonis2) using grouping variable:", group_var, "\n")
        meta <- sample_metadata
        run_permanova_betadisper(beta_results$bray, meta, group_var, "bray", out_dir)
        if (!is.null(beta_results$wunifrac)) {
            run_permanova_betadisper(beta_results$wunifrac, meta, group_var, "wunifrac", out_dir)
        }
        if (!is.null(beta_results$runifrac)) {
            run_permanova_betadisper(beta_results$runifrac, meta, group_var, "runifrac", out_dir)
        }
    }
} else {
    cat("No metadata available; skipping PERMANOVA and betadisper.\n")
}

#############################
# ORDINATION AND PLOTTING
#############################
cat("Creating ordination plots (NMDS for Bray-Curtis; static PNGs)...\n")

try(
    {
        # NMDS on Bray
        nmds <- metaMDS(beta_results$bray, k = 2, trymax = 100)
        nmds_scores <- as.data.frame(scores(nmds))
        nmds_scores$Sample <- rownames(nmds_scores)
        # Merge sample metadata if available
        if (!is.null(sample_metadata)) {
            nmds_scores <- merge(nmds_scores, sample_metadata, by.x = "Sample", by.y = "row.names", all.x = TRUE)
            rownames(nmds_scores) <- nmds_scores$Sample
        }

        # Build plot: include color aesthetic only if group_var is valid and present
        if (!is.null(sample_metadata) && !is.na(group_var) && group_var %in% colnames(nmds_scores)) {
            p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = .data[[group_var]])) +
                geom_point(size = 3, alpha = 0.8) +
                theme_minimal() +
                labs(title = "NMDS (Bray-Curtis)", color = group_var)
        } else {
            p <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2)) +
                geom_point(size = 3, alpha = 0.8) +
                theme_minimal() +
                labs(title = "NMDS (Bray-Curtis)")
        }
        ggsave(filename = file.path(fig_dir, "nmds_bray.png"), plot = p, width = 7, height = 6)
        cat("Saved NMDS plot to:", file.path(fig_dir, "nmds_bray.png"), "\n")
    },
    silent = FALSE
)

#############################
# SUMMARY
#############################
cat("Alpha/Beta pipeline finished. Outputs are in:", out_dir, "and figures in:", fig_dir, "\n")

cat("Notes and interpretation help:\n")
cat("- Alpha diversity tests: Kruskal-Wallis tests whether the distribution (median) of a diversity metric differs across groups. If significant, pairwise Wilcoxon tests provide post-hoc comparisons with BH correction.\n")
cat("- Beta diversity: Bray-Curtis quantifies compositional differences between samples; UniFrac incorporates phylogenetic relatedness.\n")
cat("- PERMANOVA (adonis2): tests whether group centroids differ in multivariate space; sensitive to differences in dispersion. betadisper is always run alongside to check homogeneity of dispersion (permanova_<label>.rds / betadisper_<label>.rds).\n")
cat("- NMDS: non-metric multidimensional scaling for visualizing sample relationships based on distance matrices.\n")

### END
