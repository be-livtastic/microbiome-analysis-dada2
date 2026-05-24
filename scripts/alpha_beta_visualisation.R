#!/usr/bin/env Rscript
options(warn = 1)
set.seed(1)

cat("Starting alpha and beta diversity visualisation workflow...\n")

choose_file <- function(cands) {
    existing <- cands[file.exists(cands)]
    if (length(existing) == 0) return(NA_character_)
    existing[1]
}

check_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE))
        stop(sprintf("Required package '%s' is not installed.", pkg))
}

for (pkg in c("phyloseq","ComplexHeatmap","circlize","ggplot2",
              "dplyr","tidyr","vegan","ape","RColorBrewer","scales")) check_pkg(pkg)

library(phyloseq); library(ComplexHeatmap); library(circlize)
library(ggplot2);  library(dplyr);          library(tidyr)
library(vegan);    library(ape);            library(RColorBrewer); library(scales)

out_fig_dir   <- file.path("outputs","figures")
out_alpha_dir <- file.path("outputs","alpha_beta")
dir.create(out_fig_dir,   recursive=TRUE, showWarnings=FALSE)
dir.create(out_alpha_dir, recursive=TRUE, showWarnings=FALSE)

asv_path <- choose_file(c(
    file.path("data","processed","dada2","ASV_paired_end_table.rds"),
    file.path("data","processed","dada2","ASV_table.rds"),
    file.path("outputs","ASV_paired_end_table.rds"),
    file.path("outputs","ASV_table.rds")))
tax_path <- choose_file(c(
    file.path("data","processed","dada2","taxonomy.rds"),
    file.path("outputs","taxonomy.rds")))
meta_path <- choose_file(c(
    file.path("data","processed","dada2","metadata.csv")))
tree_path <- choose_file(c(
    file.path("outputs","phylogeny","asv_tree_ml.rds"),
    file.path("outputs","asv_tree_ml.rds")))
alpha_metrics_path <- choose_file(c(
    file.path(out_alpha_dir,"alpha_metrics.csv"),
    file.path("outputs","alpha_beta","alpha_metrics.csv")))

if (is.na(asv_path)) stop("No ASV table found.")

cat("Loading ASV table from:", asv_path, "\n")
asv_mat <- as.matrix(readRDS(asv_path))

taxonomy_obj <- NULL
if (!is.na(tax_path)) {
    cat("Loading taxonomy from:", tax_path, "\n")
    taxonomy_obj <- tryCatch(readRDS(tax_path), error=function(e) NULL)
}

sample_metadata <- NULL
if (!is.na(meta_path)) {
    cat("Loading metadata from:", meta_path, "\n")
    sample_metadata <- tryCatch(
        read.csv(meta_path, stringsAsFactors=FALSE, row.names=1),
        error=function(e) NULL)
}

phy_tree_obj <- NULL
if (!is.na(tree_path)) {
    cat("Loading tree from:", tree_path, "\n")
    phy_tree_obj <- tryCatch(readRDS(tree_path), error=function(e) NULL)
}

# Align taxa
if (!is.null(taxonomy_obj)) {
    tax_mat     <- as.matrix(taxonomy_obj)
    shared_taxa <- intersect(colnames(asv_mat), rownames(tax_mat))
    if (length(shared_taxa) > 0) {
        asv_mat      <- asv_mat[, shared_taxa, drop=FALSE]
        taxonomy_obj <- tax_mat[shared_taxa, , drop=FALSE]
    } else {
        warning("Taxonomy does not match ASV IDs; skipping taxonomy.")
        taxonomy_obj <- NULL
    }
}

# Align samples
if (!is.null(sample_metadata)) {
    shared_samp     <- intersect(rownames(asv_mat), rownames(sample_metadata))
    cat("Shared samples:", length(shared_samp), "\n")
    asv_mat         <- asv_mat[shared_samp, , drop=FALSE]
    sample_metadata <- sample_metadata[shared_samp, , drop=FALSE]
}

# Build phyloseq
otu        <- otu_table(asv_mat, taxa_are_rows=FALSE)
components <- list(otu=otu)
if (!is.null(taxonomy_obj))  components$tax <- tax_table(as.matrix(taxonomy_obj))
if (!is.null(sample_metadata)) components$sd <- sample_data(sample_metadata)

ps     <- do.call(phyloseq, components)
ps     <- prune_samples(sample_sums(ps) > 0, ps)
ps_rel <- transform_sample_counts(ps, function(x) { s<-sum(x); if(s==0) x else x/s })
ps_rel <- prune_samples(sample_sums(ps_rel) > 0, ps_rel)

cat("Phyloseq:", nsamples(ps), "samples,", ntaxa(ps), "taxa\n")
cat("Sample data attached:", !is.null(sample_data(ps, errorIfNULL=FALSE)), "\n")

# Group variable
pick_group_var <- function(meta) {
    if (is.null(meta) || ncol(meta)==0) return(NA_character_)
    cands <- c("Group","group","host_sex","Treatment","treatment","SampleType","Condition","condition")
    direct <- intersect(cands, colnames(meta))
    if (length(direct) > 0) return(direct[1])
    is_cat <- function(v,n) { u<-unique(v[!is.na(v)]); length(u)>1 && length(u)<n }
    possible <- colnames(meta)[vapply(meta, is_cat, logical(1), n=nrow(meta))]
    if (length(possible) > 0) possible[1] else NA_character_
}
group_var <- pick_group_var(sample_metadata)
cat(if (!is.na(group_var)) paste("Using grouping variable:", group_var) else
    "No grouping variable found; sample order only.", "\n")

# Sample order
sample_depths   <- sample_sums(ps)
ordered_samples <- if (!is.null(sample_metadata) && !is.na(group_var) &&
                        group_var %in% colnames(sample_metadata)) {
    gk <- as.character(sample_metadata[sample_names(ps_rel), group_var])
    gk[is.na(gk)] <- "zzz"
    sample_names(ps_rel)[order(gk, -sample_depths[sample_names(ps_rel)], sample_names(ps_rel))]
} else {
    sample_names(ps_rel)[order(-sample_depths[sample_names(ps_rel)], sample_names(ps_rel))]
}

# Alpha diversity
alpha_metrics <- NULL
if (!is.na(alpha_metrics_path)) {
    cat("Loading alpha metrics from:", alpha_metrics_path, "\n")
    alpha_metrics <- tryCatch(read.csv(alpha_metrics_path, row.names=1, check.names=FALSE), error=function(e) NULL)
}
if (is.null(alpha_metrics)) {
    cat("Computing alpha metrics.\n")
    alpha_metrics <- estimate_richness(ps, measures=c("Observed","Shannon","Simpson"))
    if ("Simpson" %in% colnames(alpha_metrics)) alpha_metrics$Simpson <- 1 - alpha_metrics$Simpson
}
alpha_metrics$Sample <- rownames(alpha_metrics)
if (!is.null(sample_metadata))
    alpha_metrics <- merge(alpha_metrics, sample_metadata, by.x="Sample", by.y="row.names", all.x=TRUE)

alpha_long <- alpha_metrics %>%
    select(Sample, any_of(c("Observed","Shannon","Simpson","Faith_PD")), everything()) %>%
    pivot_longer(cols=any_of(c("Observed","Shannon","Simpson","Faith_PD")),
                 names_to="Metric", values_to="Value")

if (nrow(alpha_long) > 0) {
    alpha_long$Sample <- factor(alpha_long$Sample, levels=ordered_samples)
    p_alpha <- ggplot(alpha_long, aes(x=Sample, y=Value)) +
        geom_boxplot(outlier.shape=NA, fill="#d9e8fb", color="#4c78a8", width=0.7) +
        facet_wrap(~Metric, scales="free_y", ncol=2) +
        labs(title="Alpha diversity summary", x="Sample", y="Metric value") +
        theme_minimal(base_size=12) +
        theme(axis.text.x=element_text(angle=45, hjust=1),
              legend.position=if (!is.na(group_var) && group_var %in% colnames(alpha_long)) "top" else "none")
    if (!is.na(group_var) && group_var %in% colnames(alpha_long)) {
        p_alpha <- p_alpha +
            geom_point(aes(color=.data[[group_var]]), size=1.8, alpha=0.75,
                       position=position_jitter(width=0.12, height=0)) + labs(color=group_var)
    } else {
        p_alpha <- p_alpha + geom_point(size=1.8, alpha=0.75, position=position_jitter(width=0.12))
    }
    ggsave(file.path(out_fig_dir,"alpha_diversity_summary_boxplots.png"), p_alpha, width=12, height=7, dpi=300)
    cat("Saved alpha_diversity_summary_boxplots.png\n")
}

# Taxonomic barplots
top_taxon_palette <- function(vals) {
    n <- length(vals); if (n<=0) return(character())
    setNames(c(grDevices::hcl.colors(max(1,n-1), "Dark 3"), "#d9d9d9"), c(vals[seq_len(max(1,n-1))], "Other"))
}

prepare_barplot <- function(ps_obj, rank, top_n=10) {
    if (is.null(tax_table(ps_obj, errorIfNULL=FALSE))) return(NULL)
    bar_ps <- tryCatch(
        phyloseq::tax_glom(ps_obj, taxrank=rank, NArm=FALSE),
        error=function(e) { warning(sprintf("tax_glom failed: %s", conditionMessage(e))); NULL })
    if (is.null(bar_ps)) return(NULL)
    if (is.null(access(bar_ps,"sam_data",errorIfNULL=FALSE)) &&
        !is.null(access(ps_obj,"sam_data",errorIfNULL=FALSE)))
        sample_data(bar_ps) <- sample_data(ps_obj)

    tax_df <- as.data.frame(tax_table(bar_ps), stringsAsFactors=FALSE)
    if (!rank %in% colnames(tax_df)) { warning(sprintf("Rank '%s' missing.", rank)); return(NULL) }

    abund    <- taxa_sums(bar_ps)
    top_taxa <- names(sort(abund, decreasing=TRUE))[seq_len(min(top_n, length(abund)))]
    tax_df[[rank]][is.na(tax_df[[rank]]) | tax_df[[rank]]==""] <- "Unassigned"
    tax_df$PlotTaxon <- factor(ifelse(tax_df[[rank]] %in% top_taxa, tax_df[[rank]], "Other"),
                                levels=c(top_taxa,"Other"))
    tax_table(bar_ps) <- tax_table(as.matrix(tax_df))

    meta_df <- as.data.frame(sample_data(ps_obj))
    sample_order <- if (!is.na(group_var) && group_var %in% colnames(meta_df)) {
        gk <- as.character(meta_df[[group_var]]); gk[is.na(gk)] <- "zzz"
        rownames(meta_df)[order(gk, -sample_sums(bar_ps)[rownames(meta_df)], rownames(meta_df))]
    } else {
        rownames(meta_df)[order(-sample_sums(bar_ps)[rownames(meta_df)], rownames(meta_df))]
    }
    sdf <- as.data.frame(sample_data(bar_ps))
    sdf$Sample <- factor(rownames(sdf), levels=sample_order)
    sample_data(bar_ps) <- sample_data(sdf)

    p <- plot_bar(bar_ps, x="Sample", fill="PlotTaxon") +
        labs(title=paste0("Relative abundance — ", rank), x="Sample", y="Relative abundance") +
        theme_bw(base_size=11) +
        theme(axis.text.x=element_text(angle=60, hjust=1, vjust=1), legend.position="top") +
        scale_fill_manual(values=top_taxon_palette(levels(droplevels(tax_df$PlotTaxon))))
    if (!is.na(group_var) && group_var %in% sample_variables(bar_ps))
        p <- p + facet_grid(reformulate(group_var), scales="free_x", space="free_x")

    list(phyloseq=bar_ps, plot=p, top_taxa=top_taxa)
}

bar_phylum <- prepare_barplot(ps_rel, "Phylum", 10)
if (!is.null(bar_phylum)) {
    ggsave(file.path(out_fig_dir,"taxonomic_barplot_phylum_relative_abundance.png"),
           bar_phylum$plot, width=14, height=7, dpi=300)
    cat("Saved taxonomic_barplot_phylum_relative_abundance.png\n")
}
bar_genus <- prepare_barplot(ps_rel, "Genus", 10)
if (!is.null(bar_genus)) {
    ggsave(file.path(out_fig_dir,"taxonomic_barplot_genus_relative_abundance.png"),
           bar_genus$plot, width=14, height=7, dpi=300)
    cat("Saved taxonomic_barplot_genus_relative_abundance.png\n")
}

# PCoA
make_pcoa_plot <- function(dist_obj, ps_obj, title, file_name, grp=NA_character_) {
    pcoa_res <- ape::pcoa(as.dist(dist_obj))
    scores   <- as.data.frame(pcoa_res$vectors[,1:2,drop=FALSE])
    colnames(scores) <- c("Axis1","Axis2"); scores$Sample <- rownames(scores)
    if (!is.na(grp) && grp %in% sample_variables(ps_obj)) {
        meta_df <- as.data.frame(sample_data(ps_obj))
        scores  <- merge(scores, meta_df, by.x="Sample", by.y="row.names", all.x=TRUE)
    }
    ax1 <- round(100*pcoa_res$values$Relative_eig[1],1)
    ax2 <- round(100*pcoa_res$values$Relative_eig[2],1)
    p <- ggplot(scores, aes(x=Axis1, y=Axis2)) +
        geom_hline(yintercept=0, color="grey85", linewidth=0.3) +
        geom_vline(xintercept=0, color="grey85", linewidth=0.3) +
        labs(title=title, x=paste0("PCoA1 (",ax1,"%)"), y=paste0("PCoA2 (",ax2,"%)")) +
        theme_minimal(base_size=12) + coord_equal() +
        theme(legend.position=if (!is.na(grp) && grp %in% colnames(scores)) "top" else "none")
    if (!is.na(grp) && grp %in% colnames(scores)) {
        p <- p + geom_point(aes(color=.data[[grp]]), size=3, alpha=0.85) + labs(color=grp)
    } else { p <- p + geom_point(size=3, alpha=0.85) }
    ggsave(file.path(out_fig_dir,file_name), p, width=7.5, height=6.5, dpi=300)
}

bray_dist <- phyloseq::distance(ps_rel, method="bray")
make_pcoa_plot(bray_dist, ps_rel, "Bray-Curtis PCoA", "beta_pcoa_bray_curtis.png", group_var)
cat("Saved beta_pcoa_bray_curtis.png\n")

ps_tree_cand <- phy_tree(ps_rel, errorIfNULL=FALSE)
if (!is.null(phy_tree_obj) && !is.null(ps_tree_cand)) {
    shared_tips <- intersect(ps_tree_cand$tip.label, taxa_names(ps_rel))
    if (length(shared_tips) > 2) {
        ps_rel_tree <- prune_taxa(shared_tips, ps_rel)
        if (ntaxa(ps_rel_tree) > 2) {
            wunifrac_dist <- phyloseq::distance(ps_rel_tree, method="wunifrac")
            make_pcoa_plot(wunifrac_dist, ps_rel_tree, "Weighted UniFrac PCoA",
                           "beta_pcoa_weighted_unifrac.png", group_var)
            cat("Saved beta_pcoa_weighted_unifrac.png\n")
        }
    }
}

# NMDS
nmds_res    <- vegan::metaMDS(bray_dist, k=2, trymax=100, autotransform=FALSE)
nmds_scores <- as.data.frame(scores(nmds_res))
nmds_scores$Sample <- rownames(nmds_scores)
if (!is.null(sample_metadata))
    nmds_scores <- merge(nmds_scores, sample_metadata, by.x="Sample", by.y="row.names", all.x=TRUE)

nmds_plot <- ggplot(nmds_scores, aes(x=NMDS1, y=NMDS2)) +
    geom_hline(yintercept=0, color="grey85", linewidth=0.3) +
    geom_vline(xintercept=0, color="grey85", linewidth=0.3) +
    labs(title="Bray-Curtis NMDS", subtitle=paste0("Stress = ",round(nmds_res$stress,3)),
         x="NMDS1", y="NMDS2") +
    theme_minimal(base_size=12) + coord_equal() +
    theme(legend.position=if (!is.na(group_var) && group_var %in% colnames(nmds_scores)) "top" else "none")
if (!is.na(group_var) && group_var %in% colnames(nmds_scores)) {
    nmds_plot <- nmds_plot +
        geom_point(aes(color=.data[[group_var]]), size=3, alpha=0.85) + labs(color=group_var)
} else { nmds_plot <- nmds_plot + geom_point(size=3, alpha=0.85) }
ggsave(file.path(out_fig_dir,"beta_nmds_bray_curtis.png"), nmds_plot, width=7.5, height=6.5, dpi=300)
cat("Saved beta_nmds_bray_curtis.png\n")

# Heatmap
make_heatmap <- function(ps_obj, rank="Genus", top_n=30, file_name="asv_heatmap.png") {
    n_taxa <- min(top_n, ntaxa(ps_obj))
    ps_top <- prune_taxa(names(sort(taxa_sums(ps_obj), decreasing=TRUE))[seq_len(n_taxa)], ps_obj)
    ps_top <- transform_sample_counts(ps_top, function(x) { s<-sum(x); if(s==0) x else x/s })

    mat <- as.matrix(otu_table(ps_top))
    if (!taxa_are_rows(ps_top)) mat <- t(mat)
    cat("Heatmap matrix (taxa x samples):", nrow(mat), "x", ncol(mat), "\n")
    mat <- log10(mat + 1e-6)

    # Row labels -- safe length-checked assignment
    tax_df <- NULL
    if (!is.null(tax_table(ps_top, errorIfNULL=FALSE))) {
        tax_df <- as.data.frame(tax_table(ps_top), stringsAsFactors=FALSE)
        tax_df <- tax_df[rownames(mat), , drop=FALSE]
        label_rank <- if (rank %in% colnames(tax_df)) {
    rank
} else if ("Genus" %in% colnames(tax_df)) {
    "Genus"
} else {
    NULL
}
        if (!is.null(label_rank)) {
            labs <- tax_df[[label_rank]]
            labs[is.na(labs) | labs==""] <- rownames(tax_df)[is.na(labs) | labs==""]
            labs <- make.unique(as.character(labs))
            if (length(labs) == nrow(mat)) {
                rownames(mat)    <- labs
                rownames(tax_df) <- labs
            } else {
                cat("Label length mismatch; using ASV IDs\n")
            }
        }
    }

    # Column annotation
  ann_col <- NULL
    sd_obj  <- sample_data(ps_top, errorIfNULL=FALSE)
    if (!is.null(sd_obj) && !is.na(group_var) && group_var %in% sample_variables(ps_top)) {
        ann_df <- as.data.frame(sd_obj)
        # Subset to exactly the columns of mat in correct order
        ann_df <- ann_df[colnames(mat), , drop=FALSE]
        grp_vec <- as.factor(ann_df[[group_var]])
        names(grp_vec) <- colnames(mat)
        ann_col <- HeatmapAnnotation(Group = grp_vec)
    }

    # Row annotation
    ann_row <- NULL
    if (!is.null(tax_df) && all(c("Phylum","Genus") %in% colnames(tax_df))) {
        rd <- tax_df[rownames(mat), c("Phylum","Genus"), drop=FALSE]
        ann_row <- rowAnnotation(Phylum=rd$Phylum, Genus=rd$Genus)
    }

    mat_vec    <- as.vector(mat)
    col_breaks <- c(min(mat_vec,na.rm=TRUE), median(mat_vec,na.rm=TRUE), max(mat_vec,na.rm=TRUE))

    ht <- Heatmap(mat,
        name                      = "log10(rel. ab.)",
        col                       = colorRamp2(col_breaks, c("#08306b","#f7f7f7","#67000d")),
        cluster_rows              = nrow(mat) > 1,
        cluster_columns           = ncol(mat) > 1,
        clustering_method_rows    = "ward.D2",
        clustering_method_columns = "ward.D2",
        show_row_names            = TRUE,
        show_column_names         = TRUE,
        row_names_gp              = grid::gpar(fontsize=8),
        column_names_gp           = grid::gpar(fontsize=8),
        column_names_rot          = 45,
        top_annotation            = ann_col,
        left_annotation           = ann_row,
        heatmap_legend_param      = list(title="log10\nrelative abundance"))

    png(file.path(out_fig_dir, file_name), width=1800, height=1400, res=200)
    draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
    dev.off()
    cat("Saved", file_name, "\n")
}

make_heatmap(ps_rel, rank="Genus", top_n=min(30, ntaxa(ps_rel)),
             file_name="top_asv_heatmap_complexheatmap_ward_d2.png")

cat("\nVisualisation workflow complete.\n")