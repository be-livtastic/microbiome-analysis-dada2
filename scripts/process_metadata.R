#!/usr/bin/env Rscript
# Parses ffq JSON for PRJNA762524 → data/processed/dada2/metadata.csv
# Row names = SRR accession numbers, matching FASTQ filename convention:
#   SRR15862403_1.fastq.gz  →  sample name SRR15862403
# Run via: Rscript scripts/process_metadata.R
# Or called automatically by: bash scripts/fetch_metadata.sh

`%||%` <- function(a, b) if (!is.null(a)) a else b

if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite required: install.packages('jsonlite')")
}
library(jsonlite)

json_path <- file.path("outputs", "qc", "ffq_prjna762524.json")
if (!file.exists(json_path)) {
    stop("ffq JSON not found at: ", json_path, "\nRun scripts/fetch_metadata.sh first.")
}

cat("Reading ffq JSON from:", json_path, "\n")
raw <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)

if (!is.list(raw) || length(raw) == 0) {
    stop("ffq JSON is empty or malformed. Re-run scripts/fetch_metadata.sh.")
}

cat("Found", length(raw), "run records.\n")

# Each top-level element is one sequencing run.
# Extract accession + platform/instrument + all sample attributes.
meta_list <- lapply(raw, function(run) {
    acc <- run[["accession"]]
    if (is.null(acc) || nchar(acc) == 0) return(NULL)

    row <- data.frame(
        stringsAsFactors = FALSE,
        SampleAccession = run[["sample"]][["accession"]] %||% NA_character_,
        Platform        = run[["platform"]]              %||% NA_character_,
        Instrument      = run[["instrument"]]            %||% NA_character_
    )

    # Flatten sample attributes (study-specific key-value pairs)
    attrs <- run[["sample"]][["attributes"]]
    if (is.list(attrs) && length(attrs) > 0) {
        # attrs may be a named list or a list of {tag, value} objects
        if (!is.null(names(attrs))) {
            attr_df <- as.data.frame(lapply(attrs, function(v) if (is.null(v)) NA_character_ else as.character(v)),
                stringsAsFactors = FALSE)
        } else {
            # List of {tag, value} pairs
            tags   <- vapply(attrs, function(x) x[["tag"]]   %||% NA_character_, character(1))
            values <- vapply(attrs, function(x) x[["value"]] %||% NA_character_, character(1))
            attr_df <- as.data.frame(
                setNames(as.list(values), make.names(tags, unique = TRUE)),
                stringsAsFactors = FALSE
            )
        }
        row <- cbind(row, attr_df)
    }

    rownames(row) <- acc
    row
})

meta_list <- Filter(Negate(is.null), meta_list)

# Bind rows, filling any missing columns with NA
all_cols <- unique(unlist(lapply(meta_list, colnames)))
meta_df <- do.call(rbind, lapply(meta_list, function(df) {
    missing <- setdiff(all_cols, colnames(df))
    df[missing] <- NA_character_
    df[, all_cols, drop = FALSE]
}))

cat("Built metadata data frame:", nrow(meta_df), "samples,", ncol(meta_df), "columns.\n")
cat("Columns:", paste(colnames(meta_df), collapse = ", "), "\n\n")

# Validate Sample IDs against FASTQ filenames (if raw data already downloaded)
fastq_dir <- file.path("data", "raw", "fastq")
if (dir.exists(fastq_dir)) {
    fq_files   <- list.files(fastq_dir, pattern = "_1\\.fastq\\.gz$", full.names = FALSE)
    fq_samples <- sub("_1\\.fastq\\.gz$", "", fq_files)
    if (length(fq_samples) > 0) {
        in_both    <- intersect(rownames(meta_df), fq_samples)
        only_fq    <- setdiff(fq_samples, rownames(meta_df))
        only_meta  <- setdiff(rownames(meta_df), fq_samples)
        cat(sprintf("Sample ID validation: %d matched, %d FASTQ-only, %d metadata-only.\n",
            length(in_both), length(only_fq), length(only_meta)))
        if (length(only_fq) > 0) {
            warning("FASTQ files with no metadata entry: ", paste(only_fq, collapse = ", "))
        }
        if (length(only_meta) > 0) {
            cat(sprintf("  (%d metadata samples have no FASTQ yet — normal if data not fully downloaded)\n",
                length(only_meta)))
        }
    } else {
        cat("data/raw/fastq/ exists but no *_1.fastq.gz files found. Download data first.\n")
    }
} else {
    cat("data/raw/fastq/ not found; skipping FASTQ validation.\n")
    cat("Re-run this script after downloading data to confirm Sample ID alignment.\n")
}

out_path <- file.path("data", "processed", "dada2", "metadata.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(meta_df, file = out_path)
cat("\nMetadata CSV saved to:", out_path, "\n")
cat("Row names are SRR accession numbers and will match sample names produced by dada2_pipeline.R.\n")
