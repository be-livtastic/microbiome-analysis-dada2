#!/usr/bin/env Rscript
# Parses ffq SRR-level JSON → data/processed/dada2/metadata.csv
# Handles flat SRR structure where sample/experiment are string accessions

library(jsonlite)

json_path <- file.path("outputs", "qc", "ffq_prjna762524.json")
if (!file.exists(json_path)) stop("JSON not found: ", json_path)

cat("Reading JSON...\n")
raw <- jsonlite::fromJSON(json_path, simplifyVector = FALSE)
cat("Top-level keys:", length(raw), "\n")

# Inspect first record to understand structure
first_key <- names(raw)[1]
first_run <- raw[[first_key]]
cat("First record key:", first_key, "\n")
cat("Fields in first record:", paste(names(first_run), collapse = ", "), "\n")
cat("Type of 'sample' field:", class(first_run[["sample"]]), "\n")
cat("Type of 'experiment' field:", class(first_run[["experiment"]]), "\n\n")

rows <- list()

for (srr in names(raw)) {
    run <- raw[[srr]]

    # sample and experiment may be plain strings (accession IDs only)
    # or nested lists -- handle both safely
    samp_field <- run[["sample"]]
    exp_field  <- run[["experiment"]]

    samp_acc   <- if (is.character(samp_field)) samp_field else
                  if (is.list(samp_field)) samp_field[["accession"]] %||% NA_character_ else NA_character_

    platform   <- if (is.list(exp_field)) exp_field[["platform"]] %||% NA_character_ else
                  run[["platform"]] %||% NA_character_

    instrument <- if (is.list(exp_field)) exp_field[["instrument"]] %||% NA_character_ else
                  run[["instrument"]] %||% NA_character_

    exp_acc    <- if (is.character(exp_field)) exp_field else
                  if (is.list(exp_field)) exp_field[["accession"]] %||% NA_character_ else NA_character_

    run_title  <- if (is.character(run[["title"]])) run[["title"]] else NA_character_

    row <- data.frame(
        RunAccession    = srr,
        SampleAccession = samp_acc,
        ExperimentAccession = exp_acc,
        RunTitle        = run_title,
        Platform        = platform,
        Instrument      = instrument,
        stringsAsFactors = FALSE
    )

    # Run-level attributes
    run_attrs <- run[["attributes"]]
    if (is.list(run_attrs) && length(run_attrs) > 0) {
        if (!is.null(names(run_attrs))) {
            for (k in names(run_attrs)) {
                v <- run_attrs[[k]]
                row[[make.names(k)]] <- if (is.null(v)) NA_character_ else as.character(v)
            }
        } else {
            for (item in run_attrs) {
                k <- make.names(if (!is.null(item[["tag"]])) item[["tag"]] else "unknown")
                v <- item[["value"]]
                row[[k]] <- if (is.null(v)) NA_character_ else as.character(v)
            }
        }
    }

    rownames(row) <- srr
    rows[[srr]]   <- row
}

if (length(rows) == 0) stop("No records extracted.")

# Bind rows, filling missing columns with NA
all_cols <- unique(unlist(lapply(rows, colnames)))
meta_df  <- do.call(rbind, lapply(rows, function(df) {
    missing     <- setdiff(all_cols, colnames(df))
    df[missing] <- NA_character_
    df[, all_cols, drop = FALSE]
}))

cat("Built metadata:", nrow(meta_df), "rows,", ncol(meta_df), "columns\n")
cat("Columns:", paste(colnames(meta_df), collapse = ", "), "\n\n")

# Validate against FASTQs
fastq_dir <- file.path("data", "raw", "fastq")
if (dir.exists(fastq_dir)) {
    fq_samples <- sub("_1\\.fastq\\.gz$", "",
                      list.files(fastq_dir, pattern = "_1\\.fastq\\.gz$"))
    if (length(fq_samples) > 0) {
        matched <- intersect(rownames(meta_df), fq_samples)
        cat(sprintf("Matched %d / %d FASTQs to metadata\n",
                    length(matched), length(fq_samples)))
        unmatched <- setdiff(fq_samples, rownames(meta_df))
        if (length(unmatched) > 0)
            cat("Unmatched FASTQs:", paste(unmatched, collapse = ", "), "\n")
    }
}

out_path <- file.path("data", "processed", "dada2", "metadata.csv")
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
write.csv(meta_df, file = out_path)
cat("\nSaved:", out_path, "\n")
cat("Preview:\n")
print(head(meta_df[, intersect(c("RunAccession","SampleAccession","RunTitle","Platform","Instrument"), colnames(meta_df))]))
