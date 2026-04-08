# ---- LOAD REQUIRED PACKAGES ----

# Set random seed for reproducibility
# Any random operations (e.g., pseudo-pooling) will give same results each run
set.seed(1)

# ---- DEFINE PROJECT PATHS ----
# Using here::here() makes paths work regardless of working directory
# This means the script runs on any computer without modification

project_dir <- here::here() # Root directory of the R project

# Where raw FASTQ files are stored
raw_dir <- here::here("data", "raw", "fastq")

# Where filtered/processed files will be saved
filtered_dir <- here::here("data", "processed", "dada2")

# SILVA reference database for taxonomic classification
# This file assigns taxonomy (Kingdom → Species) to our sequences
reference_fasta <- here::here(
  "data", "external", "reference",
  "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
)


# Create the output directory if it doesn't exist; this is where all DADA2 results will be saved
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

fastq_files <- list.files(raw_dir)
print(fastq_files)

forward_reads <- sort(list.files(raw_dir, pattern = "_1\\.fastq$", full.names = TRUE))
reverse_reads <- sort(list.files(raw_dir, pattern = "_2\\.fastq$", full.names = TRUE))

# Basic checks to confirm files are found and paired correctly before processing
if (length(forward_reads) == 0L || length(reverse_reads) == 0L) {
  stop("No paired FASTQ files were found in: ", raw_dir)
}
if (length(forward_reads) != length(reverse_reads)) {
  stop("The number of forward and reverse FASTQ files do not match.")
}
if (length(forward_reads) != length(reverse_reads)) {
  stop("Forward and reverse read counts do not match.")
}

# Extract sample names by removing the read suffix and file extension; this will be used for naming outputs
sample_names <- sub("_1\\.fastq$", "", basename(forward_reads))
reverse_sample_names <- sub("_2\\.fastq$", "", basename(reverse_reads))

print(paste("Found", length(sample_names), "paired samples:"))
print(sample_names)

# Function: Plot quality profiles for a subset of samples
plot_quality_profiles <- function(files, n = 4L) {
  # Plot only first n files to avoid overwhelming output
  files_to_plot <- head(files, n)
  if (length(files_to_plot) > 0L) {
    plotQualityProfile(files_to_plot)
  }
}

# Function: Summarize FASTQ file contents
summarize_fastq <- function(fastq, n = 5L, nchar = 25L) {
  sequences <- as.character(sread(fastq)) # Extract DNA sequences as text
  read_count <- min(n, length(sequences)) # Don't exceed available reads

  list(
    length_distribution = table(width(fastq)), # Histogram of read lengths
    read_ids = head(as.character(id(fastq)), n), # First n read identifiers
    sequence_starts = substring(sequences[seq_len(read_count)], 1, nchar) # First nchar bases
  )
}

# Quality profiles drive the trimming decision, so inspect a few samples first.
plot_quality_profiles(forward_reads)
plot_quality_profiles(reverse_reads)

# Before filtering, confirm that forward and reverse reads are not identical (they should be different - sequencing from opposite ends).
raw_forward_fastq <- readFastq(forward_reads[[1]])
raw_reverse_fastq <- readFastq(reverse_reads[[1]])

raw_forward_summary <- summarize_fastq(raw_forward_fastq)
raw_reverse_summary <- summarize_fastq(raw_reverse_fastq)

# Print summaries to confirm files are correct and paired as expected before proceeding with filtering and DADA2 steps.
cat("\nForward read length distribution:\n")
print(raw_forward_summary$length_distribution)

cat("\nReverse read length distribution:\n")
print(raw_reverse_summary$length_distribution)

cat("\nForward read IDs (first 3):\n")
print(head(raw_forward_summary$read_ids, 3))

cat("\nReverse read IDs (first 3):\n")
print(head(raw_reverse_summary$read_ids, 3))

# If this count is high, files might be duplicated or mislabeled
paired_check_count <- min(100L, length(raw_forward_fastq), length(raw_reverse_fastq))
identical_reads <- sum(
  as.character(sread(raw_forward_fastq))[seq_len(paired_check_count)] ==
    as.character(sread(raw_reverse_fastq))[seq_len(paired_check_count)]
)
cat("\nNumber of identical forward/reverse sequences (should be ~0):", identical_reads, "\n")

use_multithread <- FALSE # Windows users should keep this FALSE for DADA2.

# These values are dataset-specific; revisit them if the primer set or quality profile changes.
trim_left <- c(19, 21)
trunc_len <- c(150, 130)

filtered_forward <- file.path(filtered_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtered_reverse <- file.path(filtered_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtered_forward) <- sample_names
names(filtered_reverse) <- sample_names

filtering_stats <- filterAndTrim(
  forward_reads, filtered_forward,
  reverse_reads, filtered_reverse,
  trimLeft = trim_left, # Remove primer sequences (19 bases for forward, 21 for reverse)
  truncLen = trunc_len, # Truncate reads at positions where quality drops (150 for forward, 130 for reverse)
  maxN = 0, # Discard reads with any ambiguous bases (N)
  maxEE = c(2, 5), # Discard reads with expected errors > 2 for forward and > 5 for reverse; this is a more nuanced quality filter than just average Q-score
  truncQ = 2, # Truncate reads at the first base with quality score ≤ 2; this helps remove low-quality tails
  rm.phix = TRUE, # Remove PhiX control sequences (Illumina spike-in)
  compress = TRUE, # Save output as .gz to save disk space
  multithread = use_multithread
)
# Display filtering results
# columns: reads.in, reads.out (for each sample)
cat("\nFiltering statistics (first 6 samples):\n")
print(head(filtering_stats))
saveRDS(filtering_stats, file.path(filtered_dir, "filtering_output.rds"))

# Check the first few filtered forward reads to confirm trimming behaved as expected.
filtered_forward_fastq <- readFastq(filtered_forward[[1]])
filtered_forward_summary <- summarize_fastq(filtered_forward_fastq)
print(filtered_forward_summary$sequence_starts)
print(filtered_forward_summary$length_distribution)

# DADA2's core innovation: learns the error profile from the data
# Why this matters: Distinguishes true biological variants from sequencing errors
# Traditional OTU methods: cluster at 97% similarity (arbitrary threshold)
# DADA2: Uses error model to infer exact sequences (single-nucleotide resolution)

cat("Learning error rates from filtered reads...\n")
cat("(This step examines many reads to build error model - be patient)\n\n")

err_forward <- learnErrors(filtered_forward, multithread = use_multithread)
err_reverse <- learnErrors(filtered_reverse, multithread = use_multithread)
# Learn the error model from the filtered reads; rerun this step after changing trimming settings.

# Black points = observed error rates
# Black line = estimated error rates (should fit points reasonably)
# Red line = expected error rates based on quality scores
plotErrors(err_forward, nominalQ = TRUE)
plotErrors(err_reverse, nominalQ = TRUE)

# Infer ASVs with the DADA2 model.
# pool = "pseudo": Compromise between speed and sensitivity
#   - "FALSE": Process each sample independently (fast, less sensitive)
#   - "TRUE": Pool all samples together (slow, most sensitive for rare variants)
#   - "pseudo": Partial pooling (good balance for most datasets)

cat("Running DADA2 algorithm on forward reads...\n")
dada_forward <- dada(filtered_forward,
  err = err_forward,
  pool = "pseudo", multithread = use_multithread
)

cat("Running DADA2 algorithm on reverse reads...\n")
dada_reverse <- dada(filtered_reverse,
  err = err_reverse,
  pool = "pseudo", multithread = use_multithread
)

# Inspect results for first sample
cat("\nDada-class object for first sample (forward reads):\n")
print(dada_forward[[1]])

# Merge paired reads before building the ASV table.
mergers <- mergePairs(dada_forward, filtered_forward, dada_reverse, filtered_reverse, verbose = TRUE)
print(head(mergers[[1]]))

# Keep only sequences that match the expected V4 amplicon length for this dataset.
sequence_table <- makeSequenceTable(mergers)
print(table(nchar(getSequences(sequence_table))))

expected_v4_lengths <- 250:256 # Update this if a different region was amplified.
v4_sequence_table <- sequence_table[, nchar(colnames(sequence_table)) %in% expected_v4_lengths]

# Remove chimeric sequences that likely arose during PCR amplification; this step is critical for accurate diversity estimates.
nonchim_sequence_table <- removeBimeraDenovo(
  v4_sequence_table,
  method = "consensus",
  multithread = use_multithread,
  verbose = TRUE
)
print(sum(nonchim_sequence_table) / sum(v4_sequence_table))
print(dim(nonchim_sequence_table))
saveRDS(nonchim_sequence_table, file.path(filtered_dir, "ASV_table.rds"))

count_unique_reads <- function(x) sum(getUniques(x))

# Create a tracking table to summarize how many reads were retained at each step of the pipeline for each sample;
# this helps identify where losses occur and if any samples had unusually high dropouts.
tracking_table <- cbind(
  filtering_stats,
  sapply(dada_forward, count_unique_reads),
  sapply(dada_reverse, count_unique_reads),
  sapply(mergers, count_unique_reads),
  rowSums(nonchim_sequence_table)
)
colnames(tracking_table) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(tracking_table) <- sample_names
print(tracking_table)

# Assign taxonomy with the SILVA reference used by this project.
taxonomy <- assignTaxonomy(
  nonchim_sequence_table,
  reference_fasta,
  multithread = use_multithread
)

taxonomy_preview <- taxonomy
rownames(taxonomy_preview) <- NULL
head(taxonomy_preview)

# Species-level assignment is not guaranteed; many ASVs will stop at genus or family.
# ---- PIPELINE COMPLETE ----
cat("\n")
cat("============================================\n")
cat("   PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("============================================\n")
cat("\nKey outputs saved in:", filtered_dir, "\n")
cat("  1. ASV_table.rds - abundance of each ASV per sample\n")
cat("  2. taxonomy.rds - taxonomic assignments (Kingdom → Species)\n")
cat("  3. filtering_output.rds - quality filtering statistics\n")
cat("\nNext steps:\n")
cat("  - Load ASV table and taxonomy into phyloseq for diversity analysis\n")
cat("  - Create bar plots of taxonomic composition\n")
cat("  - Calculate alpha/beta diversity metrics\n")
cat("  - Test for differential abundance between sample groups\n")
cat("\n")

# =============================================================================
# END OF PIPELINE
# =============================================================================
