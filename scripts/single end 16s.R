# NOTE:     Single-end provides less taxonomic resolution than paired-end
#           (shorter sequences = less discriminatory power)
#           Use paired-end pipeline when both forward/reverse available
# =============================================================================

# ---- LOAD REQUIRED PACKAGES ----

set.seed(1) # Reproducibility

# ---- DEFINE PROJECT PATHS ----
project_dir <- here::here()
raw_dir <- here::here("data", "raw", "fastq")
filtered_dir <- here::here("data", "processed", "dada2")

# Path to the directory (not specific file)
path <- raw_dir

# SILVA reference databases
# For single-end, we can optionally do two-step taxonomy:
# 1. Genus-level assignment (faster, using train set)
# 2. Species-level refinement (slower, using species database)
silva_genus_reference <- here::here(
  "data", "external", "reference",
  "silva_nr99_v138.1_train_set.fa.gz"
)
silva_species_reference <- here::here(
  "data", "external", "reference",
  "silva_species_assignment_v138.1.fa.gz"
)

# ---- VALIDATE REQUIRED FILES ----
if (!dir.exists(raw_dir)) {
  stop("Raw FASTQ directory not found: ", raw_dir)
}

if (!file.exists(silva_genus_reference)) {
  stop("Genus-level SILVA reference not found: ", silva_genus_reference)
}

if (!file.exists(silva_species_reference)) {
  stop("Species-level SILVA reference not found: ", silva_species_reference)
}

# Create output directory
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

fastq_files <- list.files(raw_dir)
print(fastq_files)

# ---- DEFINE FILE PATHS AND SAMPLE NAMES ----
fnFs <- sort(list.files(raw_dir, pattern = "_1\\.fastq$", full.names = TRUE))
if (length(fnFs) == 0L) {
  stop("No forward FASTQ files were found in: ", raw_dir)
}

sample.names <- sub("_1\\.fastq$", "", basename(fnFs))

# For single-end data, we only have forward reads. Reverse read paths are not defined.
plot_quality_profiles <- function(files, n = 4L) {
  files_to_plot <- head(files, n)
  if (length(files_to_plot) > 0L) {
    plotQualityProfile(files_to_plot)
  }
}

summarize_fastq <- function(fastq, n = 5L, nchar = 25L) {
  sequences <- as.character(sread(fastq))
  read_count <- min(n, length(sequences))

  list(
    length_distribution = table(width(fastq)),
    read_ids = head(as.character(id(fastq)), n),
    sequence_starts = substring(sequences[seq_len(read_count)], 1, nchar)
  )
}

# Inspect quality before choosing trimming settings.
plot_quality_profiles(fnFs)

raw_forward_fastq <- readFastq(fnFs[[1]])
raw_forward_summary <- summarize_fastq(raw_forward_fastq)
print(raw_forward_summary$length_distribution)
print(raw_forward_summary$read_ids)

use_multithread <- FALSE # Keep FALSE on Windows for predictable behavior.

# Primer trimming and truncation are dataset-specific; revisit these values if the run changes.
trim_left <- c(19)
trunc_len <- c(130)

filtFs <- file.path(filtered_dir, paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out_SE <- filterAndTrim(
  fnFs, filtFs,
  trimLeft = trim_left,
  truncLen = trunc_len,
  maxN = 0,
  maxEE = c(2),
  truncQ = 2,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = use_multithread
)

print(head(out_SE))
saveRDS(out_SE, file.path(filtered_dir, "filtering_output_SE.rds"))

# Check the first filtered reads to confirm the trimmed sequence start looks sensible.
filtered_forward_fastq <- readFastq(filtFs[[1]])
filtered_forward_summary <- summarize_fastq(filtered_forward_fastq)
print(filtered_forward_summary$sequence_starts)
print(filtered_forward_summary$length_distribution)

errF_SE <- learnErrors(filtFs, multithread = use_multithread)
plotErrors(errF_SE, nominalQ = TRUE)

dadaFs_SE <- dada(filtFs, err = errF_SE, pool = "pseudo", multithread = use_multithread)
print(dadaFs_SE[[1]])

# For single-end data, the ASV table comes directly from denoised forward reads.
seqtab_SE <- makeSequenceTable(dadaFs_SE)
print(table(nchar(getSequences(seqtab_SE))))
print(dim(seqtab_SE))

seqtab_SE.nochim <- removeBimeraDenovo(
  seqtab_SE,
  method = "consensus",
  multithread = use_multithread,
  verbose = TRUE
)
print(sum(seqtab_SE.nochim) / sum(seqtab_SE))

count_unique_reads <- function(x) sum(getUniques(x))
track_SE <- cbind(out_SE, sapply(dadaFs_SE, count_unique_reads), rowSums(seqtab_SE.nochim))
colnames(track_SE) <- c("input", "filtered", "denoised", "nonchim")
rownames(track_SE) <- sample.names
print(track_SE)

saveRDS(seqtab_SE.nochim, file.path(filtered_dir, "ASV_table_SE.rds"))

taxa_SE <- assignTaxonomy(
  seqtab_SE.nochim,
  silva_genus_reference,
  multithread = use_multithread,
  verbose = TRUE
)

saveRDS(taxa_SE, file.path(filtered_dir, "taxonomy_SE_genuslevel.rds"))

taxa_SE <- addSpecies(taxa_SE, silva_species_reference)

taxa_preview <- taxa_SE
rownames(taxa_preview) <- NULL
head(taxa_preview, 20)

saveRDS(taxa_SE, file.path(filtered_dir, "taxonomy_SE.rds"))
# ---- PIPELINE COMPLETE ----
cat("\n")
cat("============================================\n")
cat("   SINGLE-END PIPELINE COMPLETED!\n")
cat("============================================\n")
cat("\nKey outputs saved in:", filtered_dir, "\n")
cat("  1. ASV_table_SE.rds - abundance matrix\n")
cat("  2. taxonomy_SE.rds - taxonomic assignments (with species)\n")
cat("  3. taxonomy_SE_genuslevel.rds - taxonomic assignments (genus only)\n")
cat("  4. filtering_output_SE.rds - quality filtering stats\n")
cat("\nImportant notes about single-end data:\n")
cat("  - Shorter sequences (~130 bp) vs paired-end merged (~250 bp)\n")
cat("  - Lower taxonomic resolution (harder to distinguish closely related taxa)\n")
cat("  - Faster processing (no merging step)\n")
cat("  - Use when paired-end data unavailable or for quick exploratory analysis\n")
cat("\nNext steps: Same as paired-end pipeline\n")
cat("  - Load into phyloseq for diversity analysis\n")
cat("  - Calculate alpha/beta diversity\n")
cat("  - Test for differential abundance\n")
cat("\n")

# =============================================================================
# SINGLE-END vs PAIRED-END: Key Differences
# =============================================================================
#
# Workflow differences:
# - No mergePairs() step (only have forward reads)
# - Simpler tracking table (no "merged" column)
# - Shorter final sequences
#
# Parameter differences:
# - One truncLen value instead of two: truncLen = c(130)
# - One trimLeft value instead of two: trimLeft = c(19)
# - One maxEE value instead of two: maxEE = c(2)
#
# Biological implications:
# - Less discriminatory power (shorter = less info)
# - Species-level assignment less reliable
# - Still valid for community composition analysis
# - Good for genus/family level questions
#
# When to use single-end:
# - Paired-end data unavailable
# - Budget constraints (single-end cheaper)
# - Exploratory analysis / pilot study
# - Very short amplicons (<200 bp total)
# =============================================================================
