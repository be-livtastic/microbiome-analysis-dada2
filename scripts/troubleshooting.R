# =============================================================================
# DADA2 Pipeline Troubleshooting and Diagnostic Checks
# =============================================================================
# Purpose:  Diagnose common issues with 16S sequence processing
#           Use this when something goes wrong in the main pipeline
library(ShortRead) # FASTQ file inspection
library(Biostrings) # DNA sequence manipulation
library(here) # Portable file paths

# ---- SETUP ----
project_dir <- here::here()

# Reference database paths (for validation)
silva_genus_reference <- here::here(
    "data", "external", "reference",
    "silva_nr99_v138.1_train_set.fa.gz"
)
silva_species_reference <- here::here(
    "data", "external", "reference",
    "silva_nr99_v138.1_wSpecies_train_set.fa.gz"
)

# ---- CHECK FOR REQUIRED OBJECTS ----
cat("=== Checking for required pipeline variables ===\n")

# These objects should exist if you've run dada2_pipeline.R
required_objects <- c("fnFs", "fnRs", "filtFs", "filtRs")

missing_objects <- required_objects[!vapply(required_objects, exists,
    logical(1),
    inherits = TRUE
)]

if (length(missing_objects) > 0L) {
    stop(
        "This script expects the following objects to already exist: ",
        paste(missing_objects, collapse = ", "),
        ".\n\nRun the main DADA2 pipeline first (dada2_pipeline.R),",
        " or define these file paths manually before using these diagnostic checks."
    )
}

cat("✓ All required objects found\n\n")

# ---- HELPER FUNCTIONS ----

# Preview first few sequences from a FASTQ file
preview_reads <- function(fastq_file, n = 5L, width = 25L) {
    sequences <- as.character(sread(readFastq(fastq_file)))
    substring(sequences[seq_len(min(n, length(sequences)))], 1, width)
}

# Show read IDs (useful for checking file pairing)
show_read_ids <- function(fastq_file, n = 3L) {
    head(as.character(id(readFastq(fastq_file))), n)
}

# Show length distribution (detects truncation issues)
show_length_distribution <- function(fastq_file) {
    table(width(readFastq(fastq_file)))
}

# ---- DIAGNOSTIC CHECK 1: PAIRED-END ORIENTATION ----
cat("=== CHECK 1: Paired-End Read Orientation ===\n")
cat("Purpose: Verify forward and reverse reads are properly oriented\n\n")

# Background: In Illumina paired-end sequencing:
# - Forward reads sequence from one end of the fragment
# - Reverse reads sequence from the OTHER end
# - Reverse reads are stored in FORWARD orientation in FASTQ
# - Their reverse complement should resemble the forward reads

# Load first filtered sample
filtered_forward_fastq <- readFastq(filtFs[[1]])
filtered_reverse_fastq <- readFastq(filtRs[[1]])

# Get sequences
forward_preview <- as.character(sread(filtered_forward_fastq)[seq_len(min(5L, length(filtered_forward_fastq)))])

# Reverse complement the reverse reads
# Why? Because they're stored forward but sequence backward
reverse_rc_preview <- as.character(reverseComplement(
    sread(filtered_reverse_fastq)[seq_len(min(5L, length(filtered_reverse_fastq)))]
))

cat("Forward reads start (first 30 bp):\n")
print(substring(forward_preview, 1, 30))

cat("\nReverse reads reverse-complemented start (first 30 bp):\n")
print(substring(reverse_rc_preview, 1, 30))

cat("\nReverse reads reverse-complemented end (last 30 bp):\n")
print(substring(reverse_rc_preview, nchar(reverse_rc_preview) - 29, nchar(reverse_rc_preview)))

cat("\n✓ What to look for:\n")
cat("  - If forward and reverse RC start match → reads overlap (GOOD)\n")
cat("  - If they don't match → check truncLen parameters (may be too short)\n\n")

# ---- DIAGNOSTIC CHECK 2: FILE PAIRING ----
cat("=== CHECK 2: File Pairing Verification ===\n")
cat("Purpose: Confirm forward and reverse files match correctly\n\n")

cat("Forward read IDs:\n")
print(show_read_ids(fnFs[[1]]))

cat("\nReverse read IDs:\n")
print(show_read_ids(fnRs[[1]]))

cat("\n✓ What to look for:\n")
cat("  - Read IDs should be identical except for /1 vs /2 suffix\n")
cat("  - If completely different → files are mismatched!\n\n")

# ---- DIAGNOSTIC CHECK 3: RAW SEQUENCE CONTENT ----
cat("=== CHECK 3: Raw Sequence Inspection ===\n")
cat("Purpose: Check if primers are present in raw data\n\n")

cat("Forward read preview (first 5 reads, first 25 bp):\n")
print(preview_reads(fnFs[[1]]))

cat("\nReverse read preview (first 5 reads, first 25 bp):\n")
print(preview_reads(fnRs[[1]]))

# Check reverse reads specifically
# The reverse reads should still show the primer near the start
cat("\nFirst 10 reverse reads, first 10 bases:\n")
print(substring(as.character(sread(readFastq(fnRs[[1]])))[1:10], 1, 10))

cat("\n✓ What to look for:\n")
cat("  - Consistent sequence starts → primers present\n")
cat("  - Random starts → primers already removed or data has issues\n\n")

# ---- DIAGNOSTIC CHECK 4: FILTERED SEQUENCE CONTENT ----
cat("=== CHECK 4: Filtered Sequence Inspection ===\n")
cat("Purpose: Verify trimming removed primers correctly\n\n")

cat("Filtered forward reads (first 5, first 25 bp):\n")
print(preview_reads(filtFs[[1]]))

cat("\nFiltered reverse reads (first 5, first 25 bp):\n")
print(preview_reads(filtRs[[1]]))

cat("\n✓ What to look for:\n")
cat("  - Sequence starts should differ from raw reads\n")
cat("  - If identical to raw → trimLeft parameter not working\n\n")

# ---- DIAGNOSTIC CHECK 5: READ LENGTH DISTRIBUTIONS ----
cat("=== CHECK 5: Read Length Distribution ===\n")
cat("Purpose: Check if truncation worked as expected\n\n")

cat("Filtered forward read length distribution:\n")
print(show_length_distribution(filtFs[[1]]))

cat("\nFiltered reverse read length distribution:\n")
print(show_length_distribution(filtRs[[1]]))

cat("\n✓ What to look for:\n")
cat("  - All reads should be exactly truncLen length\n")
cat("  - If variable lengths → truncLen parameter not applied correctly\n\n")

# ---- DIAGNOSTIC CHECK 6: REFERENCE DATABASE VALIDATION ----
cat("=== CHECK 6: Reference Database Validation ===\n")
cat("Purpose: Ensure reference files are correct and accessible\n\n")

cat("Reference file sizes:\n")
cat("Genus-level database:", file.info(silva_genus_reference)$size, "bytes\n")
cat("Species-level database:", file.info(silva_species_reference)$size, "bytes\n")

cat("\nSILVA files in project directory:\n")
silva_files <- list.files(project_dir,
    pattern = "silva",
    recursive = TRUE, full.names = FALSE
)
print(silva_files)

cat("\n✓ What to look for:\n")
cat("  - Genus database should be ~130-150 MB compressed\n")
cat("  - Species database should be ~70-80 MB compressed\n")
cat("  - Files in data/external/reference/ directory\n\n")

# ---- DIAGNOSTIC CHECK 7: COMPUTE OVERLAP LENGTH ----
cat("=== CHECK 7: Paired-End Overlap Calculation ===\n")
cat("Purpose: Estimate expected overlap between forward and reverse reads\n\n")

# If you ran the main pipeline, these variables should exist
if (exists("trunc_len") && exists("expected_v4_lengths")) {
    # Expected amplicon length (from V4 region)
    amplicon_length <- mean(expected_v4_lengths)

    # Expected overlap = amplicon length - (forward truncLen + reverse truncLen)
    # Overlap must be ≥12 bp for mergePairs() to work
    expected_overlap <- amplicon_length - (trunc_len[1] + trunc_len[2])

    cat("Parameters:\n")
    cat("  Expected amplicon length:", amplicon_length, "bp\n")
    cat("  Forward truncLen:", trunc_len[1], "bp\n")
    cat("  Reverse truncLen:", trunc_len[2], "bp\n")
    cat("\nCalculated expected overlap:", expected_overlap, "bp\n")

    if (expected_overlap >= 20) {
        cat("✓ GOOD: Plenty of overlap for merging\n")
    } else if (expected_overlap >= 12) {
        cat("⚠ MARGINAL: Minimal overlap, merge rate may be low\n")
    } else {
        cat("✗ BAD: Insufficient overlap! Reduce truncLen or merging will fail\n")
    }
} else {
    cat("⚠ Variables 'trunc_len' or 'expected_v4_lengths' not found\n")
    cat("  Run main pipeline first to get these values\n")
}

cat("\n")

# ---- DIAGNOSTIC CHECK 8: QUALITY SCORE SUMMARY ----
cat("=== CHECK 8: Quality Score Summary ===\n")
cat("Purpose: Quick stats on quality after filtering\n\n")

# Read first filtered sample
filt_fq <- readFastq(filtFs[[1]])
quality_scores <- as(quality(filt_fq), "matrix")

cat("Quality score statistics (first sample):\n")
cat("  Mean quality:", round(mean(quality_scores), 1), "\n")
cat("  Median quality:", median(quality_scores), "\n")
cat("  Min quality:", min(quality_scores), "\n")
cat("  Max quality:", max(quality_scores), "\n")

cat("\n✓ What to look for:\n")
cat("  - Mean ≥30 is excellent\n")
cat("  - Mean 25-30 is good\n")
cat("  - Mean <25 suggests aggressive filtering needed\n\n")

# ---- COMMON ISSUES AND SOLUTIONS ----
cat("=== Common Issues and Solutions ===\n\n")

cat("ISSUE 1: Low filtering efficiency (<70% reads retained)\n")
cat("  Possible causes:\n")
cat("    - Quality too poor → Reduce truncLen (keep less of read)\n")
cat("    - maxEE too stringent → Increase maxEE values\n")
cat("    - Wrong quality encoding → Check FASTQ format\n")
cat("  Solution: Plot quality profiles, adjust truncLen/maxEE\n\n")

cat("ISSUE 2: Low merge rate (<50% forward reads merged)\n")
cat("  Possible causes:\n")
cat("    - Insufficient overlap → See CHECK 7 above\n")
cat("    - Poor quality at read ends → Increase truncLen\n")
cat("    - Mismatch between forward/reverse reads\n")
cat("  Solution: Calculate expected overlap, adjust truncLen upward\n\n")

cat("ISSUE 3: High chimera rate (>30% sequences removed)\n")
cat("  Possible causes:\n")
cat("    - Too many PCR cycles during library prep\n")
cat("    - Poor quality template DNA\n")
cat("    - Contamination\n")
cat("  Solution: Check wet lab protocol, consider filtering samples\n\n")

cat("ISSUE 4: Unexpected ASV lengths\n")
cat("  Possible causes:\n")
cat("    - Wrong amplicon size range (expected_v4_lengths)\n")
cat("    - Different 16S region than expected\n")
cat("    - Non-specific amplification\n")
cat("  Solution: Check length distribution, update expected_v4_lengths\n\n")

cat("ISSUE 5: Poor taxonomic assignments (many NAs)\n")
cat("  Possible causes:\n")
cat("    - Wrong reference database (not 16S)\n")
cat("    - Primers not properly removed (affecting alignment)\n")
cat("    - Novel/uncultured organisms\n")
cat("  Solution: Check primer trimming, verify SILVA version\n\n")

cat("ISSUE 6: 'Error in learnErrors: Not enough reads'\n")
cat("  Possible causes:\n")
cat("    - Too few samples (<3)\n")
cat("    - Too aggressive filtering (all reads removed)\n")
cat("  Solution: Relax filtering, increase sample count, or pool samples\n\n")

# ---- END OF DIAGNOSTICS ----
cat("\n")
cat("============================================\n")
cat("   DIAGNOSTIC CHECKS COMPLETE\n")
cat("============================================\n")
cat("\nNext steps:\n")
cat("  1. Review each diagnostic output above\n")
cat("  2. Identify which check(s) reveal issues\n")
cat("  3. Apply relevant solutions from 'Common Issues'\n")
cat("  4. Re-run main pipeline with adjusted parameters\n")
cat("  5. If still stuck, consult DADA2 tutorial or post to forum\n")
cat("\n")
