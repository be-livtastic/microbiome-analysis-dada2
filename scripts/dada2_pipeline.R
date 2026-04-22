# ---- LOAD REQUIRED PACKAGES ----

# Set random seed for reproducibility
# Any random operations (e.g., pseudo-pooling) will give same results each run
set.seed(1)

# =============================================================================
# PIPELINE MODE CONFIGURATION
# =============================================================================
# Select which primer removal strategy to use before DADA2 denoising:
#
#   "standard"  -- Removes primers via trimLeft inside filterAndTrim.
#                  Simpler; works well when primers are at a fixed, consistent
#                  position in every read (typical for 16S V4 datasets).
#
#   "cutadapt"  -- Removes primers using cutadapt (external tool, called via
#                  WSL on Windows). Recommended when:
#                    * Primers are not at a fixed position (e.g., ITS amplicons)
#                    * You want --discard-untrimmed to drop off-target reads
#                    * Primers span degenerate or variable-length regions
#                  Requires cutadapt installed in WSL: pip install cutadapt
#
# Everything downstream of primer removal (error learning, denoising, merging,
# taxonomy) is shared between both modes.
pipeline_mode <- "cutadapt"  # Change to "standard" to use the standard workflow

# ---- DEFINE PROJECT PATHS ----
# Using here::here() makes paths work regardless of working directory

project_dir <- here::here() # Root directory of the R project

# Where raw FASTQ files are stored
raw_dir <- here::here("data", "raw", "fastq") #or "NCBI_SRA_fastq"

# Where filtered/processed files will be saved
filtered_dir <- here::here("data", "processed", "dada2")

# SILVA reference database for taxonomic classification
# This file assigns taxonomy (Kingdom → Species) to our sequences
reference_fasta <- here::here(
  "data", "external", "reference",
  "silva_species_assignment_v138.1.fa.gz"
)


# =============================================================================
# DISCOVER AND PAIR FASTQ FILES
# =============================================================================

fastq_files   <- list.files(raw_dir)
print(fastq_files)

forward_reads <- sort(list.files(raw_dir, pattern = "_1\\.fastq$", full.names = TRUE))
reverse_reads <- sort(list.files(raw_dir, pattern = "_2\\.fastq$", full.names = TRUE))


# ---- PRE-CHECKS ----
# Basic checks to confirm files are found and paired correctly before processing
if (length(forward_reads) == 0L || length(reverse_reads) == 0L) {
  stop("No paired FASTQ files were found in: ", raw_dir)
}
if (length(forward_reads) != length(reverse_reads)) {
  stop("The number of forward and reverse FASTQ files do not match.")
}

# Extract sample names by removing the read suffix and file extension;
# these are used for naming all output files consistently
sample_names         <- sub("_1\\.fastq$", "", basename(forward_reads))
reverse_sample_names <- sub("_2\\.fastq$", "", basename(reverse_reads))

# Extra guard: confirm forward/reverse sample order is aligned
# sort() above should guarantee this, but explicit confirmation prevents silent mispairing
if (!all(sample_names == reverse_sample_names)) {
  stop(
    "Sample name mismatch between forward and reverse reads. ",
    "Check that files are named consistently and sort correctly."
  )
}

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

# Function: Count how many reads contain a given primer (in any orientation)
# Used to verify primer presence before cutadapt and confirm removal after
primer_hits <- function(primer, fastq_file) {
  nhits <- vcountPattern(primer, sread(readFastq(fastq_file)), fixed = FALSE)
  return(sum(nhits > 0))
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
# For NCBI SRA downloads: if identical_reads > 0, download from ENA instead
paired_check_count <- min(100L, length(raw_forward_fastq), length(raw_reverse_fastq))
identical_reads <- sum(
  as.character(sread(raw_forward_fastq))[seq_len(paired_check_count)] ==
    as.character(sread(raw_reverse_fastq))[seq_len(paired_check_count)]
)
cat("\nNumber of identical forward/reverse sequences (should be ~0):", identical_reads, "\n")

# =============================================================================
# PRIMER SEQUENCES AND ORIENTATIONS
# =============================================================================
# These are the primers used to amplify the V4 region of the 16S rRNA gene
# 515FB = GTGYCAGCMGCCGCGGTAA  (19 bp, IUPAC degenerate: Y = C/T, M = A/C)
# 806RB = GGACTACNVGGGTWTCTAAT (20 bp, IUPAC degenerate: N = any, V = A/C/G, W = A/T)

FWD_primer <- "GTGYCAGCMGCCGCGGTAA"   # 515FB: 19 bp
REV_primer <- "GGACTACNVGGGTWTCTAAT"  # 806RB: 20 bp

# Function: Generate all orientations of a primer sequence
# Primers can appear in unexpected orientations (e.g., read-through into the adapter).
# We check all four orientations to catch any that survived into the reads.
allOrients <- function(primer) {
  require(Biostrings)
  dna    <- DNAString(primer) # Biostrings works with DNAString objects, not plain character strings
  orients <- c(
    Forward    = dna,
    Complement = Biostrings::complement(dna),
    Reverse    = Biostrings::reverse(dna),
    RevComp    = Biostrings::reverseComplement(dna)
  )
  return(sapply(orients, toString)) # Convert back to character vector for printing / pattern matching
}

FWD.orients <- allOrients(FWD_primer)
REV.orients <- allOrients(REV_primer)
cat("Forward primer orientations:\n")
print(FWD.orients)
cat("\nReverse primer orientations:\n")
print(REV.orients)


# =============================================================================
# SHARED OUTPUT PATH DEFINITIONS
# =============================================================================
# These filtered file paths are written by both pipeline modes.
# In "standard" mode : filterAndTrim writes here directly from raw reads
# In "cutadapt" mode : filterAndTrim writes here from cutadapt-trimmed reads
#                      (a separate N-filter and cutadapt step runs first)

filtered_forward <- file.path(filtered_dir, paste0(sample_names, "_F_filt.fastq.gz"))
filtered_reverse <- file.path(filtered_dir, paste0(sample_names, "_R_filt.fastq.gz"))
names(filtered_forward) <- sample_names
names(filtered_reverse) <- sample_names

# Windows requires multithread = FALSE for all DADA2 steps (no forking support)
use_multithread <- FALSE


# =============================================================================
# PIPELINE BRANCH A: STANDARD MODE
# Primer removal via trimLeft inside filterAndTrim
# =============================================================================
if (pipeline_mode == "standard") {

  cat("\n--- Running STANDARD pipeline (trimLeft primer removal) ---\n\n")

  # These values are dataset-specific; revisit them if the primer set or quality profile changes.
  # BONUS: If reads are the same length, use Figaro to optimise trimming parameters
  # based on quality profiles and expected amplicon length:
  # https://github.com/Zymo-Research/figaro#figaro
  trim_left <- c(19L, 20L)    # 515F = 19 bp; 806R = 20 bp
                               # Add linker bp to REV if your reads carry a linker sequence
  trunc_len <- c(150L, 150L)  # Truncate at positions where quality drops below Q25
                               # Set per quality plots; must leave enough overlap for merging
                               # Overlap = truncLenF + truncLenR - expected_amplicon_length
                               # For V4 (~253 bp): 150 + 150 - 253 = 47 bp overlap (minimum ~20 bp needed)

  filtering_stats <- filterAndTrim(
    forward_reads, filtered_forward,
    reverse_reads, filtered_reverse,
    trimLeft    = trim_left,   # Remove primer sequences from 5' end before quality filtering
    truncLen    = trunc_len,   # Truncate reads to uniform length; also discards reads shorter than truncLen
    maxN        = 0,           # Discard reads with any ambiguous bases (N); DADA2 cannot handle Ns
    maxEE       = c(2, 5),     # Max expected errors: stricter for forward (higher quality), relaxed for reverse
    truncQ      = 2,           # Truncate at first base with quality score <= 2; removes low-quality tails
    rm.phix     = TRUE,        # Remove PhiX control sequences (Illumina spike-in; contaminant if present)
    compress    = TRUE,        # Save output as .gz to save disk space
    multithread = use_multithread
  )

  # Post-filter primer check: counts should be ~0 if trimLeft removed primers correctly
  # If counts remain high, primers may be at a variable position - consider switching to cutadapt mode
  primer_counts_post_filter <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primer_hits, fastq_file = filtered_forward[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primer_hits, fastq_file = filtered_reverse[[1]]),
    REV.ForwardReads = sapply(REV.orients, primer_hits, fastq_file = filtered_forward[[1]]),
    REV.ReverseReads = sapply(REV.orients, primer_hits, fastq_file = filtered_reverse[[1]])
  )
  cat("\nPrimer hit counts after standard filtering (should be ~0):\n")
  print(primer_counts_post_filter)

  # Hand off to shared downstream steps
  reads_for_dada_F <- filtered_forward
  reads_for_dada_R <- filtered_reverse
}


# =============================================================================
# PIPELINE BRANCH B: CUTADAPT MODE
# Primer removal via cutadapt (WSL), followed by DADA2 quality filtering
# =============================================================================
# The dataset is a set of Illumina-sequenced paired-end FASTQ files that have
# been split ("demultiplexed") by sample and from which the barcodes have
# already been removed. We therefore skip demultiplexing and go straight to
# primer removal and quality filtering.
if (pipeline_mode == "cutadapt") {

  cat("\n--- Running CUTADAPT pipeline (external primer removal via WSL) ---\n\n")

  # ------------------------------------------------------------------
  # Step B1: N-filter only (no trimming or truncation yet)
  # ------------------------------------------------------------------
  # The presence of ambiguous bases (Ns) makes accurate primer mapping difficult.
  # Pre-filter to remove N-containing reads ONLY before running cutadapt.
  # All quality trimming (truncLen, maxEE) happens AFTER cutadapt in Step B3.

  n_filtered_dir <- here::here("data", "processed", "dada2_n_filtered")
  dir.create(n_filtered_dir, recursive = TRUE, showWarnings = FALSE)

  n_filtered_forward <- file.path(n_filtered_dir, paste0(sample_names, "_F_nfilt.fastq.gz"))
  n_filtered_reverse <- file.path(n_filtered_dir, paste0(sample_names, "_R_nfilt.fastq.gz"))
  names(n_filtered_forward) <- sample_names
  names(n_filtered_reverse) <- sample_names

  n_filtering_stats <- filterAndTrim(
    forward_reads, n_filtered_forward,
    reverse_reads, n_filtered_reverse,
    maxN        = 0,   # Only filter: remove reads with ambiguous bases before primer mapping
    compress    = TRUE,
    multithread = use_multithread
  )
  cat("N-filtering statistics (reads before/after N removal):\n")
  print(head(n_filtering_stats))

  # ------------------------------------------------------------------
  # Step B2: Confirm primers are present before cutadapt
  # ------------------------------------------------------------------
  # If counts are already ~0 here, primers may already be absent from your data.
  # In that case, switch back to standard mode with trim_left = c(0, 0).
  primer_counts_pre_cut <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primer_hits, fastq_file = n_filtered_forward[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primer_hits, fastq_file = n_filtered_reverse[[1]]),
    REV.ForwardReads = sapply(REV.orients, primer_hits, fastq_file = n_filtered_forward[[1]]),
    REV.ReverseReads = sapply(REV.orients, primer_hits, fastq_file = n_filtered_reverse[[1]])
  )
  cat("\nPrimer hit counts in N-filtered reads (before cutadapt):\n")
  print(primer_counts_pre_cut)

  # ------------------------------------------------------------------
  # Step B3: Run cutadapt via WSL to remove primers
  # ------------------------------------------------------------------

  # Windows path converter: R uses C:/path/to/file, WSL expects /mnt/c/path/to/file
  # Also normalise backslashes that Windows sometimes produces
  to_wsl_path <- function(path) {
    path <- normalizePath(path, winslash = "/", mustWork = FALSE)
    gsub("^([A-Za-z]):/", "/mnt/\\L\\1/", path, perl = TRUE)
  }

  # Verify cutadapt is accessible before attempting to run the loop
  cutadapt_version <- tryCatch(
    system2("wsl", args = c("cutadapt", "--version"), stdout = TRUE, stderr = TRUE),
    error = function(e) character(0)
  )
  if (length(cutadapt_version) == 0L || grepl("not found", cutadapt_version[1], ignore.case = TRUE)) {
    stop(
      "cutadapt not found via WSL. ",
      "Install it in your WSL environment: pip install cutadapt"
    )
  }
  cat("\ncutadapt version:", cutadapt_version[1], "\n")

  # Output directory for cutadapt-trimmed (but not yet quality-filtered) files
  cutadapt_dir <- here::here("data", "processed", "dada2_cutadapt")
  dir.create(cutadapt_dir, recursive = TRUE, showWarnings = FALSE)

  fnFs.cut <- file.path(cutadapt_dir, paste0(sample_names, "_F_cut.fastq.gz"))
  fnRs.cut <- file.path(cutadapt_dir, paste0(sample_names, "_R_cut.fastq.gz"))

  # Generate reverse complements for read-through trimming.
  # cutadapt needs to handle primers appearing at BOTH ends of a read:
  #   R1 (forward read): FWD primer at 5' end (-g); RC of REV at 3' end if read-through (-a)
  #   R2 (reverse read): REV primer at 5' end (-G); RC of FWD at 3' end if read-through (-A)
  FWD.RC <- dada2:::rc(FWD_primer)
  REV.RC <- dada2:::rc(REV_primer)

  R1.flags <- paste("-g", FWD_primer, "-a", REV.RC)
  R2.flags <- paste("-G", REV_primer, "-A", FWD.RC)

  # On Linux, paths need no conversion — use them directly 
  n_filtered_forward_wsl <- to_wsl_path(n_filtered_forward)
  n_filtered_reverse_wsl <- to_wsl_path(n_filtered_reverse)
  fnFs.cut_wsl <- to_wsl_path(fnFs.cut)
  fnRs.cut_wsl <- to_wsl_path(fnRs.cut)

  
  # Run cutadapt for each sample
  for (i in seq_along(n_filtered_forward)) {
    cat("Trimming primers for sample:", sample_names[i], "\n")
    system2(
      "wsl",
      args = c(
        "cutadapt",
        R1.flags,
        R2.flags,
        "-n", "2",                      # Run adapter trimming twice to catch read-through amplicons
        "--discard-untrimmed",          # Discard read pairs where the primer was NOT found
                                        # (these are likely non-target amplification products)
  "-o", fnFs.cut_wsl[i],
        "-p", fnRs.cut_wsl[i],
        n_filtered_forward_wsl[i],
        n_filtered_reverse_wsl[i]
      ),
      stdout = TRUE,
      stderr = TRUE
    )
    
    cutadapt_logs[[i]] <- cmd
  }
  
  # ---- 11. Save logs for reproducibility ----
  log_file <- file.path(cutadapt_dir, "cutadapt_run_log.txt")
  writeLines(unlist(cutadapt_logs), con = log_file)
  
  cat("\nCutadapt run complete. Logs saved to:\n", log_file, "\n")

  # Verify primers are removed; counts should be ~0 after cutadapt
  primer_counts_post_cut <- rbind(
    FWD.ForwardReads = sapply(FWD.orients, primer_hits, fastq_file = fnFs.cut[[1]]),
    FWD.ReverseReads = sapply(FWD.orients, primer_hits, fastq_file = fnRs.cut[[1]]),
    REV.ForwardReads = sapply(REV.orients, primer_hits, fastq_file = fnFs.cut[[1]]),
    REV.ReverseReads = sapply(REV.orients, primer_hits, fastq_file = fnRs.cut[[1]])
  )
  cat("\nPrimer hit counts after cutadapt (should be ~0):\n")
  print(primer_counts_post_cut)

  # ------------------------------------------------------------------
  # Step B4: Quality filter cutadapt-trimmed files
  # ------------------------------------------------------------------
  # This filterAndTrim applies truncation and quality thresholds to the
  # primer-free reads from cutadapt. It is a SEPARATE call from the N-filter
  # above and writes to a different directory (filtered_dir, not n_filtered_dir).
  # Do not combine these two filterAndTrim calls - N-filtering must complete
  # before cutadapt runs, and quality filtering must follow cutadapt output.


  filtering_stats <- filterAndTrim(
    fnFs.cut, filtered_forward,   # Input: cutadapt-trimmed reads; output: quality-filtered reads
    fnRs.cut, filtered_reverse,,
    maxN        = 0,
    maxEE       = c(2, 5),
    truncQ      = 2,
    minLen      = 50,          # Discard reads that are too short after trimming; adjust as needed
    rm.phix     = TRUE,
    compress    = TRUE,
    multithread = use_multithread
  )

  # Hand off to shared downstream steps
  reads_for_dada_F <- filtered_forward
  reads_for_dada_R <- filtered_reverse
} # end if (pipeline_mode == "cutadapt")


# =============================================================================
# SHARED: FILTERING RESULTS
# (Both pipeline modes rejoin here via reads_for_dada_F / reads_for_dada_R)
# =============================================================================

cat("\nFiltering statistics - reads.in vs reads.out (first 6 samples):\n")
print(head(filtering_stats))
saveRDS(filtering_stats, file.path(filtered_dir, "filtering_merged_output.rds"))

# Check the first few filtered reads to confirm trimming behaved as expected
filtered_forward_fastq    <- readFastq(filtered_forward[[1]])
filtered_forward_summary  <- summarize_fastq(filtered_forward_fastq)
cat("\nFirst 25bp of filtered forward reads (primers should be absent):\n")
print(filtered_forward_summary$sequence_starts)
cat("\nFiltered forward read length distribution:\n")
print(filtered_forward_summary$length_distribution)


# =============================================================================
# SHARED: ERROR LEARNING
# =============================================================================
# DADA2's core innovation: learns the error profile from the data.
# Why this matters: Distinguishes true biological variants from sequencing errors.
# Traditional OTU methods cluster at 97% similarity (arbitrary threshold).
# DADA2 uses a statistical error model to infer exact sequences at
# single-nucleotide resolution.
#
# Important: learnErrors uses a random subset of reads to build the model.
# Always rerun this step if you change filtering settings above.

cat("\nLearning error rates from filtered reads...\n")
cat("(This step examines many reads to build the error model - be patient)\n\n")

err_forward <- learnErrors(reads_for_dada_F, multithread = use_multithread)
err_reverse  <- learnErrors(reads_for_dada_R, multithread = use_multithread)

# Visual sanity check:
#   Black dots = observed error rates per quality score
#   Black line = estimated error rates (should fit the dots reasonably well)
#   Red line   = expected error rates based purely on Q-scores
# A good fit confirms the error model has learned the run's characteristics.
# If dots diverge badly from the line, the error model is unreliable - check input read quality.
print(plotErrors(err_forward, nominalQ = TRUE))
print(plotErrors(err_reverse,  nominalQ = TRUE))
saveRDS(plotErrors(err_forward, nominalQ = TRUE), file.path(filtered_dir, "err_forward.rds"))
saveRDS(plotErrors(err_reverse,  nominalQ = TRUE), file.path(filtered_dir, "err_reverse.rds"))

# =============================================================================
# SHARED: DADA2 DENOISING
# =============================================================================
# Infer true Amplicon Sequence Variants (ASVs) by correcting sequencing errors.
#
# pool = "pseudo": Compromise between speed and sensitivity
#   "FALSE"  - Process each sample independently (fast, less sensitive to rare variants)
#   "TRUE"   - Pool all samples together (slow, most sensitive; not feasible at large scale)
#   "pseudo" - Two-pass partial pooling: first pass identifies rare variants globally,
#              second pass denoises each sample; recommended balance for most datasets
#
# For reproductive microbiome work where clinically relevant taxa (e.g., Ureaplasma)
# may appear rarely, pseudo-pooling is the correct default.

cat("\nRunning DADA2 algorithm on forward reads...\n")
dada_forward <- dada(
  reads_for_dada_F,
  err         = err_forward,
  pool        = "pseudo",
  multithread = use_multithread
)

cat("\nRunning DADA2 algorithm on reverse reads...\n")
dada_reverse <- dada(
  reads_for_dada_R,
  err         = err_reverse,
  pool        = "pseudo",
  multithread = use_multithread
)

# Inspect dada-class object for the first sample
cat("\nDada-class object for first sample (forward reads):\n")
print(dada_forward[[1]])


# =============================================================================
# SHARED: MERGE PAIRED READS
# =============================================================================
# Combine forward and reverse reads to reconstruct the full amplicon.
# Requires a minimum overlap between the two reads (typically >= 20 bp).
#
# If merge rates are very low (< 20%), diagnose using the track table below and check:
#   1. Are primers fully removed? Residual primers prevent merging.
#   2. Is truncLen leaving enough overlap?
#      Overlap = truncLenF + truncLenR - expected_amplicon_length
#      For V4 (~253 bp): 150 + 150 - 253 = 47 bp; this should be sufficient.
#   3. Was data downloaded from ENA (not NCBI SRA)?
#      NCBI SRA's fastq-dump can corrupt read orientation for some datasets.
#      ENA downloads preserve the original submitted orientation.

mergers <- mergePairs(
  dada_forward, reads_for_dada_F,
  dada_reverse, reads_for_dada_R,
  minOverlap = 7,     # reduced from default 12; matches dataset ~8bp theoretical overlap, treat results skeptically
  maxMismatch  = 0,       # keep strict on mismatches to compensate for short overlap
  verbose = TRUE
)
print(head(mergers[[1]]))

# =============================================================================
# SHARED: BUILD SEQUENCE TABLE AND LENGTH FILTER
# =============================================================================

sequence_table <- makeSequenceTable(mergers)

cat("\nASV length distribution (raw, before length filtering):\n")
print(table(nchar(getSequences(sequence_table))))

# Keep only sequences matching the expected V4 amplicon length.
# V4 (515F-806R) merged amplicons are ~253 bp; accept a narrow range to
# exclude non-specific products or chimera fragments outside this range.
# Update expected_v4_lengths if you are using a different region or primer pair.
expected_v4_lengths <- 250:256
v4_sequence_table   <- sequence_table[, nchar(colnames(sequence_table)) %in% expected_v4_lengths]

cat("\nASVs retained after length filter (", paste(range(expected_v4_lengths), collapse = "-"), "bp):",
    ncol(v4_sequence_table), "\n")


# =============================================================================
# SHARED: REMOVE CHIMERAS
# =============================================================================
# Remove chimeric sequences that likely arose during PCR amplification.
# This step is critical for accurate diversity estimates.
#
# method = "consensus": A chimera must be flagged as bimeric in a majority
# of samples to be removed; more conservative than "pooled", less so than "per-sample".
#
# If the fraction retained (printed below) is < 0.70, check that primers were
# fully removed before constructing the sequence table - residual primers are a
# common cause of inflated chimera rates.

nonchim_sequence_table <- removeBimeraDenovo(
  v4_sequence_table,
  method      = "consensus",
  multithread = use_multithread,
  verbose     = TRUE
)

cat("\nFraction of reads retained after chimera removal (expect 0.70-0.95):\n")
print(sum(nonchim_sequence_table) / sum(v4_sequence_table))

cat("\nFinal ASV table dimensions (samples x ASVs):\n")
print(dim(nonchim_sequence_table))

saveRDS(nonchim_sequence_table, file.path(filtered_dir, "ASV_table.rds"))


# =============================================================================
# SHARED: READ TRACKING TABLE
# =============================================================================
# Summary of how many reads were retained at each pipeline step, per sample.
# This is your primary diagnostic tool.
#
# How to read the track table:
#   input     -> filtered  : reads lost to quality filtering (should be 10-30%)
#   filtered  -> denoisedF : small loss expected (denoising consolidates errors)
#   denoisedF -> merged    : if this cliff is large (>50% lost), primers or overlap problem
#   merged    -> nonchim   : large loss here means high chimera rate; check primer removal
#
# Column count_unique_reads: counts distinct sequences (uniques), not total reads

count_unique_reads <- function(x) sum(getUniques(x))
# Guard: filtering may drop some samples entirely if they had very few reads.
# Use the rownames of nonchim_sequence_table as the authoritative sample list.
surviving_samples <- rownames(nonchim_sequence_table)

#Strip the "_F_cut.fastq.gz" suffix from the rownames of filtering_stats to match sample names
clean_names <- sub("_F_cut.fastq.gz", "", rownames(filtering_stats))
rownames(filtering_stats) <- clean_names

# Sanity check: Take the intersection of surviving samples and filtering_stats to ensure we only include samples that made it through filtering
common_samples <- intersect(surviving_samples, rownames(filtering_stats))

tracking_table <- cbind(
  filtering_stats[surviving_samples, , drop = FALSE],
  sapply(dada_forward[surviving_samples], count_unique_reads),
  sapply(dada_reverse[surviving_samples], count_unique_reads),
  sapply(mergers[surviving_samples],      count_unique_reads),
  rowSums(nonchim_sequence_table)
)
colnames(tracking_table) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
cat("\nRead tracking table:\n")
print(tracking_table)

saveRDS(tracking_table, file.path(filtered_dir, "tracking_table.rds"))


# =============================================================================
# SHARED: ASSIGN TAXONOMY
# =============================================================================
# Assign taxonomy using the SILVA reference database.
# k-mer based matching: fast but approximate at species level.
# Species-level assignment is not guaranteed; many ASVs will stop at genus or family.
# For species-level resolution, use addSpecies() with silva_species_assignment_v138.1.fa.gz
# (separate download from the same Zenodo record).

taxonomy <- assignTaxonomy(
  nonchim_sequence_table,
  reference_fasta,
  multithread = use_multithread
)

taxonomy_preview <- taxonomy
rownames(taxonomy_preview) <- NULL
cat("\nTaxonomy preview (first 6 ASVs):\n")
print(head(taxonomy_preview))

# Save taxonomy - required for phyloseq and all downstream analysis
saveRDS(taxonomy, file.path(filtered_dir, "taxonomy.rds"))


# =============================================================================
# PIPELINE COMPLETE
# =============================================================================
cat("\n")
cat("============================================\n")
cat("   PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("   Mode:", pipeline_mode, "\n")
cat("============================================\n")
cat("\nKey outputs saved in:", filtered_dir, "\n")
cat("  1. ASV_table.rds               - abundance of each ASV per sample\n")
cat("  2. taxonomy.rds                - taxonomic assignments (Kingdom to Species)\n")
cat("  3. filtering_merged_output.rds - quality filtering statistics\n")
cat("  4. tracking_table.rds          - read counts at each pipeline step\n")
cat("\nNext steps:\n")
cat("  - Load ASV table and taxonomy into phyloseq for diversity analysis\n")
cat("  - Create bar plots of taxonomic composition\n")
cat("  - Calculate alpha/beta diversity metrics\n")
cat("  - Test for differential abundance between sample groups (DESeq2)\n")
cat("\n")

# =============================================================================
# END OF PIPELINE
# =============================================================================
