# Purpose:  One-time setup to install all required R packages
#
# Run this script ONCE when setting up the project on a new computer
# It automatically installs only the packages you don't already have
# =============================================================================

cat("=== Microbiome Analysis Package Installer ===\n\n")
cat("This script will install required packages from:\n")
cat("  - Bioconductor (bioinformatics packages)\n")
cat("  - CRAN (general R packages)\n\n")

# ---- HELPER FUNCTION ----
# Checks which packages are missing and installs only those
# Prevents reinstalling packages you already have
install_missing <- function(packages, installer) {
  # Check each package - is it already installed?
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]

  if (length(missing) > 0L) {
    cat("Installing:", paste(missing, collapse = ", "), "\n")
    installer(missing)
  } else {
    cat("All packages already installed!\n")
  }
}

# ---- INSTALL BIOCMANAGER ----
# BiocManager is required to install Bioconductor packages
# Bioconductor = repository for bioinformatics/genomics R packages

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("Installing BiocManager (required for Bioconductor packages)...\n")
  install.packages("BiocManager")
} else {
  cat("BiocManager already installed.\n")
}

# ---- BIOCONDUCTOR PACKAGES ----
cat("\n--- Checking Bioconductor packages ---\n")

bioc_packages <- c(
  "dada2", # Core: 16S amplicon analysis pipeline
  "phyloseq", # Microbiome data visualization and analysis
  "DECIPHER", # DNA sequence alignment and tree building
  "Biostrings", # DNA/RNA/protein sequence manipulation
  "ShortRead" # FASTQ file reading and quality control
)

# What each package does:
# - dada2: Infers ASVs from sequencing data, assigns taxonomy
# - phyloseq: Combines ASV table + taxonomy + metadata for downstream analysis
# - DECIPHER: Aligns 16S sequences, builds phylogenetic trees
# - Biostrings: Handles DNA sequence objects in R
# - ShortRead: Reads FASTQ files, calculates quality metrics

install_missing(bioc_packages, function(pkgs) {
  BiocManager::install(pkgs, ask = FALSE, update = FALSE)
})

# ---- CRAN PACKAGES ----
cat("\n--- Checking CRAN packages ---\n")

cran_packages <- c(
  "ggplot2", # Plotting and data visualization
  "tidyverse", # Data manipulation (dplyr, tidyr, etc.)
  "vegan", # Ecology statistics (diversity indices, ordination)
  "ape", # Phylogenetic tree manipulation and plotting
  "pheatmap", # Heatmap visualization
  "here" # Portable file paths (works on any computer)
)

install_missing(cran_packages, function(pkgs) {
  install.packages(pkgs)
})

# ---- VERIFICATION ----
cat("\n=== Installation Complete! ===\n\n")
cat("Verifying core packages can be loaded...\n")

# Test that critical packages load without errors
test_packages <- c("dada2", "phyloseq", "ggplot2", "here")
test_results <- sapply(test_packages, function(pkg) {
  result <- tryCatch(
    {
      library(pkg, character.only = TRUE, quietly = TRUE)
      "✓ OK"
    },
    error = function(e) {
      "✗ FAILED"
    }
  )
  result
})

cat("\nPackage verification:\n")
for (i in seq_along(test_results)) {
  cat("  ", names(test_results)[i], ": ", test_results[i], "\n", sep = "")
}

# ---- NEXT STEPS ----
if (all(test_results == "✓ OK")) {
  cat("\n✓ All packages installed successfully!\n")
  cat("\nNext steps:\n")
  cat("  1. Download SILVA database (see README.md)\n")
  cat("  2. Download raw FASTQ files (or use provided test data)\n")
  cat("  3. Run: scripts/pipelines/dada2_pipeline.R\n")
} else {
  cat("\n⚠ Some packages failed to load.\n")
  cat("Try manually installing failed packages, then re-run this script.\n")
}

# =============================================================================
# TROUBLESHOOTING
# =============================================================================
# If installation fails:
#
# 1. Update R to latest version (≥4.0)
#    - Download from: https://cran.r-project.org/
#
# 2. Update RStudio to latest version
#    - Download from: https://posit.co/download/rstudio-desktop/
#
# 3. For Bioconductor package errors:
#    BiocManager::install("package_name", force = TRUE)
#
# 4. For compilation errors on Mac/Linux:
#    - Install Xcode command line tools (Mac)
#    - Install build-essential (Linux: sudo apt-get install build-essential)
#
# 5. For Windows compilation errors:
#    - Install Rtools: https://cran.r-project.org/bin/windows/Rtools/
#
# 6. Check Bioconductor version compatibility:
#    BiocManager::version()  # Should match your R version
# =============================================================================
