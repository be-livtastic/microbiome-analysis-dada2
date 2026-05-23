# =============================================================================
# One-time renv environment setup
#
# Run this script ONCE on any new machine to reproduce the project environment.
# It will generate .Rprofile, renv/activate.R, and renv.lock automatically.
#
# After the first run, collaborators cloning the repo only need:
#   renv::restore()
# =============================================================================

cat("=== Microbiome Analysis — renv Environment Setup ===\n\n")

# ---- STEP 1: Install renv ----

if (!requireNamespace("renv", quietly = TRUE)) {
  cat("Installing renv...\n")
  install.packages("renv")
} else {
  cat("renv already installed.\n")
}

# ---- STEP 2: Initialise renv ----
# Creates .Rprofile (activates renv on R startup) and renv/activate.R
# bioconductor = TRUE ensures Bioconductor packages are handled correctly

cat("\nInitialising renv project...\n")
renv::init(bioconductor = TRUE)

# ---- STEP 3: Install all project packages ----
# renv::install() wraps both install.packages() and BiocManager::install()

cat("\nInstalling Bioconductor packages...\n")

bioc_packages <- c(
  "dada2",       # 16S amplicon denoising and ASV inference
  "phyloseq",    # Microbiome data integration and visualisation
  "DECIPHER",    # Sequence alignment for phylogenetic tree construction
  "Biostrings",  # DNA/RNA/protein sequence manipulation
  "ShortRead",   # FASTQ quality control
  "phangorn",    # Phylogenetic tree construction and analysis
  "ggtree"       # Phylogenetic tree visualisation
)

cat("\nInstalling CRAN packages...\n")

cran_packages <- c(
  "ggplot2",      # Plotting
  "tidyverse",    # Data manipulation (dplyr, tidyr, etc.)
  "vegan",        # Ecology statistics: diversity indices, ordination
  "ape",          # Phylogenetic tree manipulation
  "pheatmap",     # Heatmap visualisation
  "here",         # Portable project-relative file paths
  "plotly",       # Interactive HTML plots
  "htmlwidgets",  # Save interactive plots to HTML
  "RColorBrewer", # Colour palettes for plots
  "picante"       # Phylogenetic diversity metrics (optional: used if tree available)
)

renv::install(c(bioc_packages, cran_packages))

# ---- STEP 4: Snapshot installed versions into renv.lock ----
# renv.lock records the exact version of every package installed.
# Commit this file to git so others can reproduce the environment exactly.

cat("\nGenerating renv.lock...\n")
renv::snapshot()

# ---- DONE ----

cat("\n=== Setup complete! ===\n\n")
cat("The following files were generated and should be committed to git:\n")
cat("  .Rprofile        — activates renv automatically on R startup\n")
cat("  renv/activate.R  — renv bootstrap machinery\n")
cat("  renv.lock        — locked package versions\n\n")
cat("For collaborators cloning this repo, run:\n")
cat("  renv::restore()\n\n")
cat("On startup R will print:\n")
cat("  \"Project '...' loaded. [renv x.y.z]\"\n")
cat("This confirms the isolated project library is active.\n")
