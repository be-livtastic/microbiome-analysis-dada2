library(dada2)
library(here)

filtered_dir <- here::here("data", "processed", "dada2")
outputs_dir  <- here::here("outputs")
dir.create(outputs_dir, recursive = TRUE, showWarnings = FALSE)

# Helper
safe_copy <- function(src, dst, overwrite = TRUE) {
  ok <- file.copy(src, dst, overwrite = overwrite)
  if (!all(ok)) stop(sprintf("Failed to copy %s to %s", src, dst))
  invisible(TRUE)
}

# Load the ASV table saved by the pipeline
cat("Loading nonchim sequence table...\n")
nonchim_sequence_table <- readRDS(file.path(filtered_dir, "ASV_table.rds"))
cat("ASV table dimensions:", dim(nonchim_sequence_table), "\n")

# CORRECT reference file for assignTaxonomy
reference_fasta <- here::here(
  "data", "external", "reference",
  "silva_nr99_v138.1_train_set.fa.gz"
)

species_fasta <- here::here(
  "data", "external", "reference",
  "silva_species_assignment_v138.1.fa.gz"
)

cat("Reference file:", reference_fasta, "\n")
cat("File exists:", file.exists(reference_fasta), "\n")

# Assign genus-level taxonomy
cat("\nRunning assignTaxonomy (this takes 20-60 min for 72 samples)...\n")
taxonomy <- assignTaxonomy(
  nonchim_sequence_table,
  reference_fasta,
  multithread = TRUE
)

# Add species-level assignments
cat("\nRunning addSpecies...\n")
taxonomy <- addSpecies(taxonomy, species_fasta)

taxonomy_preview <- taxonomy
rownames(taxonomy_preview) <- NULL
cat("\nTaxonomy preview (first 6 ASVs):\n")
print(head(taxonomy_preview))

saveRDS(taxonomy, file.path(filtered_dir, "taxonomy.rds"))
safe_copy(file.path(filtered_dir, "taxonomy.rds"), file.path(outputs_dir, "taxonomy.rds"))
cat("\nTaxonomy saved successfully.\n")

cat("\n============================================\n")
cat("   PIPELINE COMPLETED SUCCESSFULLY!\n")
cat("============================================\n")
cat("\nKey outputs in:", filtered_dir, "\n")
cat("  1. ASV_table.rds\n")
cat("  2. taxonomy.rds\n")
cat("  3. tracking_table.rds\n")
cat("\nNext steps:\n")
cat("  Rscript scripts/phylogeny.R\n")
cat("  Rscript scripts/alpha_beta_analysis.R\n")
