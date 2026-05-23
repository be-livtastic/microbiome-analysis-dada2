# Project Structure

This project is organized around a simple reproducible layout:

- `data/raw/fastq/` for unmodified FASTQ files
- `data/external/reference/` for SILVA reference databases and other third-party reference files
- `data/processed/dada2/` for filtering outputs, ASV tables, and taxonomy assignments
- `scripts/` for analysis, setup, and troubleshooting scripts
- `docs/` for project notes and workflow guidance
- `outputs/figures/` for plots and diagnostic images
- `outputs/tables/` for exported tables and summaries
- `outputs/qc/` for quality control outputs (seqkit read statistics, ffq metadata JSON)

Current file mapping:

- `raw_data/*.fastq` should live in `data/raw/fastq/`
- `raw_data/filtered/` should live in `data/processed/dada2/`
- `silva_*.fa.gz` should live in `data/external/reference/`
- `Rplot*.png` should live in `outputs/figures/`
- `scripts/setup/load_packages.R` is the package setup helper
- `scripts/download/ena-file-download-read_run-PRJNA762524-fastq_ftp-20260407-1506.sh` is the ENA download helper
- `scripts/os/linux_script.sh` is the Linux-specific helper
- `dada2_pipeline.R`, `single end 16s.R`, and `troubleshooting.R` are good candidates for `scripts/pipelines/` and `scripts/diagnostics/`

Path convention:

- Prefer `here::here()` or project-relative paths instead of hard-coded absolute paths.
- Keep raw data read-only.
- Write all intermediate and final outputs under `data/processed/` or `outputs/`.
- If a new file is temporary, place it in `outputs/` or a clearly named subfolder so it is easy to clean up later.

Notes:

- FigTree is a popular Java-based tree viewer for Newick files. It's not an R package; download it from the FigTree project page and use it to open trees saved to `outputs/phylogeny/`.
