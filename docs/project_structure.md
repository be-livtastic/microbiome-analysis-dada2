# Project Structure

This document describes the layout, conventions, and execution order for the `microbiome-analysis-dada2` pipeline. It is intended as a reference for anyone running, extending, or troubleshooting the analysis.

---

## Directory Layout

```text
microbiome-analysis-dada2/
│
├── data/
│   ├── raw/
│   │   └── fastq/              Raw FASTQ files (not tracked in git — download per README)
│   │                           Naming: SRR*_1.fastq.gz (forward), SRR*_2.fastq.gz (reverse)
│   │
│   ├── external/
│   │   └── reference/          SILVA v138.1 training sets (not tracked — download from Zenodo)
│   │                           silva_nr99_v138.1_train_set.fa.gz
│   │                           silva_species_assignment_v138.1.fa.gz
│   │
│   └── processed/
│       ├── dada2/              Primary DADA2 outputs (filtered reads, ASV table, taxonomy)
│       ├── dada2_n_filtered/   N-filtered reads (intermediate, cutadapt mode only)
│       └── dada2_cutadapt/     Primer-trimmed reads before quality filtering (cutadapt mode)
│
├── scripts/
│   ├── dada2_pipeline.R        Core DADA2 pipeline — run this first after data download
│   ├── phylogeny.R             Phylogenetic tree construction (NJ + ML)
│   ├── alpha_beta_analysis.R   Alpha and beta diversity metrics and statistics
│   ├── alpha_beta_visualisation.R  Publication figures (boxplots, PCoA, NMDS, heatmap)
│   ├── phylo_visualization.R   Tree figures (circular, collapsed, fan, phylogenetic heatmap)
│   ├── process_metadata.R      Parses ffq JSON → metadata CSV
│   ├── setup_renv.R            One-time renv environment initialisation
│   ├── load_packages.R         Manual package installation helper
│   ├── troubleshooting.R       Diagnostic checks — run when something fails
│   ├── single end 16s.R        Single-end pipeline variant (fallback)
│   └── os/
│       ├── fetch_metadata.sh                                    Download NCBI metadata via ffq
│       ├── seqkit_stats.sh                                      Raw read statistics
│       └── ena-file-download-read_run-PRJNA762524-*.sh          Bulk FASTQ download from ENA
│
├── outputs/
│   ├── figures/                All PNG and PDF plots (quality profiles, diversity, trees)
│   ├── phylogeny/              Phylogenetic tree objects (RDS) and Newick files
│   ├── alpha_beta/             Diversity metrics CSV files and statistical results (RDS)
│   ├── qc/                     Quality control outputs (seqkit TSV, ffq metadata JSON)
│   └── single_end/             Single-end pipeline outputs (for comparison)
│
├── docs/
│   └── project_structure.md   This file
│
├── renv/                       R package library (managed by renv — not tracked in git)
├── renv.lock                   Locked package versions (tracked in git — essential for reproducibility)
├── .Rprofile                   Sources renv/activate.R on startup
└── README.md                   Full project documentation and pipeline instructions
```

---

## Execution Order

Run scripts in this order. Each script prints a banner to the console at startup so you always know what stage is running.

| Step | Script | Language | What It Does |
| ---- | ------ | -------- | ------------ |
| 0 | `scripts/setup_renv.R` | R | Install renv and all R packages; write `renv.lock` |
| 1 | `scripts/os/ena-file-download-*.sh` | Bash | Download paired FASTQ files from ENA FTP |
| 2 | `scripts/os/seqkit_stats.sh` | Bash | Summarise read counts, lengths, GC%, N50 |
| 3a | `scripts/os/fetch_metadata.sh` | Bash | Download NCBI run metadata via ffq |
| 3b | `scripts/process_metadata.R` | R | Parse ffq JSON → `data/processed/dada2/metadata.csv` |
| 4 | `scripts/dada2_pipeline.R` | R | Core pipeline: QC → primer removal → filter → denoise → merge → taxonomy |
| 5 | `scripts/phylogeny.R` | R | Align ASVs; build NJ and ML trees |
| 6 | `scripts/alpha_beta_analysis.R` | R | Compute diversity metrics, PERMANOVA, rarefaction curves |
| 7 | `scripts/alpha_beta_visualisation.R` | R | Generate publication figures |
| 8 | `scripts/phylo_visualization.R` | R | Generate phylogenetic tree figures |

---

## Key File Paths

### Inputs (must be downloaded manually)

| File | Location | Source |
| ---- | -------- | ------ |
| Raw FASTQ | `data/raw/fastq/SRR*_{1,2}.fastq.gz` | ENA FTP (Step 1 script) |
| SILVA genus reference | `data/external/reference/silva_nr99_v138.1_train_set.fa.gz` | Zenodo 4587955 |
| SILVA species reference | `data/external/reference/silva_species_assignment_v138.1.fa.gz` | Zenodo 4587955 |

### Primary Outputs (produced by dada2_pipeline.R)

| File | Description |
| ---- | ----------- |
| `data/processed/dada2/ASV_paired_end_table.rds` | Sample × ASV count matrix (896 ASVs) |
| `data/processed/dada2/taxonomy.rds` | Kingdom → Species taxonomy (SILVA v138.1) |
| `data/processed/dada2/tracking_table.rds` | Read counts at every pipeline stage |
| `data/processed/dada2/metadata.csv` | Sample metadata with SRR accessions as row names |
| `outputs/figures/quality_profiles_raw_forward.pdf` | Forward read quality profiles |
| `outputs/figures/quality_profiles_raw_reverse.pdf` | Reverse read quality profiles |
| `outputs/figures/error_model_forward.pdf` | Forward error model diagnostic |
| `outputs/figures/error_model_reverse.pdf` | Reverse error model diagnostic |

### Downstream Outputs

| File | Produced By |
| ---- | ----------- |
| `outputs/phylogeny/asv_tree_ml.{rds,newick}` | `phylogeny.R` |
| `outputs/phylogeny/asv_tree_nj.{rds,newick}` | `phylogeny.R` |
| `outputs/alpha_beta/alpha_metrics.csv` | `alpha_beta_analysis.R` |
| `outputs/alpha_beta/permanova_bray.rds` | `alpha_beta_analysis.R` |
| `outputs/alpha_beta/bray_distance.rds` | `alpha_beta_analysis.R` |
| `outputs/figures/alpha_diversity_summary_boxplots.png` | `alpha_beta_visualisation.R` |
| `outputs/figures/beta_pcoa_bray_curtis.png` | `alpha_beta_visualisation.R` |
| `outputs/figures/beta_nmds_bray_curtis.png` | `alpha_beta_visualisation.R` |
| `outputs/figures/taxonomic_barplot_*.png` | `alpha_beta_visualisation.R` |
| `outputs/figures/full_asv_tree_circular.png` | `phylo_visualization.R` |
| `outputs/figures/family_collapsed_tree.png` | `phylo_visualization.R` |
| `outputs/figures/genus_collapsed_tree.png` | `phylo_visualization.R` |

---

## Path Conventions

- Prefer `here::here()` over hard-coded absolute paths in all R scripts. This makes paths work regardless of where the script is called from.
- Raw data is read-only. Never write to `data/raw/`.
- All intermediate and final outputs go under `data/processed/` or `outputs/`.
- Plots from `dada2_pipeline.R` go to `outputs/figures/` as PDFs. All other figures go to `outputs/figures/` as PNGs.

---

## Environment

- **R packages:** Managed by [renv](https://rstudio.github.io/renv/). Run `source("scripts/setup_renv.R")` once, then `renv::restore()` on subsequent clones.
- **cutadapt:** Required for the default pipeline mode. On Windows, install inside WSL (`pip install cutadapt`) — the pipeline converts paths automatically.
- **seqkit:** Required for `seqkit_stats.sh`. Install from [bioinf.shenwei.me/seqkit](https://bioinf.shenwei.me/seqkit/).
- **ffq:** Required for `fetch_metadata.sh`. Install with `pip install ffq`.
- **FigTree:** Optional Java-based tree viewer for exploring Newick files. Download from the [FigTree project page](http://tree.bio.ed.ac.uk/software/figtree/) and open files in `outputs/phylogeny/` directly.

---

## Windows-Specific Notes

- `use_multithread <- TRUE` in `dada2_pipeline.R` works on Linux/macOS. On Windows, set this to `FALSE` if DADA2 steps hang or crash — Windows does not support forking.
- cutadapt is called via WSL. Ensure cutadapt is installed inside your WSL distribution: `wsl pip install cutadapt`.
- Path conversion from `C:/path/` to `/mnt/c/path/` is handled automatically inside `dada2_pipeline.R`.

---

## What Is and Is Not Tracked in Git

**Tracked:**

- All scripts in `scripts/`
- `data/processed/dada2/` (ASV tables, taxonomy, metadata — computed results)
- `outputs/` (figures, phylogeny, diversity tables)
- `renv.lock` (locked package versions)
- `.Rprofile`, `renv/activate.R`
- `docs/`, `README.md`

**Not tracked (too large or user-supplied):**

- `data/raw/fastq/` — users download raw FASTQ files per README instructions
- `data/external/reference/` — SILVA databases (users acquire from Zenodo)
- `renv/library/` — compiled package binaries (rebuilt from `renv.lock`)
- `.venv/` — Python virtual environment

---

## Planned Future Work

The pipeline is complete through phylogenetic visualisation and diversity analysis. Planned additions:

- Differential abundance testing (DESeq2, LEfSe)
- Metadata integration: correlating microbiome composition with IVF outcome variables
- Co-occurrence network analysis
- Additional downstream visualisations
