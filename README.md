# 16S rRNA Microbiome Analysis Pipeline
### Reproductive Microbiome & IVF Outcomes · DADA2 · R · SILVA v138.1

**Author:** Olivia Williams  
**Dataset:** Okwelogu et al. (2021), PRJNA762524 — publicly available, analysis my own  
**Stack:** R, DADA2, cutadapt, SILVA v138.1, phyloseq  

---

## The Research Question

Infertility affects roughly one in six people globally. When couples pursue IVF, clinical attention typically focuses on embryo quality and hormonal factors — but the microbial environment of the reproductive tract may also shape outcomes in ways that are not yet fully understood.

Okwelogu et al. (2021) asked: **do the microbiome compositions of infertile couples' reproductive tracts correlate with IVF success?** Their paper sequenced samples from both partners using 16S rRNA amplicon sequencing, producing a publicly available dataset (PRJNA762524) that sits at the intersection of clinical microbiology, sequencing technology, and biostatistics.

That clinical framing is what made this dataset worth working with. Microbiome research is moving steadily into reproductive medicine, oncology, and pharmacomicrobiomics — the study of how gut and mucosal microbiomes interact with drug metabolism and therapeutic response. Building a pipeline on a clinically grounded dataset felt more useful than a synthetic tutorial.

---

## What This Project Demonstrates

This is not a tutorial replication. It is a from-scratch implementation of a full amplicon sequencing pipeline, with documented decisions, real diagnostic failures, and iterative fixes.

**Technical skills covered:**
- Raw FASTQ ingestion and quality assessment
- Primer/adapter removal with cutadapt
- Quality filtering, error modelling, and ASV denoising via DADA2
- Paired-end read merging
- Chimera detection and removal
- Taxonomic assignment against SILVA v138.1
- Read-tracking across every processing step
- Phylogenetic tree construction (NJ + ML) and visualisation
- Downstream diversity analysis: alpha (Shannon, Simpson, Hill numbers) and beta (PCoA, NMDS, PERMANOVA)

---

## The Problem I Had to Solve

This is the part of the project I am most proud of documenting.

The dataset uses 2×150 bp paired-end sequencing on the V4 16S region. My first paired-end attempts failed at the merging step: reads were producing only ~8 bp of overlap, well below the ~20 bp minimum DADA2 needs for reliable merging. No merged reads meant no usable ASVs from the paired-end pipeline.

Rather than switching to single-end and moving on, I diagnosed the issue from first principles:

1. Examined per-base quality profiles to confirm read quality was not the problem
2. Checked ASV length distributions to identify where reads were being truncated
3. Identified that `filterAndTrim()` with `trimLeft` was not reliably removing primers before truncation, leaving artefactual bases that compressed the overlap region
4. Switched to cutadapt for strict, discard-untrimmed primer removal before DADA2 processing

After applying cutadapt upstream, paired-end merging became viable and produced **896 ASVs** across the full dataset. For comparison, the single-end pipeline (which requires no merging) produced **268 ASVs** — a useful sanity check, since paired-end typically recovers more diversity.

The troubleshooting scripts and diagnostic plots are committed to the repository. The failure is documented, not hidden.

---

## Pipeline Modes

| Mode | Primer Removal | ASVs (full dataset) | Use Case |
|------|---------------|---------------------|----------|
| Paired-end (cutadapt) | cutadapt pre-step | **896** | Default — recommended |
| Paired-end (standard) | DADA2 `trimLeft` | — | When primer placement is consistent |
| Single-end | DADA2 `trimLeft` | **268** | Fallback when merging is not viable |

---

## Key Parameters

```r
filterAndTrim(
  trimLeft  = c(19, 21),   # Primer lengths for V4 region — dataset-specific
  truncLen  = c(150, 130), # Truncation at quality drop-off
  maxEE     = c(2, 5)      # Maximum expected errors per read
)
```

These are not defaults. They reflect decisions made after examining the quality profiles of this specific dataset. Applying them to new data without re-examining the profiles would produce incorrect results.

---

## 🔬 Pipeline Workflow

The full pipeline runs in nine steps. Each script prints a progress banner to the console at startup so you always know what stage is running.

### Step 0: Environment Setup

```r
source("scripts/setup_renv.R")
```

**Purpose:** Installs renv, installs all required R packages (Bioconductor + CRAN), and writes `renv.lock`. Run once on any new machine. On subsequent clones, use `renv::restore()` instead.

### Step 1: Download Sequences

```bash
# PRJNA762524 — the dataset used here
bash scripts/os/ena-file-download-read_run-PRJNA762524-fastq_ftp-*.sh

# Or place your own paired FASTQ files in data/raw/fastq/
# Naming convention: SAMPLE_1.fastq.gz / SAMPLE_2.fastq.gz
```

**Purpose:** Downloads all 72 paired-end FASTQ files from ENA FTP into `data/raw/fastq/`. The download script was generated directly from the ENA portal for PRJNA762524.

### Step 2: Read Statistics (QC)

```bash
bash scripts/os/seqkit_stats.sh
```

**Purpose:** Summarises read count, length distribution, GC content, and N50 for every raw FASTQ file before any processing. Output: `outputs/qc/seqkit_stats.tsv`. Requires `seqkit` on PATH.

### Step 3: Fetch and Process Metadata

```bash
bash scripts/os/fetch_metadata.sh
```

```r
Rscript scripts/process_metadata.R
```

**Purpose:** `fetch_metadata.sh` extracts SRR accessions from the download script, fetches run-level metadata from NCBI via `ffq`, and saves the raw JSON to `outputs/qc/`. `process_metadata.R` then parses that JSON into a tidy CSV with SRR accessions as row names, matching the FASTQ filenames. Output: `data/processed/dada2/metadata.csv`. Requires `ffq` (`pip install ffq`).

### Step 4: DADA2 Pipeline

```r
# Open scripts/dada2_pipeline.R
# Set pipeline_mode <- "cutadapt"  (default, recommended)
# Run from top to bottom in RStudio, or:
Rscript scripts/dada2_pipeline.R
```

**Purpose:** The core pipeline. Quality profiles → primer removal (cutadapt via WSL) → N-filtering → quality filtering → error learning → DADA2 denoising → paired-end merging → chimera removal → taxonomic assignment (SILVA v138.1). Quality and error-model plots are saved to `outputs/figures/` as PDFs. All key outputs are written to `data/processed/dada2/` and mirrored to `outputs/`.

**Expected runtime:** 10–30 min for a small sample set; longer on first run (error learning and taxonomy are compute-intensive).

### Step 5: Phylogenetic Tree

```r
Rscript scripts/phylogeny.R
```

**Purpose:** Aligns ASV sequences (DECIPHER), builds a neighbor-joining starting tree, and refines it with maximum likelihood optimisation (GTR+G+I model via phangorn). Saves the tree in both RDS and Newick formats. Output: `outputs/phylogeny/`.

**Expected runtime:** 30–90 min depending on ASV count and CPU.

### Step 6: Alpha & Beta Diversity Analysis

```r
Rscript scripts/alpha_beta_analysis.R
```

**Purpose:** Computes alpha diversity (Shannon, Simpson, Chao1, Hill numbers q=0,1,2, Faith's PD), performs PERMANOVA (Bray-Curtis, UniFrac), betadisper, and generates rarefaction curves. Outputs metrics tables and statistical results to `outputs/alpha_beta/`.

### Step 7: Alpha & Beta Diversity Visualisation

```r
Rscript scripts/alpha_beta_visualisation.R
```

**Purpose:** Creates publication-quality figures from the analysis outputs: alpha diversity boxplots by group, taxonomic barplots (phylum and genus level), PCoA and NMDS ordination plots, and a ComplexHeatmap of top ASVs by sample. Output: `outputs/figures/`.

### Step 8: Phylogenetic Visualisation

```r
Rscript scripts/phylo_visualization.R
```

**Purpose:** Produces tree figures using ggtree: full circular ASV tree, taxonomy-collapsed trees (family and genus level), top-50 ASV subtree, and a phylogenetic heatmap. Output: `outputs/figures/`.

---

## Reproducing This Analysis

**Requirements:** R ≥ 4.0, cutadapt (installed in WSL on Windows)

This project uses [renv](https://rstudio.github.io/renv/) to lock every R package to an exact version, so the environment is fully reproducible across machines.

### First-time setup (new machine or new collaborator)

```r
source("scripts/setup_renv.R")
```

R will print `"Project '...' loaded. [renv x.y.z]"` on every subsequent startup to confirm the isolated library is active.

### Restoring from a committed `renv.lock` (subsequent clones)

```r
renv::restore()
```

This reads `renv.lock` and installs every package at its locked version — no manual package hunting required.

### After pulling an update that added packages

```r
renv::restore()
```

Same command — syncs your local library to the updated lock file.

### Packages installed

| Source       | Packages                                                                                                                    |
|:-------------|:----------------------------------------------------------------------------------------------------------------------------|
| Bioconductor | `dada2`, `phyloseq`, `DECIPHER`, `Biostrings`, `ShortRead`, `phangorn`, `ggtree`, `ComplexHeatmap`, `microbiome`            |
| CRAN         | `ggplot2`, `tidyverse`, `vegan`, `ape`, `pheatmap`, `here`, `plotly`, `htmlwidgets`, `RColorBrewer`, `picante`              |

Note on FigTree: FigTree (<http://tree.bio.ed.ac.uk/software/figtree/>) is a Java-based standalone tree viewer useful for exploring Newick files; download and run FigTree separately to open files in `outputs/phylogeny/`.

**Reference database:** SILVA v138.1 training sets from [Zenodo](https://zenodo.org/record/4587955) → place in `data/external/reference/`

All scripts use `here::here()` for portable paths. Clone the repo, place inputs in the specified directories, and run the pipeline in order.

---

## Outputs

| File | Description |
| ---- | ----------- |
| `data/processed/dada2/ASV_paired_end_table.rds` | Sample × ASV count matrix (896 ASVs, paired-end) |
| `data/processed/dada2/ASV_table_SE.rds` | Sample × ASV count matrix (268 ASVs, single-end) |
| `data/processed/dada2/taxonomy.rds` | Kingdom → Species assignments via SILVA v138.1 |
| `data/processed/dada2/tracking_table.rds` | Read counts at each pipeline stage, per sample |
| `data/processed/dada2/metadata.csv` | Sample metadata from NCBI SRA (SRR accessions as row names) |
| `outputs/figures/quality_profiles_raw_*.pdf` | Raw read quality profiles (per-base Q scores, forward & reverse) |
| `outputs/figures/error_model_*.pdf` | DADA2 error model diagnostics (observed vs. estimated error rates) |
| `outputs/qc/seqkit_stats.tsv` | Per-file read statistics: count, length, GC%, N50 |
| `outputs/qc/ffq_prjna762524.json` | Raw project metadata from NCBI via ffq |
| `outputs/phylogeny/asv_tree_ml.{rds,newick}` | Maximum-likelihood phylogenetic tree |
| `outputs/phylogeny/asv_tree_nj.{rds,newick}` | Neighbor-joining tree (starting point for ML) |
| `outputs/alpha_beta/alpha_metrics.csv` | Shannon, Simpson, Chao1, Hill numbers per sample |
| `outputs/alpha_beta/permanova_bray.rds` | PERMANOVA results (Bray-Curtis distance) |
| `outputs/alpha_beta/bray_distance.rds` | Bray-Curtis distance matrix |
| `outputs/figures/alpha_diversity_summary_boxplots.png` | Alpha diversity by group |
| `outputs/figures/beta_pcoa_bray_curtis.png` | PCoA ordination (Bray-Curtis) |
| `outputs/figures/beta_nmds_bray_curtis.png` | NMDS ordination (Bray-Curtis) |
| `outputs/figures/taxonomic_barplot_*.png` | Relative abundance at phylum and genus level |
| `outputs/figures/full_asv_tree_circular.png` | Circular phylogenetic tree of all 896 ASVs |

---

## Scripts Reference

| Script | Language | Purpose |
| ------ | -------- | ------- |
| `scripts/setup_renv.R` | R | One-time renv environment initialisation |
| `scripts/os/ena-file-download-*.sh` | Bash | Bulk FASTQ download from ENA FTP |
| `scripts/os/seqkit_stats.sh` | Bash | Per-file read statistics (seqkit) |
| `scripts/os/fetch_metadata.sh` | Bash | NCBI metadata download via ffq |
| `scripts/process_metadata.R` | R | Parses ffq JSON → metadata CSV |
| `scripts/dada2_pipeline.R` | R | Core DADA2 pipeline (paired-end, with cutadapt mode) |
| `scripts/phylogeny.R` | R | Phylogenetic tree construction (NJ + ML) |
| `scripts/alpha_beta_analysis.R` | R | Diversity metrics, PERMANOVA, rarefaction curves |
| `scripts/alpha_beta_visualisation.R` | R | Publication figures (boxplots, PCoA, NMDS, heatmap) |
| `scripts/phylo_visualization.R` | R | Tree figures (circular, collapsed, fan, phylogenetic heatmap) |
| `scripts/load_packages.R` | R | Manual package installation helper (alternative to setup_renv.R) |
| `scripts/troubleshooting.R` | R | Diagnostic checks and debugging — run when pipeline fails |
| `scripts/single end 16s.R` | R | Single-end pipeline variant (fallback when paired-end merging fails) |

---

## What's Next

1. ~~**Alpha and beta diversity analysis** — Shannon/Simpson indices, PCoA ordination, PERMANOVA, rarefaction curves~~ ✓ (`scripts/alpha_beta_analysis.R`)
2. ~~**Phylogenetic tree construction and visualisation**~~ ✓ (`scripts/phylogeny.R`, `scripts/phylo_visualization.R`)
3. ~~**Taxonomic barplots and diversity heatmaps**~~ ✓ (`scripts/alpha_beta_visualisation.R`)
4. **Differential abundance testing** — DESeq2 and LEfSe adapted for microbiome count data
5. **Metadata integration** — Correlating microbiome composition with the clinical IVF outcome variables reported in Okwelogu et al.
6. **Co-occurrence networks** — Identifying microbial taxa that co-occur or are mutually exclusive across samples
7. **Additional downstream analyses** — as the project continues to grow

---

## Downstream Analyses

The ASV and taxonomy tables produced here feed directly into:

- **Alpha/beta diversity** (`scripts/alpha_beta_analysis.R`) — Shannon/Simpson indices, PCoA ordination, PERMANOVA
- **Differential abundance** — DESeq2 or LEfSe, for identifying taxa that differ significantly between clinical groups
- **Metadata integration** — correlating microbiome composition with the IVF outcome variables reported in Okwelogu et al.
- **Phylogenetic analysis** (`scripts/phylogeny.R`, `scripts/phylo_visualization.R`) — tree construction and visualisation

---

## Pharma & Clinical Context

Microbiome research is increasingly relevant to drug development and clinical trial design. Key areas where this type of pipeline applies:

**Reproductive medicine** — vaginal and endometrial microbiome profiles are being explored as predictors of IVF success and implantation failure. This dataset sits directly in that space.

**Biomarker discovery** — differential abundance testing on ASV tables is a standard approach for identifying microbial signatures associated with treatment response, disease state, or clinical outcome.

**Pharmacomicrobiomics** — the study of how gut and mucosal microbiomes interact with drug metabolism and therapeutic response is an emerging field where amplicon sequencing pipelines like this one provide the foundational data.

The pipeline here produces the foundational data objects (ASV table, taxonomy, phylogenetic tree) required for all of the above.

---

## References

- Okwelogu et al. (2021). Microbiome Compositions From Infertile Couples Seeking In Vitro Fertilization. *Frontiers in Cellular and Infection Microbiology*, 11:709372. [https://doi.org/10.3389/fcimb.2021.709372](https://doi.org/10.3389/fcimb.2021.709372)
- Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- Quast et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- McMurdie & Holmes (2013). phyloseq: An R package for reproducible interactive analysis of microbiome census data. *PLOS ONE*. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)

---

*Data sourced from NCBI SRA (public domain). SILVA database used under academic licence.*
