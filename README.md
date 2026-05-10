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
- Downstream diversity analysis (alpha/beta) and phylogenetic tree construction

**Pharma/clinical relevance:**
The pipeline produces analysis-ready microbiome objects suitable for differential abundance testing (e.g., DESeq2), which has direct application in microbiome-drug interaction research, clinical trial stratification by microbiome profile, and biomarker discovery in reproductive medicine.

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

## Repository Structure

```
microbiome_practice/
├── scripts/
│   ├── dada2_pipeline.R          # Main paired-end workflow (cutadapt mode default)
│   ├── single end 16s.R          # Single-end fallback pipeline
│   ├── troubleshooting.R         # Diagnostic scripts (overlap checks, quality plots)
│   ├── phylogeny.R               # Tree construction
│   ├── visualization.R           # Phylogenetic tree visualisation
│   ├── alpha_beta_analysis.R     # Diversity analysis
│   └── os/
│       ├── linux_script.sh       # NCBI SRA download (prefetch + fasterq-dump)
│       └── ena-file-download-*.sh # ENA FTP download
├── data/
│   ├── raw/fastq/                # Raw FASTQs (not tracked in Git)
│   ├── processed/dada2/          # Filtered reads, ASV table, taxonomy
│   └── external/reference/       # SILVA database (download separately)
├── outputs/
│   ├── figures/                  # Quality plots, diversity plots, tree visualisations
│   ├── phylogeny/                # Tree files (RDS, Newick)
│   └── tables/
└── docs/
    └── project_structure.md
```

---

## Outputs

| File | Description |
|------|-------------|
| `ASV_paired_end_table.rds` | Sample × ASV count matrix (896 ASVs, paired-end) |
| `ASV_table_SE.rds` | Sample × ASV count matrix (268 ASVs, single-end) |
| `taxonomy.rds` | Kingdom → Species assignments via SILVA v138.1 |
| Read tracking table | Per-sample read counts at every pipeline stage |
| Quality plots | Per-base profiles, error learning curves, ASV length distributions |
| Phylogenetic tree | Built from ASV sequences; viewable in FigTree or via `visualization.R` |

---

## Reproducing This Analysis

**Dependencies**

```r
# Bioconductor
dada2, ShortRead, Biostrings, phyloseq, ggtree

# CRAN
tidyverse, here, ggplot2, plotly, htmlwidgets
```

R ≥ 4.0 and cutadapt must be installed. On Windows, cutadapt runs via WSL.

**Reference database:** SILVA v138.1 training sets from [Zenodo](https://zenodo.org/record/4587955) — place in `data/external/reference/`

**Step 1 — Get the data**

```bash
# Option A: PRJNA762524 (the dataset used here)
bash scripts/os/ena-file-download-read_run-PRJNA762524-fastq_ftp-*.sh

# Option B: your own data
# Place paired FASTQ files in data/raw/fastq/
# Naming convention: SAMPLE_1.fastq / SAMPLE_2.fastq
```

**Step 2 — Run the pipeline**

```r
# Open scripts/dada2_pipeline.R
# Set pipeline_mode <- "cutadapt"  (default)
# Run from top to bottom, or step through interactively in RStudio
```

All scripts use `here::here()` for portable paths.

**Expected runtime:** 10–30 minutes for a small sample set; longer on first run (error learning is compute-intensive).

---

## Downstream Analyses

The ASV and taxonomy tables produced here feed directly into:

- **Alpha/beta diversity** (`scripts/alpha_beta_analysis.R`) — Shannon/Simpson indices, PCoA ordination, PERMANOVA
- **Differential abundance** — DESeq2 or LEfSe, for identifying taxa that differ significantly between clinical groups
- **Metadata integration** — correlating microbiome composition with the IVF outcome variables reported in Okwelogu et al.
- **Phylogenetic analysis** (`scripts/phylogeny.R`, `scripts/visualization.R`) — tree construction and visualisation

These next steps are in progress and will be added to this repository.

---

## Pharma & Clinical Context

Microbiome research is increasingly relevant to drug development and clinical trial design. Key areas where this type of pipeline applies:

**Pharmacomicrobiomics** — the gut microbiome modulates the metabolism of a significant proportion of clinically used drugs. Characterising microbiome composition in patient cohorts is a growing requirement in Phase I/II trial design.

**Reproductive medicine** — vaginal and endometrial microbiome profiles are being explored as predictors of IVF success and implantation failure. This dataset sits directly in that space.

**Biomarker discovery** — differential abundance testing on ASV tables is a standard approach for identifying microbial signatures associated with treatment response, disease state, or clinical outcome.

The pipeline here produces the foundational data objects (ASV table, taxonomy, phylogenetic tree) required for all of the above.

---

## References

- Okwelogu et al. (2021). Microbiome Compositions From Infertile Couples Seeking In Vitro Fertilization. *Frontiers in Cellular and Infection Microbiology*, 11:709372. [https://doi.org/10.3389/fcimb.2021.709372](https://doi.org/10.3389/fcimb.2021.709372)
- Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. [https://doi.org/10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- Quast et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*. [https://doi.org/10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219)
- McMurdie & Holmes (2013). phyloseq: An R package for reproducible interactive analysis of microbiome census data. *PLOS ONE*. [https://doi.org/10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)

---

*Data sourced from NCBI SRA (public domain). SILVA database used under academic licence.*
