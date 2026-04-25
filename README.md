# 16S rRNA Microbiome Analysis — IVF Reproductive Microbiome
## Practising with the Okwelogu et al. (2021) Dataset

**Author:** Olivia Williams  
**Stack:** R · DADA2 · cutadapt · SILVA v138.1  
**Domain:** Microbiome bioinformatics · 16S rRNA amplicon sequencing  
**Dataset:** Okwelogu et al. (2021). "Microbiome Compositions From Infertile Couples Seeking In Vitro Fertilization, Using 16S rRNA Gene Sequencing Methods: Any Correlation to Clinical Outcomes?" *Frontiers in Cellular and Infection Microbiology*, 11:709372. ([PubMed](https://pubmed.ncbi.nlm.nih.gov/34660337/) · [Full text](https://doi.org/10.3389/fcimb.2021.709372) · [PRJNA762524](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA762524))

---

## Overview

This project implements a complete **16S rRNA amplicon sequencing analysis pipeline** — from raw FASTQ files through to ASV inference, taxonomic classification, and quality-tracked outputs — using the DADA2 framework in R. The dataset is publicly available, but the scripts, decisions, and workflow are my own.

The biological question driving the original study: do microbiome compositions in the reproductive tracts of infertile couples correlate with IVF outcomes? That framing made this dataset particularly interesting to work with — it sits at the intersection of clinical microbiology, sequencing technology, and statistical biology.

---

## What I Built

### Pipeline
A reproducible, modular DADA2 workflow that handles:
- Raw data ingestion and quality assessment
- Adapter and primer removal via cutadapt
- Quality filtering, error modelling, and ASV denoising
- Paired-end merging
- Chimera removal
- Taxonomic assignment against SILVA v138.1
- Read-tracking across every processing step

### Outputs
| File | Description |
|---|---|
| `ASV_table.rds` | Sample × ASV count matrix |
| `taxonomy_SE.rds` | Kingdom → Species assignments (SILVA v138.1) |
| Quality plots | Per-base profiles, error learning curves, ASV length distributions |
| Read tracking table | Counts at each pipeline stage, per sample |

---

## Problem-Solving Record

The cleanest part of this project is not the final output — it's what happened before I got there.

The dataset uses 2×150 bp paired-end sequencing on the V4 16S region. My first attempt at paired-end merging failed: the reads produced only ~8 bp of overlap, far below the ~20 bp minimum DADA2 requires for reliable merging. I diagnosed the overlap issue directly from quality profiles and read length distributions, documented why merging was failing, and considered the tradeoffs of proceeding single-end versus resolving the upstream cause.

After revisiting the primer structure and running the raw reads through cutadapt for proper adapter and primer trimming, paired-end merging became viable. The revised paired-end pipeline is what this repository now runs.

**Pilot run:** 4 samples · 243 ASVs identified  
**Full dataset:** All samples processed; results uploaded to this repository

---

## Repository Structure

```
microbiome_practice/
├── scripts/
│   ├── pipelines/
│   │   ├── load_packages.R        # Environment setup
│   │   ├── dada2_pipeline.R       # Main paired-end workflow (cutadapt-trimmed input)
│   │   └── troubleshooting.R      # Diagnostic scripts (overlap checks, quality plots)
│   └── download/
│       └── ena-file-download-*.sh # SRA data retrieval
├── data/
│   ├── raw/fastq/                 # Raw FASTQ files (not tracked in Git)
│   ├── processed/dada2/           # Filtered reads, ASV table, taxonomy
│   └── external/reference/        # SILVA database (download separately)
├── outputs/
│   ├── figures/
│   └── tables/
└── docs/
    └── project_structure.md
```

---

## Key Parameters

```r
filterAndTrim(
  trimLeft  = c(19, 21),   # Primer lengths (V4 region, dataset-specific)
  truncLen  = c(150, 130), # Truncation at quality drop-off
  maxEE     = c(2, 5)      # Maximum expected errors per read
)
```

These are not defaults — they reflect deliberate decisions based on the quality profiles of this specific dataset. Adapting them to new data requires re-examining quality plots, not just copying the numbers.

---

## Reproducing This Analysis

**Dependencies:** R ≥ 4.0, cutadapt, and the following R packages:

```r
# Bioconductor
dada2, ShortRead, Biostrings, phyloseq

# CRAN  
tidyverse, here, ggplot2
```

**Reference database:** SILVA v138.1 training sets from [Zenodo](https://zenodo.org/record/4587955) → place in `data/external/reference/`

All scripts use `here::here()` for portable paths. Clone the repo, place inputs in the specified directories, and run `scripts/pipelines/dada2_pipeline.R`.

---

## What's Next

1. **Alpha and beta diversity analysis** — Shannon/Simpson indices, PCoA ordination, PERMANOVA
2. **Differential abundance testing** — DESeq2 and LEfSe adapted for microbiome count data
3. **Visualisation** — Taxonomic bar plots, co-occurrence networks, per-group composition heatmaps
4. **Metadata integration** — Correlating microbiome composition with the clinical IVF outcome variables reported in Okwelogu et al.

---

## References

- Callahan et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*. https://doi.org/10.1038/nmeth.3869  
- Quast et al. (2013). The SILVA ribosomal RNA gene database project. *Nucleic Acids Research*. https://doi.org/10.1093/nar/gks1219  
- McMurdie & Holmes (2013). phyloseq: An R package for reproducible interactive analysis of microbiome census data. *PLOS ONE*. https://doi.org/10.1371/journal.pone.0061217

---

*Data sourced from NCBI SRA (public domain). SILVA database used under academic licence.*


### In Case Anyone Wants to re-use this pipeline

**Option A: Use my test dataset (PRJNA762524)**
```bash
# Download from NCBI SRA using provided script
bash scripts/download/ena-file-download-read_run-PRJNA762524-fastq_ftp-*.sh
```

**Option B: Use your own 16S data**
- Place paired-end FASTQ files in `data/raw/fastq/`
- Files must follow naming convention: `SAMPLE_1.fastq`, `SAMPLE_2.fastq`

### Run the Pipeline

```r
# Open R project in RStudio
# Open scripts/pipelines/dada2_pipeline.R
# Run the entire script or step-through interactively

# Key parameters to adjust for your data:
# - trimLeft: primer length to remove (line 83)
# - truncLen: where to truncate reads based on quality (line 84)
# - expected_v4_lengths: amplicon size range (line 133)
```

**Expected runtime:** 10-30 minutes for 4 samples (longer for more samples or first run)

---

## 🔬 Pipeline Workflow

### Step 1: Quality Assessment
```r
plotQualityProfile(forward_reads)   # Visualize base quality scores
plotQualityProfile(reverse_reads)
```
**Decision point:** Determine where quality drops below Q30 → sets truncation length

### Step 2: Filter and Trim
```r
filterAndTrim(
  trimLeft = c(19, 21),   # Remove primers (dataset-specific)
  truncLen = c(150, 130), # Trim at quality drop-off
  maxEE = c(2, 5)         # Maximum expected errors
)
```
**Output:** Filtered FASTQ files with low-quality reads removed

### Step 3: Learn Error Rates
```r
errF <- learnErrors(filtFs, multithread = TRUE)
```
**Purpose:** DADA2 learns the sequencing error profile to distinguish true variants from errors

### Step 4: Denoise with DADA2
```r
dadaFs <- dada(filtFs, err = errF, pool = "pseudo")
```
**Purpose:** Infer true Amplicon Sequence Variants (ASVs) by correcting sequencing errors

### Step 5: Merge Paired Reads
```r
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs)
```
**Purpose:** Combine forward and reverse reads to reconstruct full amplicon

### Step 6: Remove Chimeras
```r
seqtab.nochim <- removeBimeraDenovo(seqtab)
```
**Purpose:** Detect and remove PCR chimeras (artificial sequences)

### Step 7: Assign Taxonomy
```r
taxa <- assignTaxonomy(seqtab.nochim, "silva_nr99_v138.1_train_set.fa.gz")
```
**Output:** Taxonomic assignments from Kingdom → Species

---

## 📊 Key Outputs

### 1. ASV Table (`ASV_table.rds`)
- Matrix: Samples × ASVs
- Each cell = read count for that ASV in that sample
- Used for downstream diversity analysis

### 2. Taxonomy Table (`taxonomy_SE.rds`)
- Taxonomic assignments for each ASV
- Levels: Kingdom, Phylum, Class, Order, Family, Genus, Species
- Based on SILVA v138.1 reference database

### 3. Quality Tracking Table
```r
         input  filtered  denoisedF  denoisedR  merged  nonchim
Sample1  45821     42156      41234      40987   38456    36892
Sample2  52134     48923      47821      47234   44123    42567
...
```
Shows read counts at each processing step

### 4. Quality Plots
- Per-base quality profiles (before/after filtering)
- Error rate learning curves
- ASV length distribution

---

## 🔧 Customization Guide

**For different 16S regions:**
1. Update `trimLeft` values (primer lengths)
2. Update `truncLen` based on quality profiles
3. Update `expected_v4_lengths` to your amplicon size range

**For different sequencing platforms:**
- Ion Torrent: Use `trimLeft = 0`, adjust `maxEE` upward
- PacBio/Nanopore: Use alternative long-read pipelines (not DADA2)

**For single-end data:**
Use `scripts/pipelines/single_end_16s.R` instead of main pipeline

---

## 🎓 What I Learned

**Technical Skills:**
- How 16S rRNA sequencing differs from whole genome sequencing
- Quality control for amplicon data (different from RNA-seq QC)
- ASV vs OTU approaches to microbiome analysis
- Importance of chimera removal in PCR-based methods

**Bioinformatics Concepts:**
- Error model learning in sequencing data
- Reference database limitations (genus-level vs species-level)
- Paired-end read merging requirements
- Batch effects in microbiome studies

**R Programming:**
- Portable project structure with `here` package
- Bioconductor workflow integration
- Quality tracking through multi-step pipelines
- Function-based data validation

---

## 🔮 Future Extensions

Planned additions to this project:

1. **Alpha/Beta Diversity Analysis**
   - Shannon diversity, Simpson index
   - PCoA ordination, PERMANOVA testing
   - Rarefaction curves

2. **Differential Abundance Testing**
   - DESeq2 for microbiome data
   - LEfSe analysis
   - Indicator species analysis

3. **Visualization**
   - Bar plots of taxonomic composition
   - Heatmaps of abundant taxa
   - Network analysis of co-occurrence

4. **Metadata Integration**
   - Sample grouping (treatment vs control)
   - Clinical variables correlation
   - Longitudinal time-series analysis

---

## 📝 License

This project is open source and available for educational purposes.  
Data sourced from NCBI SRA (public domain).  
SILVA database used under academic license.

---

## Path Convention

The analysis scripts now use `here::here()` so paths stay relative to the project root. Keep new raw inputs read-only, and write any intermediate or final results into `data/processed/` or `outputs/`.

## Notes

- The contents of the old `raw_data/` folder have been migrated into `data/raw/fastq/` and `data/processed/dada2/`; the legacy folder may still exist as an empty shell.
- SILVA references live in `data/external/reference/`.
- If you add a new script, place it under `scripts/` and keep it focused on one step of the workflow.
