# Master Guide: 16S rRNA Microbiome Analysis
## Practising with the Okwelogu et al. (2021) IVF Microbiome Dataset

**Purpose:** Build the biological knowledge and bioinformatics skills to confidently discuss, analyse, and interpret 16S rRNA microbiome data 

**Dataset:** Okwelogu et al. (2021). "Microbiome Compositions From Infertile Couples Seeking In Vitro Fertilization, Using 16S rRNA Gene Sequencing Methods: Any Correlation to Clinical Outcomes?" *Frontiers in Cellular and Infection Microbiology*, 11:709372.

- **PubMed:** https://pubmed.ncbi.nlm.nih.gov/34660337/
- **Full text (free):** https://doi.org/10.3389/fcimb.2021.709372
- **Raw data:** NCBI SRA, Project Accession PRJNA762524
- **Direct SRA link:** https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA762524

---

> **Author:** Olivia Williams  
> **Technology:** R, DADA2, SILVA database  
> **Domain:** Microbiome bioinformatics, 16S rRNA gene sequencing  
> **Project Type:** Reproducible bioinformatics workflow for amplicon sequencing data

---

## 📋 Project Overview

This project implements a complete **16S rRNA gene amplicon sequencing analysis workflow** using the industry-standard DADA2 pipeline in R. Starting from raw FASTQ files, the pipeline performs quality control, denoising, taxonomic classification, and generates publication-ready microbiome composition data.

**What this pipeline does:**
- Processes raw Illumina paired-end sequencing data from microbiome samples
- Performs quality filtering and adapter trimming
- Infers Amplicon Sequence Variants (ASVs) using DADA2's error-correcting algorithm
- Assigns taxonomy using the SILVA reference database (v138.1)
- Tracks read counts through each processing step for quality assurance

---

## 🎯 Learning Objectives & Skills Demonstrated

**Bioinformatics Skills:**
- Raw sequencing data (FASTQ) processing and quality assessment
- Understanding 16S rRNA gene amplicon sequencing workflow
- Reference database usage (SILVA) for taxonomic classification
- ASV inference vs traditional OTU clustering approaches

**R Programming:**
- Bioconductor package management (dada2, ShortRead, Biostrings)
- Reproducible project structure with `here` package
- Error handling and data validation
- Function-based workflow design

**Best Practices:**
- Version control ready (Git/GitHub)
- Portable file paths (works on any machine)
- Documented decision points (trim lengths, quality thresholds)
- Comprehensive quality tracking throughout pipeline

---

## 📂 Project Structure

```
microbiome_practice/
├── README.md                    # This file
├── scripts/
│   ├── pipelines/
│   │   ├── load_packages.R      # One-time setup: install required packages
│   │   ├── dada2_pipeline.R     # Main paired-end analysis workflow
│   │   ├── single_end_16s.R     # Alternative: single-end only workflow
│   │   └── troubleshooting.R    # Diagnostic checks for data quality issues
│   └── download/
│       └── ena-file-download-*.sh  # Automated SRA data download script
├── data/
│   ├── raw/fastq/               # Downloaded FASTQ files (not tracked in Git)
│   ├── processed/dada2/         # Filtered reads, ASV table, taxonomy
│   └── external/reference/      # SILVA database files (download separately)
├── outputs/
│   ├── figures/                 # Quality plots, error rate plots
│   └── tables/                  # ASV abundance tables, taxonomy assignments
└── docs/
    └── project_structure.md     # Detailed methodology documentation
```

---

## 🚀 Quick Start

### Prerequisites

**Required software:**
- R (≥ 4.0)
- RStudio (recommended)

**Required R packages:**
- Bioconductor: `dada2`, `ShortRead`, `Biostrings`, `phyloseq`
- CRAN: `tidyverse`, `here`, `ggplot2`

### Step 1: Install Dependencies

```r
# Run this once to install all required packages
source("scripts/pipelines/load_packages.R")
```

This script automatically installs missing packages from both CRAN and Bioconductor.

### Step 2: Download Reference Database

Download SILVA v138.1 training sets and place in `data/external/reference/`:

```bash
# Required files:
silva_nr99_v138.1_train_set.fa.gz
silva_nr99_v138.1_wSpecies_train_set.fa.gz  
silva_species_assignment_v138.1.fa.gz
```

**Download from:** https://zenodo.org/record/4587955

### Step 3: Obtain Raw Data

**Option A: Use my test dataset (PRJNA762524)**
```bash
# Download from NCBI SRA using provided script
bash scripts/download/ena-file-download-read_run-PRJNA762524-fastq_ftp-*.sh
```

**Option B: Use your own 16S data**
- Place paired-end FASTQ files in `data/raw/fastq/`
- Files must follow naming convention: `SAMPLE_1.fastq`, `SAMPLE_2.fastq`

### Step 4: Run the Pipeline

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

## 📖 Dataset Information

**Project:** PRJNA762524  
**Organism:** Human gut microbiome  
**Sequencing:** Illumina paired-end (2×250 bp)  
**Target:** 16S rRNA V4 region  
**Samples:** 4 samples (can be expanded)

**SRA Accessions:**
- SRR15862403
- SRR15862405  
- SRR15862407
- SRR15862409

---

## 🐛 Troubleshooting

### Common Issues

**1. "Reference FASTA not found"**
- Download SILVA files from Zenodo
- Place in `data/external/reference/` folder
- Verify file names match exactly

**2. "No paired FASTQ files found"**
- Check files are in `data/raw/fastq/`
- Verify naming: `*_1.fastq` and `*_2.fastq`
- Ensure files are uncompressed (or use `.fastq.gz`)

**3. "Error in learnErrors: Not enough reads"**
- Increase number of samples (need ≥3 for error learning)
- Or reduce `truncLen` (retaining more reads after filtering)

**4. Low merge rate (<50%)**
- Increase `truncLen` for overlap
- Check primer trimming removed all adapters
- Use `troubleshooting.R` to inspect read sequences

**For detailed diagnostics:**
```r
source("scripts/pipelines/troubleshooting.R")
```

---

## 📚 References & Resources

**DADA2 Documentation:**
- Official tutorial: https://benjjneb.github.io/dada2/tutorial.html
- Paper: Callahan et al. (2016) Nature Methods

**SILVA Database:**
- Website: https://www.arb-silva.de/
- Training sets: https://zenodo.org/record/4587955

**Microbiome Analysis in R:**
- phyloseq package: https://joey711.github.io/phyloseq/
- Statistical analysis: vegan, microbiome, DESeq2

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

## 🤝 Acknowledgments

- **DADA2 developers** (Benjamin Callahan et al.) for the excellent pipeline
- **SILVA team** for maintaining the reference database
- **Bioconductor community** for microbiome analysis tools

---

## Path Convention

The analysis scripts now use `here::here()` so paths stay relative to the project root. Keep new raw inputs read-only, and write any intermediate or final results into `data/processed/` or `outputs/`.

## Notes

- The contents of the old `raw_data/` folder have been migrated into `data/raw/fastq/` and `data/processed/dada2/`; the legacy folder may still exist as an empty shell.
- SILVA references live in `data/external/reference/`.
- If you add a new script, place it under `scripts/` and keep it focused on one step of the workflow.
