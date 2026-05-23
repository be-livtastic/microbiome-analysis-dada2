#!/usr/bin/env bash
# Usage: bash scripts/seqkit_stats.sh  (run from project root)
# Runs seqkit stats on all raw FASTQ files and saves a tab-separated summary.
# Requires: seqkit on PATH (https://bioinf.shenwei.me/seqkit/)

set -euo pipefail

mkdir -p outputs/qc

echo "Running seqkit stats on data/raw/fastq/*.fastq.gz ..."
seqkit stats -a data/raw/fastq/*.fastq.gz | tee outputs/qc/seqkit_stats.tsv

echo ""
echo "Saved read statistics to outputs/qc/seqkit_stats.tsv"
