#!/usr/bin/env bash
# Usage: bash scripts/fetch_metadata.sh  (run from project root)
# Fetches run-level metadata for PRJNA762524 via ffq and converts it to a validated
# CSV at data/processed/dada2/metadata.csv.
# Requires: ffq (pip install ffq), Rscript on PATH, jsonlite R package.

set -euo pipefail

if ! command -v ffq &>/dev/null; then
    echo "ffq not found; installing via pip..."
    pip install ffq
fi

mkdir -p outputs/qc

echo "Fetching metadata for PRJNA762524..."
ffq PRJNA762524 > outputs/qc/ffq_prjna762524.json
echo "Raw metadata JSON saved to outputs/qc/ffq_prjna762524.json"

echo "Processing metadata to CSV..."
Rscript scripts/process_metadata.R

echo "Done. Metadata saved to data/processed/dada2/metadata.csv"
