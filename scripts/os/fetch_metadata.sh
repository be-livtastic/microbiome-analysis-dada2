#!/usr/bin/env bash
# Fetches run-level metadata for PRJNA762524 via ffq using SRR accessions
set -euo pipefail

cd "$(dirname "$0")/../.." # run from project root regardless of where called from

if ! command -v ffq &>/dev/null; then
    echo "ffq not found; installing..."
    pip install ffq
fi

mkdir -p outputs/qc

echo "Extracting SRR accessions from download script..."
# Extract only complete SRR accessions (8 digits minimum)
grep -oP 'SRR\d{8,}' scripts/os/ena-file-download-read_run-PRJNA762524-fastq_ftp-20260407-1506.sh \
    | sort -u > outputs/qc/srr_accessions.txt

wc -l outputs/qc/srr_accessions.txt   # should now be 72
cat outputs/qc/srr_accessions.txt     # verify all look like SRR15862xxx

echo "Found $(wc -l < outputs/qc/srr_accessions.txt) unique SRR accessions"

echo "Fetching metadata via ffq..."
ffq $(cat outputs/qc/srr_accessions.txt | tr '\n' ' ') > outputs/qc/ffq_prjna762524.json

echo "Raw metadata JSON saved to outputs/qc/ffq_prjna762524.json"

echo "Processing metadata to CSV..."
Rscript scripts/process_metadata.R

echo "Done. Metadata saved to data/processed/dada2/metadata.csv"