#!/bin/bash
# SpaMFC Gene Correlation Analysis Example (CLI)
#
# This script demonstrates how to perform gene correlation analysis
# on spatial transcriptomics data.
#
# Usage:
#   bash run_corr.sh

set -e

# Configuration
INPUT_FILE="./data/input.h5ad"
OUTPUT_DIR="./results/correlation"

# Target genes (markers of interest)
TARGET_GENES="COL1A1,DCN,FAP,THY1,PDGFRB"

# DE genes (differentially expressed genes)
DE_GENES="POSTN,SPP1,IL6,CCL2,TNF"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "SpaMFC CLI Example: Gene Correlation Analysis"
echo "=========================================="
echo "Input: ${INPUT_FILE}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    echo "Please prepare your data file first."
    exit 1
fi

echo "Step 1: Running Spearman correlation analysis..."

SpaMFC corr \
    --input "${INPUT_FILE}" \
    --target-genes "${TARGET_GENES}" \
    --de-genes "${DE_GENES}" \
    --method spearman \
    --p-adjust fdr_bh \
    --threshold-p 0.05 \
    --min-corr 0.1 \
    --output "${OUTPUT_DIR}"

echo ""
echo "Step 2: Running Pearson correlation analysis..."

SpaMFC corr \
    --input "${INPUT_FILE}" \
    --target-genes "${TARGET_GENES}" \
    --de-genes "${DE_GENES}" \
    --method pearson \
    --p-adjust fdr_bh \
    --threshold-p 0.05 \
    --min-corr 0.1 \
    --output "${OUTPUT_DIR}/pearson"

echo ""
echo "=========================================="
echo "Gene correlation analysis completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
