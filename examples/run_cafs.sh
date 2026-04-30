#!/bin/bash
# SpaMFC CAFs Analysis Example (CLI)
# 
# This script demonstrates how to analyze CAFs (Cancer Associated Fibroblasts)
# subtypes using the SpaMFC command line interface.
#
# Features used:
# - spatial: spatial neighborhood features
# - expression: gene expression features
# - CNV: disabled (CAFs are not malignant)
# - niche: disabled
#
# Usage:
#   bash run_cafs.sh

set -e

# Configuration
INPUT_FILE="./data/input.h5ad"
CELLTYPE_COL="anno_cell2location_res"
CELLTYPE="CAFs"
OUTPUT_DIR="./results/cafs"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "SpaMFC CLI Example: CAFs Analysis"
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

# Show data info
echo "Step 1: Checking data..."
SpaMFC info \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}"

echo ""
echo "Step 2: Running CAFs analysis..."

# Run CAFs analysis
SpaMFC run \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}" \
    --celltype "${CELLTYPE}" \
    --features spatial,expression \
    --fusion-method adaptive \
    --output "${OUTPUT_DIR}" \
    --n-clusters 5 \
    --nmf-components 5 \
    --nmf-runs 10 \
    --enable-marker-analysis \
    --enable-enrichment \
    --enable-unification \
    --save-plots

echo ""
echo "=========================================="
echo "CAFs analysis completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
