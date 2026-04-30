#!/bin/bash
# SpaMFC Multi-Cell Type Analysis Example (CLI)
#
# This script demonstrates how to analyze multiple cell types
# with different feature configurations for each type.
#
# Usage:
#   bash run_multi.sh

set -e

# Configuration
INPUT_FILE="./data/input.h5ad"
CELLTYPE_COL="anno_cell2location_res"
OUTPUT_DIR="./results/multi"

# Cell types to analyze
CELLTYPES="Malignant cells,CAFs,ILC,ECs"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "SpaMFC CLI Example: Multi-Cell Type Analysis"
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
echo "Step 2: Running multi-cell type analysis..."

# Run multi-cell type analysis
SpaMFC run-multi \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}" \
    --celltypes "${CELLTYPES}" \
    --features spatial,expression \
    --fusion-method adaptive \
    --output "${OUTPUT_DIR}"

echo ""
echo "=========================================="
echo "Multi-cell type analysis completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
