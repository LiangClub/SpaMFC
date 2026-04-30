#!/bin/bash
# SpaMFC Malignant Cell Analysis Example (CLI)
#
# This script demonstrates how to analyze malignant cell subtypes
# using all available features (spatial + CNV + expression).
#
# Features used:
# - spatial: spatial neighborhood features
# - expression: gene expression features
# - CNV: enabled (CNV features for malignant cells)
# - niche: disabled
#
# Usage:
#   bash run_malignant.sh

set -e

# Configuration
INPUT_FILE="./data/input.h5ad"
CELLTYPE_COL="anno_cell2location_res"
CELLTYPE="Malignant cells"
OUTPUT_DIR="./results/malignant"

# Tumor marker genes for CNV filtering
MARKERS="EPCAM KRT8 KRT18 KRT19 MUC1 CEACAM5"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "SpaMFC CLI Example: Malignant Cell Analysis"
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

# Check if CNV key exists in the data
echo "Step 1: Checking data..."
SpaMFC info \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}"

echo ""
echo "Step 2: Running malignant cell analysis with CNV..."

# Run malignant cell analysis with all features
SpaMFC run \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}" \
    --celltype "${CELLTYPE}" \
    --features spatial,cnv,expression \
    --fusion-method adaptive \
    --output "${OUTPUT_DIR}" \
    --n-clusters 5 \
    --nmf-components 5 \
    --nmf-runs 10 \
    --markers ${MARKERS} \
    --enable-marker-analysis \
    --enable-enrichment \
    --enable-unification \
    --save-plots

echo ""
echo "=========================================="
echo "Malignant cell analysis completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
