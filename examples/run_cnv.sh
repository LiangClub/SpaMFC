#!/bin/bash
# SpaMFC CNV Inference Analysis Example (CLI)
#
# This script demonstrates how to perform CNV inference
# for malignant cell analysis.
#
# Note: CNV analysis requires:
# 1. Gene annotation file (GTF or TSV format)
# 2. Reference cells (normal cells for comparison)
#
# Usage:
#   bash run_cnv.sh

set -e

# Configuration
INPUT_FILE="./data/input.h5ad"
GTF_FILE="./data/genes.gtf"
CELLTYPE_COL="anno_cell2location_res"
OUTPUT_DIR="./results/cnv"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "SpaMFC CLI Example: CNV Inference Analysis"
echo "=========================================="
echo "Input: ${INPUT_FILE}"
echo "GTF: ${GTF_FILE}"
echo "Output: ${OUTPUT_DIR}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    echo "Please prepare your data file first."
    exit 1
fi

# Check if GTF file exists
if [ ! -f "${GTF_FILE}" ]; then
    echo "Error: GTF file not found: ${GTF_FILE}"
    echo "Please provide a gene annotation file."
    echo "The GTF file should contain gene symbols and chromosome positions."
    echo "Example format:"
    echo "  gene_symbol  chromosome  start  end"
    echo "  COL1A1       chr17    48284435  48287193"
    echo "  DCN          chr12    91037540  91049780"
    exit 1
fi

echo "Step 1: Checking data..."
SpaMFC info \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}"

echo ""
echo "Step 2: Running CNV inference analysis..."

SpaMFC cnv \
    --input "${INPUT_FILE}" \
    --gtf "${GTF_FILE}" \
    --reference-key "${CELLTYPE_COL}" \
    --reference-cat "Normal cells,Fibroblasts" \
    --method infercnv \
    --output "${OUTPUT_DIR}" \
    --window-size 100 \
    --step-size 10 \
    --resolution 0.5 \
    --n-pcs 30 \
    --plot-heatmap \
    --plot-umap

echo ""
echo "=========================================="
echo "CNV inference analysis completed!"
echo "Results saved to: ${OUTPUT_DIR}"
echo "=========================================="
