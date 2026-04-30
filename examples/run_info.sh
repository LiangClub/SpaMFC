#!/bin/bash
# SpaMFC Data Info Example (CLI)
#
# This script demonstrates how to check data information
# before running analysis.
#
# Usage:
#   bash run_info.sh [INPUT_FILE]

set -e

# Configuration
INPUT_FILE="${1:-./data/input.h5ad}"
CELLTYPE_COL="${2:-anno_cell2location_res}"

echo "=========================================="
echo "SpaMFC CLI Example: Data Information"
echo "=========================================="
echo "Input: ${INPUT_FILE}"
echo "Cell Type Column: ${CELLTYPE_COL}"
echo ""

# Check if input file exists
if [ ! -f "${INPUT_FILE}" ]; then
    echo "Error: Input file not found: ${INPUT_FILE}"
    echo "Usage: bash run_info.sh [INPUT_FILE] [CELLTYPE_COL]"
    exit 1
fi

# Show data info
SpaMFC info \
    --input "${INPUT_FILE}" \
    --celltype-col "${CELLTYPE_COL}"

echo ""
echo "=========================================="
echo "Note: Use the cell types shown above for analysis."
echo "=========================================="
