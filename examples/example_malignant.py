"""
Example: Malignant Cell Subtype Analysis

This example demonstrates how to analyze malignant cell subtypes
using all available features (spatial + CNV + expression).

Feature selection (user configurable):
- use_spatial: True (spatial neighborhood features)
- use_cnv: True (CNV features for malignant cells)
- use_expression: True (gene expression features)
"""

import scanpy as sc
import sys
sys.path.insert(0, "../src")

from SpaMFC import SpaMFCPipeline

def main():
    print("=" * 60)
    print("SpaMFC Example: Malignant Cell Subtype Analysis")
    print("=" * 60)
    
    adata = sc.read_h5ad("./data/input.h5ad")
    
    print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    pipeline = SpaMFCPipeline("./configs/malignant_config.yaml")
    
    pipeline.set_feature_usage(
        use_spatial=True,
        use_cnv=True,
        use_expression=True,
        use_niche=False
    )
    
    pipeline.print_config()
    
    tumor_markers = [
        "EPCAM", "KRT8", "KRT18", "KRT19",
        "MUC1", "CEACAM5", "VIL1"
    ]
    
    adata = pipeline.run(adata, "Malignant cells", markers=tumor_markers)
    
    subtype_col = "Malignant cells_subtype"
    unified_col = "Malignant cells_subtype_unified"
    
    print("\nSubtype distribution:")
    print(adata.obs[subtype_col].value_counts())
    
    if unified_col in adata.obs:
        print("\nUnified subtype distribution:")
        print(adata.obs[unified_col].value_counts())
    
    results = pipeline.get_results("Malignant cells")
    
    print("\nMarker genes for each subtype:")
    for subtype, markers in results.get("markers", {}).items():
        print(f"  {subtype}: {markers[:5]}")
    
    adata.write_h5ad("./results/malignant/malignant_subtype_annotated.h5ad")
    
    print("\nResults saved to ./results/malignant/")
    print("=" * 60)

if __name__ == "__main__":
    main()