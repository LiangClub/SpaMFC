"""
Example: Multi-Cell Type Analysis

This example demonstrates how to analyze multiple cell types
with different feature configurations for each type.

Feature selection varies by cell type:
- Malignant cells: spatial + CNV + expression
- CAFs: spatial + expression (no CNV)
- Immune cells: expression only (no spatial, no CNV)
"""

import scanpy as sc

from SpaMFC import SpaMFCPipeline

def main():
    print("=" * 60)
    print("SpaMFC Example: Multi-Cell Type Analysis")
    print("=" * 60)
    
    adata = sc.read_h5ad("./data/input.h5ad")
    
    print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    pipeline = SpaMFCPipeline("./configs/default_config.yaml")
    
    celltypes = [
        "Malignant cells",
        "CAFs",
        "ILC",
        "ECs",
        "Epithelial Cell"
    ]
    
    tumor_markers = {
        "Malignant cells": ["EPCAM", "KRT8", "KRT18", "KRT19"]
    }
    
    adata = pipeline.run_multi_celltypes(adata, celltypes, tumor_markers)
    
    print("\n" + "=" * 60)
    print("Subtype Summary for All Cell Types")
    print("=" * 60)
    
    for celltype in celltypes:
        subtype_col = f"{celltype}_subtype"
        
        if subtype_col in adata.obs:
            print(f"\n{celltype}:")
            print(adata.obs[subtype_col].value_counts())
    
    adata.write_h5ad("./results/all_celltypes_subtype_annotated.h5ad")
    
    pipeline.save_results("./results/summary/")
    
    print("\nAll results saved to ./results/")
    print("=" * 60)

if __name__ == "__main__":
    main()