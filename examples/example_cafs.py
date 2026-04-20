"""
Example: CAFs Subtype Analysis

This example demonstrates how to analyze CAFs (Cancer Associated Fibroblasts)
subtypes using spatial and expression features (no CNV).

Feature selection (user configurable):
- use_spatial: True (spatial neighborhood features)
- use_cnv: False (CAFs are not malignant, no CNV)
- use_expression: True (gene expression features)
"""

import scanpy as sc
import sys
sys.path.insert(0, "../src")

from SpaMFC import SpaMFCPipeline

def main():
    print("=" * 60)
    print("SpaMFC Example: CAFs Subtype Analysis")
    print("=" * 60)
    
    adata = sc.read_h5ad("./data/input.h5ad")
    
    print(f"Loaded data: {adata.shape[0]} cells, {adata.shape[1]} genes")
    
    pipeline = SpaMFCPipeline("./configs/cafs_config.yaml")
    
    pipeline.set_feature_usage(
        use_spatial=True,
        use_cnv=False,
        use_expression=True,
        use_niche=False
    )
    
    pipeline.print_config()
    
    adata = pipeline.run(adata, "CAFs")
    
    subtype_col = "CAFs_subtype"
    unified_col = "CAFs_subtype_unified"
    
    print("\nSubtype distribution:")
    print(adata.obs[subtype_col].value_counts())
    
    if unified_col in adata.obs:
        print("\nUnified subtype distribution:")
        print(adata.obs[unified_col].value_counts())
    
    results = pipeline.get_results("CAFs")
    
    print("\nMarker genes for each subtype:")
    for subtype, markers in results.get("markers", {}).items():
        print(f"  {subtype}: {markers[:5]}")
    
    adata.write_h5ad("./results/cafs/cafs_subtype_annotated.h5ad")
    
    print("\nResults saved to ./results/cafs/")
    print("=" * 60)

if __name__ == "__main__":
    main()