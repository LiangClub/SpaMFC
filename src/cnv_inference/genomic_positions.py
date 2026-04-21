"""
SpaMFC Genomic Positions Module

Handles gene genomic position information for CNV inference.
Functions:
1. Load GTF files
2. Add chromosome information to AnnData
3. Validate genomic positions
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict
import warnings
from pathlib import Path

try:
    import infercnvpy as cnv
    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False
    warnings.warn("infercnvpy not installed")


def load_gtf_file(gtf_file: str) -> pd.DataFrame:
    gtf_path = Path(gtf_file)
    
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    gtf_data = []
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            if parts[2] != 'gene':
                continue
            
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            
            attributes = parts[8]
            gene_name = None
            
            for attr in attributes.split(';'):
                attr = attr.strip()
                if attr.startswith('gene_name'):
                    gene_name = attr.split('"')[1]
                    break
                elif attr.startswith('gene_id'):
                    gene_name = attr.split('"')[1]
            
            if gene_name:
                gtf_data.append({
                    'gene_name': gene_name,
                    'chromosome': chrom,
                    'start': start,
                    'end': end
                })
    
    gtf_df = pd.DataFrame(gtf_data)
    
    return gtf_df


def add_genomic_positions(
    adata,
    gtf_file: str,
    inplace: bool = True
):
    if not INFERCNVPY_AVAILABLE:
        raise ImportError("infercnvpy is required. Install with: pip install infercnvpy")
    
    cnv.io.genomic_position_from_gtf(gtf_file, adata)
    
    return adata


def validate_genomic_positions(adata) -> Dict[str, Any]:
    validation_result = {
        'has_chromosome': False,
        'has_start': False,
        'has_end': False,
        'genes_with_positions': 0,
        'total_genes': adata.n_vars,
        'chromosomes': [],
        'missing_genes': []
    }
    
    if 'chromosome' in adata.var.columns:
        validation_result['has_chromosome'] = True
        validation_result['chromosomes'] = list(adata.var['chromosome'].dropna().unique())
    
    if 'start' in adata.var.columns:
        validation_result['has_start'] = True
    
    if 'end' in adata.var.columns:
        validation_result['has_end'] = True
    
    if validation_result['has_chromosome']:
        genes_with_pos = adata.var['chromosome'].notna().sum()
        validation_result['genes_with_positions'] = genes_with_pos
        
        missing_mask = adata.var['chromosome'].isna()
        if missing_mask.any():
            validation_result['missing_genes'] = list(adata.var_names[missing_mask][:20])
    
    validation_result['valid'] = (
        validation_result['has_chromosome'] and 
        validation_result['genes_with_positions'] > 0
    )
    
    return validation_result


def get_chromosome_order(adata) -> List[str]:
    if 'chromosome' not in adata.var.columns:
        return []
    
    chromosomes = adata.var['chromosome'].dropna().unique()
    
    chr_order = []
    for i in range(1, 23):
        chr_name = f'chr{i}'
        if chr_name in chromosomes:
            chr_order.append(chr_name)
    
    for chr_name in chromosomes:
        if chr_name not in chr_order and not chr_name.startswith('chr'):
            chr_order.append(chr_name)
    
    for chr_name in ['chrX', 'chrY', 'chrM']:
        if chr_name in chromosomes:
            chr_order.append(chr_name)
    
    return chr_order


def filter_genes_by_chromosome(
    adata,
    chromosomes: List[str],
    inplace: bool = False
):
    if 'chromosome' not in adata.var.columns:
        warnings.warn("Chromosome information not found")
        return adata
    
    mask = adata.var['chromosome'].isin(chromosomes)
    
    if inplace:
        adata._inplace_subset_var(mask)
        return adata
    else:
        return adata[:, mask].copy()