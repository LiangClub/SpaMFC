"""
SpaMFC Genomic Positions Module

Handles gene genomic position information for CNV inference.
Functions:
1. Load GTF files
2. Add chromosome information to AnnData
3. Validate genomic positions
4. Support multiple CNV methods (infercnv, copykat)
"""

import numpy as np
import pandas as pd
from typing import Optional, List, Dict, Any
import warnings
from pathlib import Path

try:
    import infercnvpy as cnv
    INFERCNVPY_AVAILABLE = True
except ImportError:
    INFERCNVPY_AVAILABLE = False
    warnings.warn("infercnvpy not installed")

try:
    import scanpy as sc
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    warnings.warn("scanpy not installed")


def add_gene_location(
    adata,
    gene_file: Optional[str] = None,
    method: str = "infercnv",
    species: str = "human",
    inplace: bool = True,
    verbose: bool = True
):
    """
    Dynamic gene location processing for CNV inference
    
    Supports multiple CNV methods:
    - copykat: CopyKAT automatically infers gene positions, no explicit file needed
    - infercnv: Requires gene annotation file
    
    Parameters:
        adata: AnnData object
        gene_file: Path to gene annotation file (TSV format: gene_symbol, chromosome, start, end)
        method: CNV method ('infercnv' or 'copykat')
        species: Species name ('human' or 'mouse')
        inplace: Whether to modify adata inplace
        verbose: Print progress information
    
    Returns:
        AnnData with genomic position information in .var
    """
    if not inplace:
        adata = adata.copy()
    
    if method == 'copykat':
        if verbose:
            print(f"[GeneLocation] CopyKAT automatically infers {species} gene positions")
        return adata
    
    if {'chromosome', 'start', 'end'}.issubset(adata.var.columns):
        if verbose:
            n_genes = adata.var['chromosome'].notna().sum()
            print(f"[GeneLocation] Detected existing chromosome position info ({n_genes} genes)")
        return adata
    
    if gene_file is None:
        raise ValueError("InferCNV requires gene annotation file. Please provide --gtf parameter.")
    
    gene_path = Path(gene_file)
    if not gene_path.exists():
        raise FileNotFoundError(f"Gene annotation file not found: {gene_file}")
    
    if verbose:
        print(f"[GeneLocation] Loading gene positions from: {gene_file}")
    
    try:
        gene_df = pd.read_csv(
            gene_file,
            sep='\t',
            header=None,
            names=['gene_symbol', 'chromosome', 'start', 'end'],
            dtype={'chromosome': str, 'start': int, 'end': int}
        )
    except Exception as e:
        raise ValueError(f"Failed to read gene annotation file: {str(e)}")
    
    if verbose:
        print(f"[GeneLocation] Loaded {len(gene_df)} gene entries")
    
    original_var = adata.var.copy()
    
    adata.var = adata.var.merge(
        gene_df.set_index('gene_symbol'),
        how='left',
        left_index=True,
        right_index=True
    )
    
    n_matched = adata.var['chromosome'].notna().sum()
    n_total = adata.n_vars
    if verbose:
        print(f"[GeneLocation] Matched {n_matched}/{n_total} genes ({n_matched/n_total*100:.1f}%)")
    
    valid_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
    
    adata = adata[:, adata.var['chromosome'].notna()]
    
    adata.var['chromosome'] = pd.Categorical(
        adata.var['chromosome'],
        categories=valid_chromosomes,
        ordered=True
    )
    
    adata.var = adata.var.sort_values(['chromosome', 'start'])
    adata = adata[:, adata.var.index]
    
    if verbose:
        n_final = adata.n_vars
        print(f"[GeneLocation] Final: {n_final} genes with valid positions")
        chromosomes = adata.var['chromosome'].unique()
        print(f"[GeneLocation] Chromosomes: {list(chromosomes)}")
    
    return adata


def load_gtf_file(gtf_file: str) -> pd.DataFrame:
    gtf_path = Path(gtf_file)
    
    if not gtf_path.exists():
        raise FileNotFoundError(f"GTF file not found: {gtf_file}")
    
    gtf_data = []
    skipped_lines = 0
    parse_errors = 0
    
    try:
        with open(gtf_file, 'r', encoding='utf-8') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#'):
                    continue
                
                try:
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        skipped_lines += 1
                        continue
                    
                    if parts[2] != 'gene':
                        continue
                    
                    chrom = parts[0]
                    try:
                        start = int(parts[3])
                        end = int(parts[4])
                    except ValueError:
                        parse_errors += 1
                        continue
                    
                    attributes = parts[8]
                    gene_name = None
                    
                    for attr in attributes.split(';'):
                        attr = attr.strip()
                        if attr.startswith('gene_name'):
                            try:
                                gene_name = attr.split('"')[1]
                            except IndexError:
                                continue
                            break
                        elif attr.startswith('gene_id'):
                            try:
                                gene_name = attr.split('"')[1]
                            except IndexError:
                                continue
                    
                    if gene_name:
                        gtf_data.append({
                            'gene_name': gene_name,
                            'chromosome': chrom,
                            'start': start,
                            'end': end
                        })
                except Exception as e:
                    parse_errors += 1
                    continue
    except (IOError, OSError, UnicodeDecodeError) as e:
        raise ValueError(f"Failed to read GTF file '{gtf_file}': {e}")
    
    if skipped_lines > 0:
        warnings.warn(f"Skipped {skipped_lines} lines with insufficient columns in GTF file")
    if parse_errors > 0:
        warnings.warn(f"Failed to parse {parse_errors} lines in GTF file")
    
    if not gtf_data:
        warnings.warn(f"No gene entries found in GTF file '{gtf_file}'")
    
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