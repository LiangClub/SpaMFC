"""
Functional Enrichment Analysis Module

Performs GO and KEGG pathway enrichment analysis for marker genes.
Uses gseapy's enrichr function.

Users can configure:
- gene_sets: Gene sets for enrichment ("GO_Biological_Process_2023", "KEGG_2021_Human")
- pval_threshold: P-value threshold for significant enrichment
- top_n: Number of top enrichment results to return
"""

import pandas as pd
from typing import Dict, List, Optional
import warnings

try:
    import gseapy as gp
except ImportError:
    warnings.warn("gseapy not installed, enrichment analysis will not work")


class EnrichmentAnalyzer:
    """Functional enrichment analyzer"""
    
    def __init__(
        self,
        gene_sets: Optional[List[str]] = None,
        pval_threshold: float = 0.05,
        top_n: int = 20
    ):
        """
        Initialize enrichment analyzer
        
        Parameters:
            gene_sets: List of gene sets for enrichment
            pval_threshold: P-value threshold
            top_n: Number of top results
        """
        self.gene_sets = gene_sets or [
            "GO_Biological_Process_2023",
            "KEGG_2021_Human"
        ]
        self.pval_threshold = pval_threshold
        self.top_n = top_n
    
    def analyze(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> Dict[str, Dict]:
        """
        Perform enrichment analysis for each subtype
        
        Parameters:
            markers_dict: Dictionary of marker genes for each subtype
        
        Returns:
            Dictionary of enrichment results for each subtype
        """
        try:
            import gseapy as gp
        except ImportError:
            warnings.warn("gseapy not installed, returning empty results")
            return {}
        
        enrichment_results = {}
        
        for subtype, genes in markers_dict.items():
            if len(genes) < 5:
                warnings.warn(f"Subtype {subtype} has fewer than 5 genes ({len(genes)}), skipping enrichment")
                enrichment_results[subtype] = {}
                continue
            
            subtype_results = {}
            
            for gene_set in self.gene_sets:
                try:
                    result = gp.enrichr(
                        gene_list=genes,
                        gene_sets=gene_set,
                        cutoff=self.pval_threshold,
                        no_plot=True
                    )
                    
                    if result.results is not None and len(result.results) > 0:
                        sig_results = result.results[
                            result.results["Adjusted P-value"] < self.pval_threshold
                        ]
                        if len(sig_results) > 0:
                            subtype_results[gene_set] = sig_results.head(self.top_n)
                        else:
                            subtype_results[gene_set] = pd.DataFrame()
                    else:
                        subtype_results[gene_set] = pd.DataFrame()
                    
                except Exception as e:
                    warnings.warn(f"Enrichment failed for {subtype} with {gene_set}: {e}")
                    subtype_results[gene_set] = pd.DataFrame()
            
            enrichment_results[subtype] = subtype_results
        
        return enrichment_results
    
    def analyze_single(
        self,
        genes: List[str],
        gene_set: Optional[str] = None
    ) -> pd.DataFrame:
        """
        Perform enrichment analysis for a single gene list
        
        Parameters:
            genes: List of genes
            gene_set: Gene set to use (default: first in self.gene_sets)
        
        Returns:
            DataFrame of enrichment results
        """
        try:
            import gseapy as gp
        except ImportError:
            warnings.warn("gseapy not installed")
            return pd.DataFrame()
        
        if len(genes) < 5:
            warnings.warn("Fewer than 5 genes provided")
            return pd.DataFrame()
        
        gene_set = gene_set or self.gene_sets[0]
        
        try:
            result = gp.enrichr(
                gene_list=genes,
                gene_sets=gene_set,
                cutoff=self.pval_threshold,
                no_plot=True
            )
            
            if result.results is not None:
                return result.results.head(self.top_n)
            else:
                return pd.DataFrame()
            
        except (ImportError, ValueError, RuntimeError, ConnectionError) as e:
            warnings.warn(f"Enrichment failed: {e}")
            return pd.DataFrame()
    
    def get_top_pathways(
        self,
        enrichment_results: Dict[str, Dict],
        gene_set: str = "KEGG_2021_Human",
        top_n: int = 5
    ) -> Dict[str, List[str]]:
        """
        Get top pathways for each subtype
        
        Parameters:
            enrichment_results: Dictionary of enrichment results
            gene_set: Gene set to extract pathways from
            top_n: Number of top pathways
        
        Returns:
            Dictionary of top pathways for each subtype
        """
        top_pathways = {}
        
        for subtype, results in enrichment_results.items():
            if gene_set in results and len(results[gene_set]) > 0:
                pathways = results[gene_set]["Term"].tolist()[:top_n]
                top_pathways[subtype] = pathways
            else:
                top_pathways[subtype] = []
        
        return top_pathways
    
    def get_enrichment_dataframe(
        self,
        enrichment_results: Dict[str, Dict],
        gene_set: str = "KEGG_2021_Human"
    ) -> pd.DataFrame:
        """
        Convert enrichment results to DataFrame
        
        Parameters:
            enrichment_results: Dictionary of enrichment results
            gene_set: Gene set to extract
        
        Returns:
            DataFrame of enrichment results
        """
        rows = []
        
        for subtype, results in enrichment_results.items():
            if gene_set in results and len(results[gene_set]) > 0:
                df = results[gene_set]
                for _, row in df.iterrows():
                    rows.append({
                        "subtype": subtype,
                        "term": row["Term"],
                        "pval": row["Adjusted P-value"],
                        "zscore": row.get("Z-score", 0),
                        "combined_score": row.get("Combined Score", 0),
                        "genes": row.get("Genes", "")
                    })
        
        return pd.DataFrame(rows)
    
    def compare_pathways(
        self,
        enrichment_results: Dict[str, Dict],
        gene_set: str = "KEGG_2021_Human"
    ) -> pd.DataFrame:
        """
        Compare pathways across subtypes
        
        Parameters:
            enrichment_results: Dictionary of enrichment results
            gene_set: Gene set to compare
        
        Returns:
            DataFrame of pathway comparison
        """
        pathway_sets = {}
        
        for subtype, results in enrichment_results.items():
            if gene_set in results and len(results[gene_set]) > 0:
                pathways = set(results[gene_set]["Term"].tolist())
                pathway_sets[subtype] = pathways
        
        all_pathways = set()
        for pathways in pathway_sets.values():
            all_pathways.update(pathways)
        
        comparison_rows = []
        for pathway in all_pathways:
            row = {"pathway": pathway}
            for subtype, pathways in pathway_sets.items():
                row[subtype] = 1 if pathway in pathways else 0
            comparison_rows.append(row)
        
        return pd.DataFrame(comparison_rows)