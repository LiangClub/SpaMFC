"""
Report Generation Module

Generates comprehensive reports for subtype analysis including:
- Subtype statistics
- Marker gene lists
- Functional enrichment results
- Niche analysis results
- Cross-sample unification results

Users can configure:
- output_format: Report format (xlsx, html, pdf)
- include_plots: Whether to include plots in report
"""

import os
import html
import pandas as pd
from pathlib import Path
from typing import Optional, Dict, List
from datetime import datetime
import warnings


class ReportGenerator:
    """Report generator for subtype analysis"""
    
    def __init__(
        self,
        output_format: str = "xlsx",
        include_plots: bool = True
    ):
        """
        Initialize report generator
        
        Parameters:
            output_format: Output format (xlsx, html)
            include_plots: Whether to include plots
        """
        self.output_format = output_format
        self.include_plots = include_plots
    
    def generate(
        self,
        adata,
        celltype: str,
        markers_dict: Dict[str, List[str]],
        enrichment_results: Dict[str, Dict],
        output_dir: str,
        mapping_table: Optional[pd.DataFrame] = None,
        niche_profiles: Optional[Dict] = None
    ):
        """
        Generate comprehensive report
        
        Parameters:
            adata: AnnData object
            celltype: Cell type name
            markers_dict: Dictionary of marker genes
            enrichment_results: Dictionary of enrichment results
            output_dir: Output directory
            mapping_table: Subtype mapping table
            niche_profiles: Niche profiles dictionary
        """
        output_path = Path(output_dir) / celltype.replace(" ", "_") / "reports"
        output_path.mkdir(parents=True, exist_ok=True)
        
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        if self.output_format == "xlsx":
            self._generate_xlsx_report(
                adata, celltype, markers_dict, enrichment_results,
                output_path, mapping_table, niche_profiles, timestamp
            )
        elif self.output_format == "html":
            self._generate_html_report(
                adata, celltype, markers_dict, enrichment_results,
                output_path, mapping_table, niche_profiles, timestamp
            )
        else:
            warnings.warn(f"Unknown format: {self.output_format}")
    
    def _generate_xlsx_report(
        self,
        adata,
        celltype: str,
        markers_dict: Dict[str, List[str]],
        enrichment_results: Dict[str, Dict],
        output_path: Path,
        mapping_table: Optional[pd.DataFrame],
        niche_profiles: Optional[Dict],
        timestamp: str
    ):
        """Generate Excel report"""
        report_path = output_path / f"{celltype}_report_{timestamp}.xlsx"
        
        with pd.ExcelWriter(report_path, engine="openpyxl") as writer:
            subtype_col = f"{celltype}_subtype"
            
            if subtype_col in adata.obs:
                subtype_counts = adata.obs[subtype_col].value_counts()
                counts_df = pd.DataFrame({
                    "subtype": subtype_counts.index,
                    "count": subtype_counts.values
                })
                counts_df.to_excel(writer, sheet_name="Subtype Counts", index=False)
            
            markers_df = self._markers_to_dataframe(markers_dict)
            if not markers_df.empty:
                markers_df.to_excel(writer, sheet_name="Marker Genes", index=False)
            
            enrichment_df = self._enrichment_to_dataframe(enrichment_results)
            if not enrichment_df.empty:
                enrichment_df.to_excel(writer, sheet_name="Enrichment", index=False)
            
            if mapping_table is not None and not mapping_table.empty:
                mapping_table.to_excel(writer, sheet_name="Subtype Mapping", index=False)
            
            if niche_profiles is not None:
                niche_df = self._niche_to_dataframe(niche_profiles)
                if not niche_df.empty:
                    niche_df.to_excel(writer, sheet_name="Niche Profiles", index=False)
            
            summary_df = self._generate_summary(
                adata, celltype, markers_dict, enrichment_results
            )
            summary_df.to_excel(writer, sheet_name="Summary", index=False)
    
    def _generate_html_report(
        self,
        adata,
        celltype: str,
        markers_dict: Dict[str, List[str]],
        enrichment_results: Dict[str, Dict],
        output_path: Path,
        mapping_table: Optional[pd.DataFrame],
        niche_profiles: Optional[Dict],
        timestamp: str
    ):
        """Generate HTML report"""
        report_path = output_path / f"{celltype}_report_{timestamp}.html"
        
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>SpaMFC Report - {celltype}</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 20px; }}
        h1 {{ color: #2c3e50; }}
        h2 {{ color: #34495e; border-bottom: 1px solid #bdc3c7; }}
        table {{ border-collapse: collapse; width: 100%; margin: 10px 0; }}
        th, td {{ border: 1px solid #bdc3c7; padding: 8px; text-align: left; }}
        th {{ background-color: #3498db; color: white; }}
        tr:nth-child(even) {{ background-color: #ecf0f1; }}
        .summary {{ background-color: #f9f9f9; padding: 15px; border-radius: 5px; }}
    </style>
</head>
<body>
    <h1>SpaMFC Subtype Analysis Report</h1>
    <p>Cell Type: <strong>{html.escape(str(celltype))}</strong></p>
    <p>Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    
    <div class="summary">
        <h2>Summary</h2>
        <p>Total subtypes identified: {len(markers_dict)}</p>
        <p>Total marker genes analyzed: {sum(len(g) for g in markers_dict.values())}</p>
    </div>
    
    <h2>Marker Genes</h2>
    {self._markers_to_html(markers_dict)}
    
    <h2>Functional Enrichment</h2>
    {self._enrichment_to_html(enrichment_results)}
    
    {self._mapping_to_html(mapping_table) if mapping_table else ""}
    
    {self._niche_to_html(niche_profiles) if niche_profiles else ""}
    
</body>
</html>
"""
        
        try:
            with open(report_path, "w", encoding="utf-8") as f:
                f.write(html_content)
        except (IOError, OSError, PermissionError) as e:
            warnings.warn(f"Failed to write HTML report: {e}")
    
    def _markers_to_dataframe(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> pd.DataFrame:
        """Convert markers dictionary to DataFrame"""
        rows = []
        
        for subtype, genes in markers_dict.items():
            for i, gene in enumerate(genes):
                rows.append({
                    "subtype": subtype,
                    "rank": i + 1,
                    "gene": gene
                })
        
        return pd.DataFrame(rows)
    
    def _enrichment_to_dataframe(
        self,
        enrichment_results: Dict[str, Dict]
    ) -> pd.DataFrame:
        """Convert enrichment results to DataFrame"""
        rows = []
        
        for subtype, results in enrichment_results.items():
            for gene_set, df in results.items():
                if df is not None and len(df) > 0:
                    for _, row in df.iterrows():
                        rows.append({
                            "subtype": subtype,
                            "gene_set": gene_set,
                            "term": row.get("Term", ""),
                            "pval": row.get("Adjusted P-value", ""),
                            "zscore": row.get("Z-score", ""),
                            "genes": row.get("Genes", "")
                        })
        
        return pd.DataFrame(rows)
    
    def _niche_to_dataframe(
        self,
        niche_profiles: Dict[str, Dict]
    ) -> pd.DataFrame:
        """Convert niche profiles to DataFrame"""
        rows = []
        
        for subtype, profile in niche_profiles.items():
            composition = profile.get("composition", {})
            dominant = profile.get("dominant_niches", [])
            primary = profile.get("primary_niche", "")
            
            rows.append({
                "subtype": subtype,
                "primary_niche": primary,
                "dominant_niches": ",".join(dominant),
                "niche_diversity": len(dominant)
            })
        
        return pd.DataFrame(rows)
    
    def _generate_summary(
        self,
        adata,
        celltype: str,
        markers_dict: Dict[str, List[str]],
        enrichment_results: Dict[str, Dict]
    ) -> pd.DataFrame:
        """Generate summary DataFrame"""
        subtype_col = f"{celltype}_subtype"
        
        total_cells = 0
        n_subtypes = len(markers_dict)
        
        if subtype_col in adata.obs:
            subtype_mask = adata.obs[subtype_col].notna()
            total_cells = subtype_mask.sum()
        
        total_markers = sum(len(g) for g in markers_dict.values())
        
        n_enriched = sum(
            1 for r in enrichment_results.values()
            if r is not None and any(len(df) > 0 for df in r.values())
        )
        
        summary_data = {
            "item": [
                "Cell Type",
                "Total Cells Analyzed",
                "Number of Subtypes",
                "Total Marker Genes",
                "Subtypes with Enrichment"
            ],
            "value": [
                celltype,
                total_cells,
                n_subtypes,
                total_markers,
                n_enriched
            ]
        }
        
        return pd.DataFrame(summary_data)
    
    def _markers_to_html(
        self,
        markers_dict: Dict[str, List[str]]
    ) -> str:
        """Convert markers to HTML table"""
        if not markers_dict:
            return "<p>No marker genes identified</p>"
        
        html_table = "<table><tr><th>Subtype</th><th>Top Markers</th></tr>"
        
        for subtype, genes in markers_dict.items():
            top_genes = genes[:10]
            escaped_subtype = html.escape(str(subtype))
            escaped_genes = html.escape(', '.join(str(g) for g in top_genes))
            html_table += f"<tr><td>{escaped_subtype}</td><td>{escaped_genes}</td></tr>"
        
        html_table += "</table>"
        return html_table
    
    def _enrichment_to_html(
        self,
        enrichment_results: Dict[str, Dict]
    ) -> str:
        """Convert enrichment to HTML"""
        if not enrichment_results:
            return "<p>No enrichment results</p>"
        
        html_table = ""
        
        for subtype, results in enrichment_results.items():
            html_table += f"<h3>{html.escape(str(subtype))}</h3>"
            
            for gene_set, df in results.items():
                if df is not None and len(df) > 0:
                    html_table += f"<h4>{html.escape(str(gene_set))}</h4>"
                    html_table += "<table><tr><th>Term</th><th>P-value</th></tr>"
                    
                    for _, row in df.head(5).iterrows():
                        term = html.escape(str(row.get('Term', '')))
                        pval = html.escape(str(row.get('Adjusted P-value', '')))
                        html_table += f"<tr><td>{term}</td><td>{pval}</td></tr>"
                    
                    html_table += "</table>"
        
        return html_table
    
    def _mapping_to_html(
        self,
        mapping_table: pd.DataFrame
    ) -> str:
        """Convert mapping to HTML"""
        if mapping_table is None or mapping_table.empty:
            return ""
        
        html_table = "<h2>Subtype Mapping</h2>"
        html_table += "<table><tr><th>Original</th><th>Unified</th><th>Sample</th></tr>"
        
        for _, row in mapping_table.iterrows():
            orig = html.escape(str(row.get('original_subtype', '')))
            unified = html.escape(str(row.get('unified_subtype', '')))
            sample = html.escape(str(row.get('sample', '')))
            html_table += f"<tr><td>{orig}</td><td>{unified}</td><td>{sample}</td></tr>"
        
        html_table += "</table>"
        return html_table
    
    def _niche_to_html(
        self,
        niche_profiles: Dict[str, Dict]
    ) -> str:
        """Convert niche profiles to HTML"""
        if niche_profiles is None:
            return ""
        
        html_table = "<h2>Niche Profiles</h2>"
        html_table += "<table><tr><th>Subtype</th><th>Primary Niche</th><th>Dominant Niches</th></tr>"
        
        for subtype, profile in niche_profiles.items():
            primary = html.escape(str(profile.get("primary_niche", "")))
            dominant = html.escape(",".join(str(n) for n in profile.get("dominant_niches", [])))
            escaped_subtype = html.escape(str(subtype))
            html_table += f"<tr><td>{escaped_subtype}</td><td>{primary}</td><td>{dominant}</td></tr>"
        
        html_table += "</table>"
        return html_table