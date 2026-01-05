"""
Coordinate Mapper for VarPredBrowser

Maps coordinates between genomic positions and protein residues.
Uses precomputed mapping tables for fast lookups.
"""

import json
import math
import polars as pl
from pathlib import Path
from typing import Dict, List, Optional, Any


class CoordinateMapper:
    """
    Maps coordinates between genomic positions and protein residues.
    Uses precomputed mapping tables for fast lookups.
    """

    def __init__(self):
        self.protein_maps: Dict[str, pl.DataFrame] = {}
        self.structure_metadata: Dict[str, Dict] = {}

    def load_protein_map(self, gene_symbol: str, map_path: Path) -> bool:
        """Load a precomputed protein mapping table for a gene."""
        if not map_path.exists():
            return False

        self.protein_maps[gene_symbol.upper()] = pl.read_parquet(map_path)
        return True

    def load_structure_metadata(self, metadata_path: Path) -> bool:
        """Load structure metadata JSON."""
        if not metadata_path.exists():
            return False

        with open(metadata_path, 'r') as f:
            metadata = json.load(f)

        gene_symbol = metadata.get('gene_symbol', '').upper()
        if gene_symbol:
            self.structure_metadata[gene_symbol] = metadata

        return True

    def genomic_to_protein(self, gene: str, chrom: str, pos: int) -> Optional[Dict]:
        """
        Map a genomic position to protein residue.

        Returns:
            Dict with protein_residue, codon_position, ref_aa, or None if not found
        """
        gene_upper = gene.upper()
        if gene_upper not in self.protein_maps:
            return None

        protein_map = self.protein_maps[gene_upper]
        result = protein_map.filter(
            (pl.col('chrom') == chrom) &
            (pl.col('pos') == pos)
        )

        if len(result) == 0:
            return None

        row = result.to_dicts()[0]
        return {
            'protein_residue': row['protein_residue'],
            'codon_position': row['codon_position'],
            'ref_aa': row['ref_aa'],
            'cds_offset': row['cds_offset']
        }

    def protein_to_genomic(self, gene: str, residue: int) -> List[Dict]:
        """
        Map a protein residue to genomic positions (returns all 3 codon positions).

        Returns:
            List of dicts with chrom, pos, codon_position for each base in the codon
        """
        gene_upper = gene.upper()
        if gene_upper not in self.protein_maps:
            return []

        protein_map = self.protein_maps[gene_upper]
        results = protein_map.filter(pl.col('protein_residue') == residue)

        if len(results) == 0:
            return []

        return [
            {
                'chrom': row['chrom'],
                'pos': row['pos'],
                'codon_position': row['codon_position'],
                'ref_aa': row['ref_aa']
            }
            for row in results.to_dicts()
        ]

    def get_protein_range(self, gene: str, start_residue: int, end_residue: int) -> List[Dict]:
        """Get all genomic positions for a range of protein residues."""
        gene_upper = gene.upper()
        if gene_upper not in self.protein_maps:
            return []

        protein_map = self.protein_maps[gene_upper]
        results = protein_map.filter(
            (pl.col('protein_residue') >= start_residue) &
            (pl.col('protein_residue') <= end_residue)
        ).sort(['protein_residue', 'codon_position'])

        return results.to_dicts()

    def get_structure_metadata(self, gene: str) -> Optional[Dict]:
        """Get structure metadata for a gene."""
        return self.structure_metadata.get(gene.upper())

    def has_gene(self, gene: str) -> bool:
        """Check if a gene has protein mapping loaded."""
        return gene.upper() in self.protein_maps


def sanitize_float(value: Any) -> Optional[float]:
    """Convert NaN, Inf, -Inf to None for JSON serialization."""
    if value is None:
        return None
    try:
        if math.isnan(value) or math.isinf(value):
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def extract_constraint_pred(constraint_array: List, model_index: int) -> Optional[float]:
    """
    Extract prediction value from constraint_preds array.

    Args:
        constraint_array: List of constraint prediction tuples
        model_index: 0 for Constraint_1000, 1 for Core_1000, 2 for Complete_1000

    Returns:
        Average prediction value across all variants, or None if no data
    """
    if not constraint_array or len(constraint_array) == 0:
        return None

    # Map model_index to struct field name
    model_fields = ['_0', '_1', '_2']  # Constraint_1000, Core_1000, Complete_1000
    if model_index < 0 or model_index >= len(model_fields):
        return None

    field_name = model_fields[model_index]

    # Extract pred values from all variants and average them
    pred_values = []
    for variant in constraint_array:
        if variant and field_name in variant:
            model_data = variant[field_name]
            if model_data and '0' in model_data:  # '0' is the pred value
                pred_val = model_data['0']
                if pred_val is not None:
                    pred_values.append(float(pred_val))

    if len(pred_values) == 0:
        return None

    # Return average prediction across variants
    return sum(pred_values) / len(pred_values)


def extract_constraint_variants(variant_array: List) -> List[tuple]:
    """
    Extract variant prediction data from Constraint/Core/Complete columns.

    Args:
        variant_array: List of variant structs [{_0: allele, _1: pred, _2: n_pred}, ...]

    Returns:
        List of tuples: [(allele, pred), ...] sorted by pred value
    """
    if not variant_array or len(variant_array) == 0:
        return []

    variants = []
    for v in variant_array:
        if v and '_0' in v and '_1' in v:
            allele = v['_0']
            pred = float(v['_1']) if v['_1'] is not None else None
            if allele and pred is not None:
                variants.append((allele, pred))

    # Sort by pred value (lowest to highest for bottom-to-top stacking)
    variants.sort(key=lambda x: x[1])

    return variants
