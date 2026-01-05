#!/usr/bin/env python3
"""
FastAPI backend for VarPredBrowser

Provides coordinate mapping and data retrieval endpoints for the compressed
genome viewer with protein structure visualization support.
"""

import math
import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import polars as pl
import pyBigWig
from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse

from .config import get_config, get_data_dir, get_bigwig_dir, get_project_root
from .models import (
    FilterInfo, Position, WindowResponse, GeneSearchResult, PositionSearchResult
)
from .coordinate_mapper import (
    CoordinateMapper, sanitize_float, extract_constraint_variants
)
from .track_tree import (
    TRACK_TREE, FILTERS, CONSTRAINT_STACKED_FIELDS,
    simplify_track_name, categorize_track
)


# Global data cache
axis_tables: Dict[str, pl.DataFrame] = {}
gene_indexes: Dict[str, pl.DataFrame] = {}
coord_mapper = CoordinateMapper()

# Legacy - kept for compatibility
CHROMOSOME = 'all'


app = FastAPI(
    title="VarPredBrowser API",
    description="API for exploring genomic data in compressed coordinate space with protein structure visualization",
    version="0.1.0"
)

# Enable CORS for frontend development
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # In production, specify your frontend domain
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


# =============================================================================
# Startup: Load data into memory
# =============================================================================

@app.on_event("startup")
async def load_data():
    """Load axis tables and gene indexes for all filters into memory on startup."""
    global axis_tables, gene_indexes, coord_mapper

    config = get_config()
    data_dir = get_data_dir()

    print(f"Loading data from {data_dir}...")

    for filter_id in FILTERS.keys():
        # Try all-chromosomes file first, then fall back to chr1-specific
        axis_file = data_dir / f'{filter_id}.parquet'
        if not axis_file.exists():
            axis_file = data_dir / f'{filter_id}_chr1.parquet'

        if not axis_file.exists():
            print(f"Warning: Axis table not found for {filter_id}")
            continue

        print(f"Loading {axis_file.name}...")
        axis_tables[filter_id] = pl.read_parquet(axis_file)
        print(f"  Loaded {filter_id}: {len(axis_tables[filter_id]):,} positions, {len(axis_tables[filter_id].columns)} columns")

        # Load gene index
        gene_file = data_dir / f'gene_index_{filter_id}.parquet'
        if not gene_file.exists():
            gene_file = data_dir / f'gene_index_{filter_id}_chr1.parquet'

        if gene_file.exists():
            gene_indexes[filter_id] = pl.read_parquet(gene_file)
            print(f"  Gene index: {len(gene_indexes[filter_id]):,} genes")

    # Load protein coordinate maps
    print("\nLoading protein coordinate maps...")
    protein_map_files = list(data_dir.glob('*_protein_map.parquet'))
    for map_file in protein_map_files:
        gene_name = map_file.stem.replace('_protein_map', '').upper()
        if coord_mapper.load_protein_map(gene_name, map_file):
            print(f"  Loaded protein map for {gene_name}")

    # Load structure metadata
    structure_metadata_file = data_dir / 'structure_metadata.json'
    if structure_metadata_file.exists():
        if coord_mapper.load_structure_metadata(structure_metadata_file):
            print(f"  Loaded structure metadata")

    print(f"\nReady! Loaded {len(axis_tables)} filter(s), {len(coord_mapper.protein_maps)} protein map(s)")


# =============================================================================
# Health & Status Endpoints
# =============================================================================

@app.get("/api/health")
async def health():
    """Health check endpoint."""
    return {
        "status": "ok",
        "message": "VarPredBrowser API",
        "data_loaded": len(axis_tables) > 0,
        "filters_loaded": list(axis_tables.keys())
    }


@app.get("/api/filters", response_model=List[FilterInfo])
async def get_filters():
    """Get available filters."""
    if not axis_tables:
        raise HTTPException(status_code=503, detail="Data not loaded")

    return [
        FilterInfo(
            id=filter_id,
            name=f"{FILTERS[filter_id]['name']} (all chromosomes)",
            description=FILTERS[filter_id]['description'],
            total_positions=len(axis_tables[filter_id]),
            chromosome="all"
        )
        for filter_id in axis_tables.keys()
    ]


@app.get("/api/track-tree")
async def get_track_tree():
    """Get hierarchical track tree for frontend."""
    return TRACK_TREE


@app.get("/api/tracks")
async def get_tracks():
    """Get all available tracks from BigWig files and axis table columns."""
    tracks = []
    bigwig_dir = get_bigwig_dir()

    # Get tracks from BigWig directory
    if bigwig_dir.exists():
        for bw_file in sorted(bigwig_dir.glob("*.bw")):
            track_id = bw_file.stem
            tracks.append({
                "id": track_id,
                "name": simplify_track_name(track_id),
                "source": "bigwig",
                "category": categorize_track(track_id)
            })

    # Also add columns from axis tables
    if axis_tables:
        all_columns = set()
        for axis_table in axis_tables.values():
            all_columns.update(axis_table.columns)

        for col in all_columns:
            if col not in ['filtered_idx', 'chrom', 'pos', 'gene_symbol']:
                if not any(t['id'] == col for t in tracks):
                    tracks.append({
                        "id": col,
                        "name": simplify_track_name(col),
                        "source": "axis_table",
                        "category": categorize_track(col)
                    })

    return {"tracks": tracks, "total": len(tracks)}


# =============================================================================
# Genomic Data Endpoints
# =============================================================================

@app.get("/api/filtered-window", response_model=WindowResponse)
async def get_filtered_window(
    filter_id: str = Query(..., description="Filter ID"),
    start: int = Query(..., ge=0, description="Start position in filtered coordinates"),
    end: int = Query(..., description="End position in filtered coordinates"),
):
    """Get positions in a window of compressed coordinates."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]

    # Clamp to valid range
    max_idx = len(axis_table) - 1
    start = max(0, min(start, max_idx))
    end = max(start, min(end, max_idx))

    # Get window
    window_df = axis_table.filter(
        (pl.col('filtered_idx') >= start) &
        (pl.col('filtered_idx') <= end)
    )

    positions = window_df.to_dicts()

    # Compute real coordinate ranges
    real_ranges = []
    if len(window_df) > 0:
        real_ranges.append({
            "chrom": CHROMOSOME,
            "start": int(window_df['pos'].min()),
            "end": int(window_df['pos'].max())
        })

    return WindowResponse(
        filter_id=filter_id,
        window={
            "filtered_start": start,
            "filtered_end": end,
            "num_positions": len(window_df)
        },
        positions=positions,
        real_coordinate_ranges=real_ranges
    )


@app.get("/api/track-data")
async def get_track_data(
    track_id: str = Query(..., description="Track/column name"),
    filter_id: str = Query(..., description="Filter ID"),
    filtered_start: int = Query(..., ge=0),
    filtered_end: int = Query(...),
):
    """Get track values for a compressed coordinate window."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]
    bigwig_dir = get_bigwig_dir()

    # Get window positions
    window_df = axis_table.filter(
        (pl.col('filtered_idx') >= filtered_start) &
        (pl.col('filtered_idx') <= filtered_end)
    )

    # Check if track_id is in axis_table
    if track_id in window_df.columns:
        # For constraint columns, return full variant data
        if track_id in ['Constraint', 'Core', 'Complete']:
            values = []
            for row in window_df.select(['filtered_idx', track_id]).to_dicts():
                variant_data = extract_constraint_variants(row[track_id])
                values.append({
                    "filtered_idx": row["filtered_idx"],
                    "variants": [
                        {"allele": allele, "pred": sanitize_float(pred)}
                        for allele, pred in variant_data
                    ]
                })
            return {"track_id": track_id, "values": values}
        else:
            values = [
                {"filtered_idx": row["filtered_idx"], "value": sanitize_float(row[track_id])}
                for row in window_df.select(['filtered_idx', track_id]).to_dicts()
            ]
            return {"track_id": track_id, "values": values}

    # Try to load from BigWig
    bw_path = bigwig_dir / f'{track_id}.bw'
    if not bw_path.exists():
        raise HTTPException(status_code=404, detail=f"Track not found: {track_id}")

    try:
        bw = pyBigWig.open(str(bw_path))
        values = []
        for row in window_df.iter_rows(named=True):
            chrom = row['chrom']
            pos = row['pos']
            filtered_idx = row['filtered_idx']

            try:
                val = bw.values(chrom, pos, pos + 1)[0]
                sanitized_val = sanitize_float(val)
                values.append({
                    "filtered_idx": filtered_idx,
                    "value": sanitized_val if sanitized_val is not None else 0.0
                })
            except:
                values.append({"filtered_idx": filtered_idx, "value": 0.0})

        bw.close()
        return {"track_id": track_id, "values": values}

    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error reading BigWig: {str(e)}")


@app.get("/api/genes-in-window")
async def genes_in_window(
    filter_id: str = Query(..., description="Filter ID"),
    start: int = Query(..., description="Start position in compressed coordinates"),
    end: int = Query(..., description="End position in compressed coordinates"),
):
    """Get unique gene names visible in the compressed coordinate window."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]

    window = axis_table.filter(
        (pl.col('filtered_idx') >= start) &
        (pl.col('filtered_idx') <= end)
    )

    if len(window) == 0:
        return {"genes": [], "count": 0}

    genes = window.select('gene_symbol').unique().filter(
        pl.col('gene_symbol').is_not_null()
    ).sort('gene_symbol')

    gene_list = [row['gene_symbol'] for row in genes.to_dicts()]
    return {"genes": gene_list, "count": len(gene_list)}


@app.get("/api/axis-labels")
async def get_axis_labels(
    filter_id: str = Query(...),
    filtered_start: int = Query(...),
    filtered_end: int = Query(...),
    density: int = Query(default=10, description="Approximate number of labels"),
):
    """Get axis labels for compressed coordinate window."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]

    window_size = filtered_end - filtered_start
    step = max(1, window_size // density)

    labels = []
    for filtered_idx in range(filtered_start, filtered_end, step):
        result = axis_table.filter(pl.col('filtered_idx') == filtered_idx)
        if len(result) > 0:
            row = result.to_dicts()[0]
            labels.append({
                "filtered_idx": filtered_idx,
                "chrom": row['chrom'],
                "pos": row['pos'],
                "label": f"{row['chrom']}:{row['pos']:,}",
                "type": "position"
            })

    return {"labels": labels}


# =============================================================================
# Search Endpoints
# =============================================================================

@app.get("/api/search/gene", response_model=List[GeneSearchResult])
async def search_gene(
    gene: str = Query(..., description="Gene symbol (e.g., BRCA1)"),
    filter_id: str = Query(..., description="Filter ID"),
):
    """Search for a gene and get its compressed coordinates."""
    if filter_id not in gene_indexes:
        raise HTTPException(status_code=400, detail=f"Gene index not available for filter: {filter_id}")

    gene_index = gene_indexes[filter_id]

    results = gene_index.filter(
        pl.col('gene_symbol').str.to_uppercase() == gene.upper()
    )

    if len(results) == 0:
        return []

    return [GeneSearchResult(**row) for row in results.to_dicts()]


@app.get("/api/search/gene/autocomplete")
async def autocomplete_gene(
    query: str = Query(..., min_length=1, description="Partial gene symbol"),
    filter_id: str = Query(..., description="Filter ID"),
    limit: int = Query(10, ge=1, le=50, description="Maximum results"),
):
    """Autocomplete gene symbols with prefix matching."""
    if filter_id not in gene_indexes:
        raise HTTPException(status_code=400, detail=f"Gene index not available for filter: {filter_id}")

    gene_index = gene_indexes[filter_id]

    query_upper = query.upper()
    results = gene_index.filter(
        pl.col('gene_symbol').str.to_uppercase().str.starts_with(query_upper)
    ).head(limit)

    return [row['gene_symbol'] for row in results.to_dicts()]


@app.get("/api/search/position", response_model=PositionSearchResult)
async def search_position(
    chrom: str = Query(..., description="Chromosome (e.g., chr1)"),
    pos: int = Query(..., description="Position"),
    filter_id: str = Query(..., description="Filter ID"),
):
    """Convert real genomic position to compressed coordinate."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]

    result = axis_table.filter(
        (pl.col('chrom') == chrom) &
        (pl.col('pos') == pos)
    )

    if len(result) == 0:
        return PositionSearchResult(
            query={"chrom": chrom, "pos": pos},
            filter_id=filter_id,
            result=None
        )

    row = result.to_dicts()[0]
    return PositionSearchResult(
        query={"chrom": chrom, "pos": pos},
        filter_id=filter_id,
        result={
            "in_filter": True,
            "filtered_idx": row['filtered_idx'],
            "gene_symbol": row.get('gene_symbol'),
        }
    )


@app.get("/api/gene-aa-lookup")
async def gene_aa_lookup(
    filter_id: str = Query(..., description="Filter ID"),
    gene: str = Query(..., description="Gene symbol"),
    aa_start: int = Query(..., description="AA start position"),
    aa_end: int = Query(..., description="AA end position"),
):
    """Find genomic coordinates and compressed coordinates for a gene + AA range."""
    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Unknown filter: {filter_id}")

    axis_table = axis_tables[filter_id]

    if 'aa_pos' not in axis_table.columns:
        raise HTTPException(
            status_code=400,
            detail="AA position data not available in this dataset"
        )

    results = axis_table.filter(
        (pl.col('gene_symbol').str.to_uppercase() == gene.upper()) &
        (pl.col('aa_pos') >= aa_start) &
        (pl.col('aa_pos') <= aa_end)
    ).sort('aa_pos')

    if len(results) == 0:
        raise HTTPException(
            status_code=404,
            detail=f"No positions found for gene {gene} in AA range {aa_start}-{aa_end}"
        )

    first_row = results.row(0, named=True)
    last_row = results.row(len(results) - 1, named=True)

    return {
        "gene": gene,
        "aa_start": aa_start,
        "aa_end": aa_end,
        "compressed_start": first_row['filtered_idx'],
        "compressed_end": last_row['filtered_idx'],
        "genomic_start": first_row['pos'],
        "genomic_end": last_row['pos'],
        "chrom": first_row['chrom'],
        "positions_found": len(results)
    }


# =============================================================================
# Protein API Endpoints
# =============================================================================

@app.get("/api/protein/{gene}")
async def get_protein(
    gene: str,
    include_domains: bool = Query(True),
    include_structure: bool = Query(True),
):
    """Get protein information for a gene."""
    gene_upper = gene.upper()

    metadata = coord_mapper.get_structure_metadata(gene_upper)

    if not metadata:
        if not coord_mapper.has_gene(gene_upper):
            raise HTTPException(
                status_code=404,
                detail=f"Gene {gene} not found. Currently only SCN2A is supported."
            )
        return {
            "gene_symbol": gene_upper,
            "message": "Gene found but no structure metadata available"
        }

    response = {
        "gene_symbol": metadata.get('gene_symbol'),
        "gene_name": metadata.get('gene_name'),
        "uniprot_accession": metadata.get('uniprot_accession'),
        "transcript_id": metadata.get('transcript_id'),
        "protein_length": metadata.get('protein_length'),
        "chromosome": metadata.get('chromosome'),
        "strand": metadata.get('strand'),
    }

    if include_domains:
        response["domains"] = metadata.get('domains', [])

    if include_structure:
        response["structures"] = metadata.get('structures', {})

    return response


@app.get("/api/protein/{gene}/residues")
async def get_protein_residues(
    gene: str,
    start: int = Query(1, ge=1, description="Start residue (1-based)"),
    end: int = Query(None, description="End residue (defaults to protein length)"),
    filter_id: str = Query("mis_count_gt0", description="Filter ID for constraint data"),
    include_plddt: bool = Query(True, description="Include pLDDT scores"),
    include_constraints: bool = Query(True, description="Include constraint metrics"),
):
    """Get per-residue data for a protein region."""
    gene_upper = gene.upper()

    if not coord_mapper.has_gene(gene_upper):
        raise HTTPException(
            status_code=404,
            detail=f"Gene {gene} not found. Currently only SCN2A is supported."
        )

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    protein_length = metadata.get('protein_length', 2005) if metadata else 2005

    if end is None:
        end = protein_length

    start = max(1, min(start, protein_length))
    end = max(start, min(end, protein_length))

    positions = coord_mapper.get_protein_range(gene_upper, start, end)

    if not positions:
        raise HTTPException(
            status_code=404,
            detail=f"No positions found for {gene} residues {start}-{end}"
        )

    # Get constraint data from axis table
    if filter_id in axis_tables and include_constraints:
        axis_table = axis_tables[filter_id]
        gene_data = axis_table.filter(
            pl.col('gene_symbol').str.to_uppercase() == gene_upper
        )
        pos_to_data = {}
        for row in gene_data.to_dicts():
            pos_to_data[(row['chrom'], row['pos'])] = row
    else:
        pos_to_data = {}

    # Get pLDDT scores
    plddt_by_residue = {}
    if include_plddt and metadata:
        plddt_by_residue = metadata.get('plddt_by_residue', {})

    # Build residue-level response
    residues = []
    current_residue = None
    residue_data = None

    for pos_info in positions:
        residue_num = pos_info['protein_residue']

        if residue_num != current_residue:
            if residue_data is not None:
                residues.append(residue_data)

            current_residue = residue_num
            residue_data = {
                'residue': residue_num,
                'amino_acid': pos_info['ref_aa'],
                'genomic_positions': [],
                'plddt': plddt_by_residue.get(str(residue_num)) or plddt_by_residue.get(residue_num),
            }

        pos_key = (pos_info['chrom'], pos_info['pos'])
        genomic_info = {
            'chrom': pos_info['chrom'],
            'pos': pos_info['pos'],
            'codon_position': pos_info['codon_position'],
        }

        if pos_key in pos_to_data:
            constraint_data = pos_to_data[pos_key]
            genomic_info['filtered_idx'] = constraint_data.get('filtered_idx')
            genomic_info['clinvar_count'] = constraint_data.get('clinvar.clinvar_count', 0)
            genomic_info['clinvar_labels'] = constraint_data.get('clinvar.clinvar_label_list')
            genomic_info['alphamissense'] = sanitize_float(
                constraint_data.get('dbnsfp.max_AlphaMissense_am_pathogenicity')
            )
            genomic_info['mtr'] = sanitize_float(
                constraint_data.get('dbnsfp.max_RGC_MTR_MTR')
            )
            genomic_info['mis_count'] = constraint_data.get('rgc_mis_count')
            genomic_info['mis_oe_21bp'] = sanitize_float(
                constraint_data.get('rgc_mis_exomes_XX_XY_21bp_oe_af0epos00')
            )

        residue_data['genomic_positions'].append(genomic_info)

    if residue_data is not None:
        residues.append(residue_data)

    return {
        "gene_symbol": gene_upper,
        "range": {"start": start, "end": end},
        "protein_length": protein_length,
        "residue_count": len(residues),
        "residues": residues
    }


@app.get("/api/protein/{gene}/residue-scores")
async def get_residue_scores(
    gene: str,
    field: str = Query(..., description="Field ID to aggregate"),
    aggregation: str = Query("max", description="Aggregation method: max, min, or mean"),
    filter_id: str = Query("mis_count_gt0", description="Filter ID for constraint data"),
):
    """Get aggregated scores per residue for a specific field."""
    gene_upper = gene.upper()

    if not coord_mapper.has_gene(gene_upper):
        raise HTTPException(status_code=404, detail=f"Gene {gene} not found")

    if filter_id not in axis_tables:
        raise HTTPException(status_code=400, detail=f"Filter {filter_id} not available")

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    protein_length = metadata.get('protein_length', 2005) if metadata else 2005

    positions = coord_mapper.get_protein_range(gene_upper, 1, protein_length)

    if not positions:
        return {"gene_symbol": gene_upper, "field": field, "scores": {}, "range": [None, None]}

    axis_table = axis_tables[filter_id]
    gene_data = axis_table.filter(
        pl.col('gene_symbol').str.to_uppercase() == gene_upper
    )

    pos_to_data = {}
    for row in gene_data.to_dicts():
        pos_to_data[(row['chrom'], row['pos'])] = row

    is_constraint_stacked = field in CONSTRAINT_STACKED_FIELDS
    residue_values = {}

    for pos_info in positions:
        residue_num = pos_info['protein_residue']
        pos_key = (pos_info['chrom'], pos_info['pos'])

        if pos_key in pos_to_data:
            constraint_data = pos_to_data[pos_key]
            value = constraint_data.get(field)

            if value is None:
                continue

            if is_constraint_stacked:
                if isinstance(value, list):
                    for variant in value:
                        if isinstance(variant, dict):
                            pred = variant.get('_1')
                            if pred is not None and not (isinstance(pred, float) and (pred != pred)):
                                if residue_num not in residue_values:
                                    residue_values[residue_num] = []
                                residue_values[residue_num].append(float(pred))
            else:
                if not (isinstance(value, float) and (value != value)):
                    if residue_num not in residue_values:
                        residue_values[residue_num] = []
                    residue_values[residue_num].append(float(value))

    scores = {}
    all_values = []

    for residue, values in residue_values.items():
        if aggregation == 'max':
            score = max(values)
        elif aggregation == 'min':
            score = min(values)
        else:
            score = sum(values) / len(values)

        scores[residue] = round(score, 4)
        all_values.append(score)

    if is_constraint_stacked:
        value_range = [0.0, 1.0]
    else:
        value_range = [min(all_values), max(all_values)] if all_values else [None, None]

    return {
        "gene_symbol": gene_upper,
        "field": field,
        "aggregation": aggregation,
        "residue_count": len(scores),
        "range": value_range,
        "scores": scores,
        "is_constraint_stacked": is_constraint_stacked
    }


@app.get("/api/protein/{gene}/coordinate-map")
async def get_coordinate_map(
    gene: str,
    chrom: str = Query(None, description="Chromosome"),
    pos: int = Query(None, description="Genomic position"),
    residue: int = Query(None, description="Protein residue"),
):
    """Map between genomic and protein coordinates."""
    gene_upper = gene.upper()

    if not coord_mapper.has_gene(gene_upper):
        raise HTTPException(
            status_code=404,
            detail=f"Gene {gene} not found. Currently only SCN2A is supported."
        )

    if chrom and pos:
        result = coord_mapper.genomic_to_protein(gene_upper, chrom, pos)
        if result is None:
            raise HTTPException(
                status_code=404,
                detail=f"Position {chrom}:{pos} not found in {gene} CDS"
            )
        return {
            "query": {"chrom": chrom, "pos": pos},
            "gene_symbol": gene_upper,
            "mapping": result
        }

    elif residue:
        results = coord_mapper.protein_to_genomic(gene_upper, residue)
        if not results:
            raise HTTPException(
                status_code=404,
                detail=f"Residue {residue} not found in {gene}"
            )
        return {
            "query": {"residue": residue},
            "gene_symbol": gene_upper,
            "genomic_positions": results
        }

    else:
        raise HTTPException(
            status_code=400,
            detail="Provide either (chrom, pos) or residue parameter"
        )


# =============================================================================
# Structure API Endpoints
# =============================================================================

@app.get("/api/structure/{gene}")
async def get_structure(gene: str):
    """Get structure information for a gene."""
    gene_upper = gene.upper()
    data_dir = get_data_dir()

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    if not metadata:
        raise HTTPException(
            status_code=404,
            detail=f"No structure data available for {gene}. Currently only SCN2A is supported."
        )

    structures = metadata.get('structures', {})

    response = {
        "gene_symbol": gene_upper,
        "uniprot_accession": metadata.get('uniprot_accession'),
        "protein_length": metadata.get('protein_length'),
        "available_structures": list(structures.keys()),
    }

    if 'alphafold' in structures:
        af = structures['alphafold']
        response['alphafold'] = {
            'id': af.get('id'),
            'version': af.get('version'),
            'coverage': af.get('coverage'),
            'mean_plddt': af.get('mean_plddt'),
            'min_plddt': af.get('min_plddt'),
            'max_plddt': af.get('max_plddt'),
            'confidence_distribution': af.get('confidence_distribution'),
            'file_url': f"/api/structure/{gene}/file/alphafold"
        }

    response['domains'] = metadata.get('domains', [])

    return response


@app.get("/api/structure/{gene}/file/{structure_type}")
async def get_structure_file(gene: str, structure_type: str):
    """Get the actual structure file (PDB format) for a gene."""
    gene_upper = gene.upper()
    data_dir = get_data_dir()

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    if not metadata:
        raise HTTPException(
            status_code=404,
            detail=f"No structure data available for {gene}"
        )

    structures = metadata.get('structures', {})

    if structure_type not in structures:
        raise HTTPException(
            status_code=404,
            detail=f"Structure type {structure_type} not available for {gene}"
        )

    structure_info = structures[structure_type]
    file_path = data_dir / structure_info.get('file', '')

    if not file_path.exists():
        raise HTTPException(
            status_code=404,
            detail=f"Structure file not found: {file_path}"
        )

    return FileResponse(
        path=file_path,
        media_type="chemical/x-pdb",
        filename=file_path.name
    )


@app.get("/api/structure/{gene}/plddt")
async def get_plddt_scores(
    gene: str,
    start: int = Query(1, ge=1),
    end: int = Query(None),
):
    """Get pLDDT scores for a range of residues."""
    gene_upper = gene.upper()

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    if not metadata:
        raise HTTPException(
            status_code=404,
            detail=f"No structure data available for {gene}"
        )

    plddt_by_residue = metadata.get('plddt_by_residue', {})
    protein_length = metadata.get('protein_length', len(plddt_by_residue))

    if end is None:
        end = protein_length

    scores = []
    for residue in range(start, end + 1):
        score = plddt_by_residue.get(str(residue)) or plddt_by_residue.get(residue)
        if score is not None:
            scores.append({'residue': residue, 'plddt': score})

    return {
        "gene_symbol": gene_upper,
        "range": {"start": start, "end": end},
        "scores": scores
    }


@app.get("/api/structure/{gene}/residue-data")
async def get_structure_residue_data(
    gene: str,
    residues: str = Query(..., description="Comma-separated list of residue numbers"),
    filter_id: str = Query("mis_count_gt0"),
):
    """Get combined structure and constraint data for specific residues."""
    gene_upper = gene.upper()

    if not coord_mapper.has_gene(gene_upper):
        raise HTTPException(status_code=404, detail=f"Gene {gene} not found")

    try:
        residue_list = [int(r.strip()) for r in residues.split(',')]
    except ValueError:
        raise HTTPException(
            status_code=400,
            detail="Invalid residue format. Provide comma-separated integers."
        )

    metadata = coord_mapper.get_structure_metadata(gene_upper)
    plddt_by_residue = metadata.get('plddt_by_residue', {}) if metadata else {}

    axis_table = axis_tables.get(filter_id)
    gene_data = {}
    if axis_table is not None:
        gene_rows = axis_table.filter(
            pl.col('gene_symbol').str.to_uppercase() == gene_upper
        )
        for row in gene_rows.to_dicts():
            gene_data[(row['chrom'], row['pos'])] = row

    residue_data = []
    for residue in residue_list:
        positions = coord_mapper.protein_to_genomic(gene_upper, residue)
        if not positions:
            continue

        clinvar_count = 0
        clinvar_labels = []
        max_alphamissense = None
        filtered_indices = []

        for pos in positions:
            key = (pos['chrom'], pos['pos'])
            if key in gene_data:
                data = gene_data[key]
                clinvar_count += data.get('clinvar.clinvar_count', 0) or 0
                labels = data.get('clinvar.clinvar_label_list')
                if labels:
                    if isinstance(labels, list):
                        clinvar_labels.extend(labels)
                    else:
                        clinvar_labels.append(labels)

                am = data.get('dbnsfp.max_AlphaMissense_am_pathogenicity')
                if am is not None and not math.isnan(am):
                    if max_alphamissense is None or am > max_alphamissense:
                        max_alphamissense = am

                if data.get('filtered_idx') is not None:
                    filtered_indices.append(data['filtered_idx'])

        residue_data.append({
            'residue': residue,
            'amino_acid': positions[0]['ref_aa'] if positions else None,
            'plddt': plddt_by_residue.get(str(residue)) or plddt_by_residue.get(residue),
            'clinvar_count': clinvar_count,
            'clinvar_labels': list(set(clinvar_labels)),
            'max_alphamissense': sanitize_float(max_alphamissense),
            'filtered_indices': filtered_indices,
            'genomic_positions': positions
        })

    return {
        "gene_symbol": gene_upper,
        "residue_count": len(residue_data),
        "residues": residue_data
    }


# =============================================================================
# Static Files & Frontend
# =============================================================================

# Get project root for static file mounting
project_root = get_project_root()
frontend_dist = project_root / "browser" / "frontend" / "dist"
frontend_dev = project_root / "browser" / "frontend"

# Mount frontend - prefer dist (built) over dev
if frontend_dist.exists():
    app.mount("/", StaticFiles(directory=str(frontend_dist), html=True), name="frontend")
elif frontend_dev.exists() and (frontend_dev / "index.html").exists():
    app.mount("/", StaticFiles(directory=str(frontend_dev), html=True), name="frontend")


# =============================================================================
# Main Entry Point
# =============================================================================

if __name__ == "__main__":
    import uvicorn
    config = get_config()
    host = config['browser'].get('host', '0.0.0.0')
    port = config['browser'].get('port', 8000)
    uvicorn.run(app, host=host, port=port, log_level="info")
