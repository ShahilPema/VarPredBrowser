"""
Pydantic models for VarPredBrowser API

Defines request/response models for the FastAPI endpoints.
"""

from pydantic import BaseModel
from typing import List, Optional, Dict, Any


class FilterInfo(BaseModel):
    """Information about a data filter."""
    id: str
    name: str
    description: str
    total_positions: int
    chromosome: str


class Position(BaseModel):
    """A single genomic position."""
    filtered_idx: int
    chrom: str
    pos: int
    gene_symbol: Optional[str] = None
    any_count: Optional[int] = None
    mis_count: Optional[int] = None


class WindowResponse(BaseModel):
    """Response for a window query in compressed coordinates."""
    filter_id: str
    window: Dict[str, Any]
    positions: List[Dict[str, Any]]
    real_coordinate_ranges: List[Dict[str, Any]]


class GeneSearchResult(BaseModel):
    """Result of a gene search."""
    gene_symbol: str
    chrom: str
    pos_start: int
    pos_end: int
    filtered_idx_start: int
    filtered_idx_end: int
    num_positions: int


class PositionSearchResult(BaseModel):
    """Result of a position search."""
    query: Dict[str, Any]
    filter_id: str
    result: Optional[Dict[str, Any]]


class VariantPrediction(BaseModel):
    """A single variant prediction."""
    allele: str
    pred: Optional[float]


class TrackValue(BaseModel):
    """A track value at a filtered position."""
    filtered_idx: int
    value: Optional[float] = None
    variants: Optional[List[VariantPrediction]] = None


class TrackDataResponse(BaseModel):
    """Response for track data query."""
    track_id: str
    values: List[TrackValue]


class ResidueData(BaseModel):
    """Data for a single protein residue."""
    residue: int
    amino_acid: Optional[str]
    plddt: Optional[float]
    genomic_positions: List[Dict[str, Any]]


class ProteinResiduesResponse(BaseModel):
    """Response for protein residues query."""
    gene_symbol: str
    range: Dict[str, int]
    protein_length: int
    residue_count: int
    residues: List[Dict[str, Any]]


class ResidueScoresResponse(BaseModel):
    """Response for residue scores query."""
    gene_symbol: str
    field: str
    aggregation: str
    residue_count: int
    range: List[Optional[float]]
    scores: Dict[int, float]
    is_constraint_stacked: bool = False


class CoordinateMappingResponse(BaseModel):
    """Response for coordinate mapping query."""
    query: Dict[str, Any]
    gene_symbol: str
    mapping: Optional[Dict[str, Any]] = None
    genomic_positions: Optional[List[Dict[str, Any]]] = None


class StructureInfo(BaseModel):
    """Information about a protein structure."""
    gene_symbol: str
    uniprot_accession: Optional[str]
    protein_length: Optional[int]
    available_structures: List[str]
    alphafold: Optional[Dict[str, Any]] = None
    domains: List[Dict[str, Any]] = []


class PlddtScore(BaseModel):
    """pLDDT score for a residue."""
    residue: int
    plddt: float


class PlddtResponse(BaseModel):
    """Response for pLDDT scores query."""
    gene_symbol: str
    range: Dict[str, int]
    scores: List[PlddtScore]
