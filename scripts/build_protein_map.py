#!/usr/bin/env python3
"""
Build genomic-to-protein coordinate mapping for a gene.

This script creates a mapping table that links each genomic position
to its corresponding protein residue for the MANE Select transcript.

Example for SCN2A:
- Gene: SCN2A (sodium voltage-gated channel alpha subunit 2)
- MANE Select transcript: ENST00000375437
- UniProt: Q99250
- Protein length: 2005 amino acids

Usage:
    python scripts/build_protein_map.py [--gene GENE] [--transcript TRANSCRIPT] [--uniprot UNIPROT]
"""

import argparse
import polars as pl
import requests
from pathlib import Path
from typing import List, Dict, Tuple

# Try to import from browser config, fall back to local path
try:
    import sys
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from browser.backend.config import get_data_dir
    DATA_DIR = get_data_dir()
except ImportError:
    DATA_DIR = Path(__file__).parent.parent / 'data'

# Default configuration (SCN2A)
DEFAULT_GENE_SYMBOL = 'SCN2A'
DEFAULT_TRANSCRIPT_ID = 'ENST00000375437'
DEFAULT_UNIPROT_ID = 'Q99250'

# Genetic code for codon to amino acid translation
CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def fetch_ensembl_exons(transcript_id: str) -> Tuple[Dict, List[Dict]]:
    """
    Fetch exon information from Ensembl REST API for a transcript.
    Returns list of exons with their genomic coordinates.
    """
    # Remove version if present
    transcript_base = transcript_id.split('.')[0]

    url = f"https://rest.ensembl.org/lookup/id/{transcript_base}?expand=1"
    headers = {"Content-Type": "application/json"}

    print(f"Fetching transcript info from Ensembl: {transcript_base}")
    response = requests.get(url, headers=headers)
    response.raise_for_status()

    data = response.json()

    # Get transcript details
    print(f"  Gene: {data.get('display_name', 'N/A')}")
    print(f"  Biotype: {data.get('biotype', 'N/A')}")
    print(f"  Strand: {'+' if data.get('strand', 1) == 1 else '-'}")

    # Fetch exons
    exon_url = f"https://rest.ensembl.org/overlap/id/{transcript_base}?feature=exon"
    exon_response = requests.get(exon_url, headers=headers)
    exon_response.raise_for_status()

    exons = [e for e in exon_response.json() if e.get('Parent') == transcript_base]
    exons.sort(key=lambda x: x['start'])

    return data, exons


def fetch_ensembl_cds(transcript_id: str) -> Tuple[str, List[Dict]]:
    """
    Fetch CDS (coding sequence) information from Ensembl.
    Returns the CDS sequence and CDS region coordinates.
    """
    transcript_base = transcript_id.split('.')[0]

    # Get CDS sequence
    seq_url = f"https://rest.ensembl.org/sequence/id/{transcript_base}?type=cds"
    headers = {"Content-Type": "text/plain"}

    print(f"Fetching CDS sequence...")
    seq_response = requests.get(seq_url, headers=headers)
    seq_response.raise_for_status()
    cds_sequence = seq_response.text.strip()

    print(f"  CDS length: {len(cds_sequence)} bp")
    print(f"  Protein length: {len(cds_sequence) // 3} aa")

    # Get CDS regions (genomic coordinates)
    headers = {"Content-Type": "application/json"}
    cds_url = f"https://rest.ensembl.org/overlap/id/{transcript_base}?feature=cds"
    cds_response = requests.get(cds_url, headers=headers)
    cds_response.raise_for_status()

    cds_regions = [c for c in cds_response.json() if c.get('Parent') == transcript_base]
    cds_regions.sort(key=lambda x: x['start'])

    print(f"  CDS regions: {len(cds_regions)}")

    return cds_sequence, cds_regions


def build_genomic_to_protein_map(
    cds_sequence: str,
    cds_regions: List[Dict],
    strand: int,
    chrom: str
) -> pl.DataFrame:
    """
    Build a mapping from genomic positions to protein residues.

    For each genomic position in the CDS:
    - Map to CDS offset (0-based position in coding sequence)
    - Calculate codon position (1, 2, or 3)
    - Calculate protein residue (1-based)
    - Determine reference amino acid
    """

    mappings = []

    # Build list of all CDS genomic positions in transcript order
    cds_positions = []
    for region in cds_regions:
        start = region['start']
        end = region['end']
        positions = list(range(start, end + 1))
        if strand == -1:
            positions = positions[::-1]  # Reverse for minus strand
        cds_positions.extend([(pos, region['seq_region_name']) for pos in positions])

    # For minus strand, we need to reverse the entire list
    if strand == -1:
        cds_positions = cds_positions[::-1]

    print(f"  Total CDS positions: {len(cds_positions)}")

    # Map each position
    for cds_offset, (genomic_pos, seq_region) in enumerate(cds_positions):
        codon_position = (cds_offset % 3) + 1  # 1, 2, or 3
        protein_residue = (cds_offset // 3) + 1  # 1-based

        # Get the codon for this position
        codon_start = (cds_offset // 3) * 3
        if codon_start + 3 <= len(cds_sequence):
            codon = cds_sequence[codon_start:codon_start + 3]
            ref_aa = CODON_TABLE.get(codon.upper(), 'X')
        else:
            ref_aa = 'X'

        mappings.append({
            'chrom': f'chr{seq_region}' if not seq_region.startswith('chr') else seq_region,
            'pos': genomic_pos,
            'cds_offset': cds_offset,
            'codon_position': codon_position,
            'protein_residue': protein_residue,
            'ref_aa': ref_aa,
            'strand': '+' if strand == 1 else '-'
        })

    df = pl.DataFrame(mappings)
    return df


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    """Fetch protein sequence from UniProt."""
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta"
    response = requests.get(url)
    response.raise_for_status()

    # Parse FASTA
    lines = response.text.strip().split('\n')
    sequence = ''.join(lines[1:])
    return sequence


def main(gene_symbol: str, transcript_id: str, uniprot_id: str):
    print("=" * 70)
    print(f"Building Genomic-to-Protein Coordinate Map for {gene_symbol}")
    print("=" * 70)

    # Create data directory if needed
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # Fetch transcript info from Ensembl
    print("\n[Step 1] Fetching transcript information from Ensembl...")
    transcript_info, exons = fetch_ensembl_exons(transcript_id)
    strand = transcript_info.get('strand', 1)
    chrom = transcript_info.get('seq_region_name', '2')

    # Fetch CDS
    print("\n[Step 2] Fetching CDS sequence and regions...")
    cds_sequence, cds_regions = fetch_ensembl_cds(transcript_id)

    # Verify against UniProt
    print("\n[Step 3] Verifying against UniProt sequence...")
    uniprot_seq = fetch_uniprot_sequence(uniprot_id)
    print(f"  UniProt sequence length: {len(uniprot_seq)} aa")

    # Translate CDS and compare
    translated = ''
    for i in range(0, len(cds_sequence) - 2, 3):
        codon = cds_sequence[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        translated += aa

    print(f"  Translated CDS length: {len(translated)} aa")

    if translated == uniprot_seq:
        print("  Sequences match!")
    else:
        # Check for minor differences
        matches = sum(1 for a, b in zip(translated, uniprot_seq) if a == b)
        print(f"  Warning: Sequences differ: {matches}/{min(len(translated), len(uniprot_seq))} positions match")

    # Build mapping
    print("\n[Step 4] Building genomic-to-protein mapping...")
    mapping_df = build_genomic_to_protein_map(cds_sequence, cds_regions, strand, chrom)

    # Add gene symbol and transcript info
    mapping_df = mapping_df.with_columns([
        pl.lit(gene_symbol).alias('gene_symbol'),
        pl.lit(transcript_id).alias('transcript_id'),
        pl.lit(uniprot_id).alias('uniprot_accession')
    ])

    # Save
    output_path = DATA_DIR / f'{gene_symbol.lower()}_protein_map.parquet'
    mapping_df.write_parquet(output_path)
    print(f"\nSaved mapping to: {output_path}")
    print(f"  Total positions: {len(mapping_df)}")
    print(f"  Protein residues: {mapping_df['protein_residue'].max()}")

    # Show sample
    print("\nSample mapping (first 10 positions):")
    print(mapping_df.head(10))

    print(f"\nSample mapping (around residue 100):")
    print(mapping_df.filter(pl.col('protein_residue') == 100))

    return mapping_df


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Build genomic-to-protein coordinate mapping')
    parser.add_argument('--gene', default=DEFAULT_GENE_SYMBOL, help='Gene symbol')
    parser.add_argument('--transcript', default=DEFAULT_TRANSCRIPT_ID, help='Ensembl transcript ID')
    parser.add_argument('--uniprot', default=DEFAULT_UNIPROT_ID, help='UniProt accession')

    args = parser.parse_args()
    main(args.gene, args.transcript, args.uniprot)
